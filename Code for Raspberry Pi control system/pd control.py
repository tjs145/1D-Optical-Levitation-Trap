import time as t
import spidev
import numpy as np
import matplotlib.pyplot as plt
import os

adc = spidev.SpiDev()
adc.open(0,0)
adc.max_speed_hz = 8000000

dac = spidev.SpiDev()
dac.open(0,1)
dac.max_speed_hz = 50000000


#Settings:
wait_time = 10
control_time = 20
P = 50000
D = 600

channel1_raw = []
channel2_raw = []
times=[]
outputs = []

#Set the DAC output to the median voltage readble to the Pockels cell driver
out=1650
dac.writebytes([(0b0001 << 4) + (out>>8), 0b11111111 & out])

t.wait(1) #Plenty of time to let the particle reach equilibrium

start = t.time()

while t.time() - start < wait_time:
    #During wait time, channel data is recorded
    out = adc.xfer2([1, 0b10000000, 0])
    ch1 = ((out[1] & 0x0f)<<8) + out[2]
    out = adc.xfer2([1, 0b11000000, 0])
    ch2 = ((out[1] & 0x0f)<<8) + out[2]
    channel1_raw.append(ch1)
    channel2_raw.append(ch2)
    times.append(t.time())
    outputs.append(0)

#Calculate mean position from wait time so it can be sued as the reference position.
ch1=np.mean(channel1_raw) 
ch2=np.mean(channel2_raw)
ref = ((1899-ch1-1849+ch2)/(1.06*(1899+1849-ch1-ch2))) #Reference position. Ranges from -1 to 1 along length of PSD active region


while t.time()-start < sample_time:
    out = adc.xfer2([1, 0b10000000, 0])
    ch1 = ((out[1] & 0x0f)<<8) + out[2]
    out = adc.xfer2([1, 0b11000000, 0])
    ch2 = ((out[1] & 0x0f)<<8) + out[2]
    channel1_raw.append(ch1)
    channel2_raw.append(ch2)
    ch1 = np.mean(channel1_raw[-10:])
    ch2 = np.mean(channel2_raw[-10:])
    relative_position = ((1899-ch1-1849+ch2)/(1.06*(1899+1849-ch1-ch2)))-ref
    pre_ch1 = np.mean(channel1_raw[-20:-10])
    pre_ch2 = np.mean(channel2_raw[-20:-10])
    dt = np.mean(times[-10:])-np.mean(times[-20:-10])
    previous_position = ((1899-pre_ch1-1849+pre_ch2)/(1.06*(1899+1849-pre_ch1-pre_ch2)))-ref
    derivative = (relative_position-previous_position)/dt
    out = int(1650-P*relative_position-D*derivative) 
    if out > 4095:
        out=4095   #Doesn't try to output more than the maximum value
    elif out < 0:
        out=0 #Doesn't try to output a negative voltage
    dac.writebytes([(0b0001 << 4) + (out>>8), 0b11111111 & out])
    outputs.append(out)
    times.append(t.time())

channel1 = (3.3/4095)*np.array(channel1_raw)
channel2 = (3.3/4095)*np.array(channel2_raw)

outputs = (3.3/4095)*np.array(outputs)

times = np.array(times)-start

#save arrays to files
np.savetxt('channel1.txt', channel1)
np.savetxt('channel2.txt', channel2)
np.savetxt('output.txt', outputs)
np.savetxt('times.txt', times)

position = -1*(1.75/1.06)*((1.53-channel1-1.49+channel2)/(1.53+1.49-channel1-channel2))


#Calculate moving averages of data to make plots less noisy
s = 10
n=len(position)-s

av_p = np.zeros((n))
av_o = np.zeros((n))
av_t = np.zeros((n))
for i in range(n):
    av_p[i] = np.mean(position[i:i+s])
    av_o[i] = np.mean(outputs[i:i+s])
    av_t[i] = np.mean(times[i:i+s])

#Plot averaged position data and DAC output data    

plt.plot(av_t, av_p)
plt.figure(2)
plt.plot(av_t, av_o, 'r')

plt.show()

