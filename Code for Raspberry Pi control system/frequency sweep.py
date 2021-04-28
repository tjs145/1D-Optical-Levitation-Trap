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
sample_time = 120

def frequency(t):
    #Increases every 10s (except 0.1Hz which lasts for 20 so that its properly sampled)
    if t>110:
        return 1000
    elif t>100:
        return 500
    elif t>90:
        return 200
    elif t>80:
        return 120
    elif t>70:
        return 70
    elif t>60:
        return 30
    elif t>50:
        return 10
    elif t>40:
        return 5
    elif t>30:
        return 1
    elif t>20:
        return 0.5
    else:
        return 0.1



channel1_raw = [] #Store raw data from channels
channel2_raw = []
times = [] #Stores time data
outputs=[] #Stores DAC outputs

def angular_frequency(t):
    return 2*np.pi*frequency(t)

start = t.time()

while t.time()-start < sample_time:
    out = adc.xfer2([1, 0b10000000, 0])#First byte contains start bit, second contains mode (1) and channel (0), third byte added because only receives number of bytes sent and needs to receive 3 bytes.
    raw = ((out[1] & 0x0f)<<8) + out[2] #Raw output from channel 1
    channel1_raw.append(raw)
    out = adc.xfer2([1, 0b11000000, 0]) #First byte contains start bit, second contains mode (1) and channel (1), third byte added because only receives number of bytes sent and needs to receive 3 bytes.
    raw = ((out[1] & 0x0f)<<8) + out[2] #Raw output from channel 2
    channel2_raw.append(raw)
    out = int(1000*(1.65+1.65*np.sin((t.time()-start)*angular_frequency(t.time()-start)))) #Calculates raw output to be sent to DAC
    dac.writebytes([(0b0001 << 4) + (out>>8), 0b11111111 & out]) #Sends data to DAC
    times.append(t.time()) 
    outputs.append(out)
    

os.chdir('frequency sweep')
np.savetxt('channel 1.txt', channel1_raw, delimiter=',')
np.savetxt('channel 2.txt', channel2_raw, delimiter=',')
np.savetxt('Times.txt', times, delimiter=',')


channel1 = (3.3/4095)*np.array(channel1_raw)
channel2 = (3.3/4095)*np.array(channel2_raw)
position = (1.75/1.06)*((1.53-channel1-1.49+channel2)/(1.53+1.49-channel1-channel2))
times = np.array(times) - start

plt.plot(times, outputs)
plt.figure(2)
plt.plot(times, position)

plt.show()
