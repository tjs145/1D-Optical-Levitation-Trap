import time as t
import matplotlib.pyplot as plt
import spidev

total_samples=1000000
count=0
voltages=[]
times=[]
vref=3.3 #ref voltage always 3.3v for raspberry pi

adc = spidev.SpiDev()
adc.open(0,0)
adc.max_speed_hz = 8000000 #Any faster than this causes glitches. 

start_time=t.time()

while True:#count<total_samples:
    out = adc.xfer2([0b1, 0b10, 0]) #First byte contains start bit, second contains mode (1) and channel (0), third byte added because only receives number of bytes sent and needs to receive 3 bytes.
    raw = ((out[1] & 0x0f)<<8) + out[2] #combines the 12 data bits (last 4 of second byte and all 8 of 3rd byte)
    voltage = vref*(raw/4095)
    voltages.append(voltage)
    times.append(t.time())
    count+=1
    
end_time=t.time()

total_time = end_time - start_time

sample_rate = total_samples/total_time

print(sample_rate)

for i in range(len(times)):
    times[i]-=start_time

plt.plot(times[100:200], voltages[100:200])
plt.plot(times[100:200], voltages[100:200], 'rx')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.show()

