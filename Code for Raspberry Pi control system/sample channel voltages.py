import time as t
import matplotlib.pyplot as plt
import spidev
import numpy as np

total_samples=1000
count=0
voltages=[] #Stores voltages from channel 1
voltages2 = [] #Stores voltages from channel 2
times=[] #Stores times at which samples are taken
vref=3.3 #ref voltage always 3.3v for raspberry pi

adc = spidev.SpiDev()
adc.open(0,0) #Opens the SPI interface
adc.max_speed_hz = 8000000   #Highest reliable SPI clock rate (9MHz caused glitches)

start_time=t.time()

while count<total_samples:
    out = adc.xfer2([1, 0b10000000, 0]) #First byte contains start bit, second contains mode (1) and channel (0), third byte added because only receives number of bytes sent and needs to receive 3 bytes.
    raw = ((out[1] & 0x0f)<<8) + out[2] #combines the 12 data bits (last 4 of second byte and all 8 of 3rd byte)
    voltage = vref*(raw/4095)
    voltages.append(voltage)
    
    out = adc.xfer2([1, 0b11000000, 0])#First byte contains start bit, second contains mode (1) and channel (0), third byte added because only receives number of bytes sent and needs to receive 3 bytes.
    raw = ((out[1] & 0x0f)<<8) + out[2] #combines the 12 data bits (last 4 of second byte and all 8 of 3rd byte)
    voltage = vref*(raw/4095)
    voltages2.append(voltage)
    
    times.append(t.time())
    count+=1
    
end_time=t.time()

total_time = end_time - start_time

sample_rate = total_samples/total_time

print(sample_rate)
print(np.mean(voltages)*3) #Multiplying by 3 accounts for the 1/3 voltage divider

for i in range(len(times)):
    times[i]-=start_time   #Shifts the time values so they start at 0

plt.plot(times, voltages, 'b')
plt.plot(times, voltages, 'bx')
plt.plot(times, voltages2, 'r')
plt.plot(times, voltages2, 'rx')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.show()

