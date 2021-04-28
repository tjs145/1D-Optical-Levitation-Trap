import spidev
import time as t

dac = spidev.SpiDev()
dac.open(0,1)
dac.max_speed_hz = 50000000

bytes_high=[0b00011111, 0b11111111] #4 configuration bits (DAC A or B?, Ignored, Gain selection, Shutdown control) plus 12 data bits (111111111111 for max voltage)
bytes_low=[0b00010000, 0b00000000]

total_outputs=100000

start=t.time()

for i in range(total_outputs):
    dac.writebytes(bytes_high)
    t.sleep(0.005)
    dac.writebytes(bytes_low)
    t.sleep(0.005)
    
end=t.time()

total_time=end-start
sample_rate = 2*total_outputs/total_time

print(sample_rate)
    
dac.close()