# Code for Raspberry Pi control system

These codes are designed for controlling the ADCDAC Pi Zero from AB Electronics.

'adc sample rate test.py' takes 100000 voltage readings from one of the ADC channels as quickly as possible, measures the time taken, and outputs the sample  rate. It also plots a subset of the samples taken so any glitchy behaviour can be identified. 

'dac sample rate test.py' does the same thing, but instead outputs the voltages through the DAC rather than reading through the ADC. It also doesn't produce a plot, but the outputs can be investigated seperately using an oscilloscope. 

'frequency sweep.py' uses the DAC to output sine waves of discretely increasing frequencies, which can be connected to the Pockels cell driver and used to drive oscillations in the trapped particle. It also uses the ADC to read voltages from the position detector, which can then be used to calculate the position of the particle at a given time. The particles amplitude of motion at different frequencies can theoretically be used to make estimations of it's size.

'pd control.py' implements a PD control algorithm. For the first 10 seconds, position data is read by the ADC, but no voltages are applied to the Pockels cell driver by the DAC. The particles average position during this time is calculated and used as the reference position for the PD algorithm. The PD algorithm involves using the DAC to apply a voltage across the Pockels cell driver which is a weighted sum of the trapped particles displacement from the reference position, and the derivative of it's motion (it's velocity).

'sample channel voltages.py' simply uses the ADC to take 1000 voltage readings from each channel, and plots the readings. This is useful for diagnosing noise in the system. 

'timed sample.py' is similar. It takes samples from both channels over 5 seconds, plots them, and saves the data to three separate .txt files (one file contains channel 1 data, one for channel 2 data, and one contains the times at which the samples were taken.)
