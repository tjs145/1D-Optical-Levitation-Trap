from trap_model import Trapped_particle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
os.chdir('various radii')


file='default_2.txt' 
s = Trapped_particle(1e2, 0.4, parameter_filename=file)  
file2='default_4.txt' 
s2 = Trapped_particle(1e2, 0.4, parameter_filename=file2)  
file3='default_6.txt' 
s3 = Trapped_particle(1e2, 0.4, parameter_filename=file3)  

plt.figure(1) #Fig 1 shows position with time
space_ax = plt.subplot(1,1,1)

plt.figure(2) #Fig 2 shows power spectrum
amplitude_ax = plt.subplot(1,1,1)

s.run_ivp(1000,0,3000,10000)
s2.run_ivp(1000,0,3000,10000)
s3.run_ivp(1000,0,3000,10000)

equilibrium=np.mean(s.positions[9000:])
equilibrium2=np.mean(s2.positions[9000:])
equilibrium3=np.mean(s3.positions[9000:])

#Define various laser power functions of varying frequencies 

lasers=[]

class laser:
    def __init__(self, f):
        self.f = f
    def power(self, t):
        return 0.4 + 0.01*np.sin(2*np.pi*self.f*t/1000)

frequencies = np.logspace(-1.5,2,num=20,base=10)

for f in frequencies:
    lasers.append(laser(f))

#Create simulation instances for each laser power function
simulations = []
simulations2 = []
simulations3 = []

for l in lasers:
    simulations.append(Trapped_particle(1e5, l.power, parameter_filename=file))
    simulations2.append(Trapped_particle(1e5, l.power, parameter_filename=file2))
    simulations3.append(Trapped_particle(1e5, l.power, parameter_filename=file3))

amplitudes=[]
amplitudes2=[]
amplitudes3=[]

for s in simulations:
    s.run_ivp(equilibrium, 0, 30000,50000)
    
    amplitudes.append(np.ptp(s.positions)/2)

    space_ax.plot(s.times, s.positions)

amplitude_ax.plot(frequencies, amplitudes, 'rx', label='2$\mu$m')

for s in simulations2:
    s.run_ivp(equilibrium2, 0, 30000,50000)
    
    amplitudes2.append(np.ptp(s.positions)/2)

    space_ax.plot(s.times, s.positions)

amplitude_ax.plot(frequencies, amplitudes2, 'gx', label='4$\mu$m')

for s in simulations3:
    s.run_ivp(equilibrium3, 0, 30000,50000)
    
    amplitudes3.append(np.ptp(s.positions)/2)

    space_ax.plot(s.times, s.positions)

amplitude_ax.plot(frequencies, amplitudes3, 'bx', label='6$\mu$m')


########### Fit response function to data and plot

def A(f, fc, damp, a0):
    return a0*(fc/np.sqrt((fc-(2*np.pi/damp)*f**2)**2+f**2))
    
plot_frequencies = np.logspace(-1.5,2,num=1000,base=10)

popt, pcov = curve_fit(A, frequencies, amplitudes, p0=[0.3, 3412, 23])

print(popt)

amplitude_ax.plot(plot_frequencies, A(plot_frequencies, popt[0], popt[1], popt[2]), 'r-', label='2$\mu$m')

popt, pcov = curve_fit(A, frequencies, amplitudes2, p0=[0.3, 3412, 23])

print(popt)

amplitude_ax.plot(plot_frequencies, A(plot_frequencies, popt[0], popt[1], popt[2]), 'g-', label='4$\mu$m')

popt, pcov = curve_fit(A, frequencies, amplitudes3, p0=[0.3, 3412, 23])

print(popt)

amplitude_ax.plot(plot_frequencies, A(plot_frequencies, popt[0], popt[1], popt[2]), 'b-', label='6$\mu$m')

######

amplitude_ax.set_xscale('log')
amplitude_ax.set_xlabel('Laser modulation frequency (Hz)')
amplitude_ax.set_ylabel('Amplitude ($\mu$m)')
amplitude_ax.legend()
plt.show()
