import os

os.chdir('various radii')  

files = os.listdir()

from trap_model import Trapped_particle
import matplotlib.pyplot as plt
import numpy as np

initial_laser_power = 0.3  #Change to relevant powers
decreased_laser_power = 0.15

def lp(t):  
    if t<1000:
        return initial_laser_power
    else:
        return decreased_laser_power

simulations=[]

for file in files:
    simulations.append(Trapped_particle(1, initial_laser_power, parameter_filename=file))  #Initialise low pressure model so the (approximate) equilibrium position can be quickly found

plt.figure(1) #Fig 1 shows positions with time

space_ax = plt.subplot(1,1,1)

for s in simulations:
    s.run_ivp(1000,0,6000,10000)


#Find equilibrium point for initial laser power and re-run at 1 atm with changing laser power function

    
plt.figure(2) #Fig 2 plots equilibrium change for various radii    
eq_ax = plt.subplot(1,1,1)
    
for s in simulations:
    s.power=lp
    s.p = 1e5
    equilibrium=np.mean(s.positions[9000:])
    s.run_ivp(equilibrium, 0, 5000, 10000)
    space_ax.plot(s.times, s.positions, label=str(int(round(s.normalized_radius)))+'$\mu$m')  #Change 

space_ax.set_xlabel('Time (ms)', fontsize='large')
space_ax.set_ylabel('Displacement ($\mu$m)', fontsize='large')
space_ax.tick_params(labelsize='large')
space_ax.legend()

for s in simulations:
    z0 = s.positions[0]
    z1 = s.positions[-1]
    dz = z0-z1
    if s.numerical_pressure_option:
        eq_ax.plot(int(round(s.normalized_radius)), dz, 'b*')
    else:
        eq_ax.plot(int(round(s.normalized_radius)), dz, 'bo')


eq_ax.set_ylabel('Change in equilibrium position ($\mu$m)', color='b', fontsize='large')
eq_ax.set_xlabel('Particle radius ($\mu$m)', fontsize='large') 
eq_ax.set_title('Change in equilibrium position for laser power decrease from 300mW to 150mW') #Change
plt.show()

