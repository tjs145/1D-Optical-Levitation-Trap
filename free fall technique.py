import os

os.chdir('various radii')

files = os.listdir()

from trap_model import Trapped_particle
import matplotlib.pyplot as plt
import numpy as np

def lp(t):
    if t>5000 and t<5100:
        return 0
    else:
        return 0.3

simulations=[]

for file in files:
    simulations.append(Trapped_particle(1e5, lp, parameter_filename=file))

plt.figure(1) #Fig 1 shows position with time
space_ax = plt.subplot(1,1,1)

for s in simulations:
    s.run_ivp(1000,0,7000,10000)
    space_ax.plot(s.times, s.positions, label=str(round(s.normalized_radius))+'$\mu$m')


space_ax.set_xlabel('Time (ms)', fontsize='large')
space_ax.set_ylabel('Displacement ($\mu$m)', fontsize='large')
space_ax.tick_params(labelsize='large')
space_ax.legend()

plt.show()
