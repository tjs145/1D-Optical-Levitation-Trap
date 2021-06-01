# 1D-Optical-Levitation-Trap

## Basic Operation

First initialise a 'trapped particle' object, defined in 'trap_model.py':

    t = Trapped_particle(1e5, 0.4, parameter_filename='default_parameters.txt', length_scale=1e-6, time_scale=1e-3) 
    
Arguments:

Pressure in Pa (Can be constant or function of time)

Laser power in W (Can be constant or function of time)

parameter_filename - The name of the file containing various model parameters (must be stored in current working directory), see 'default_parameters.txt' for an example.

length_scale and time_scale are the characteristic scales of the system, used to normalize the units to reduce floating point error. 
Floating point error is negligible and the code has only been tested using the default values so we wouldn't recommend changing these.


Run a simulation with the trapped particle using 'run_ivp':

    t.run_ivp(1000,0,3000,10000)

Arguments:

z0 - Starting position relative to trapping beam focus in normalised units (so microns when length_scale=1e-6)

v0 - Starting velocity in normalised units (so default is micrometers per millisecond)

t_end - How long to simulate the particle motion for in normalised units (default is milliseconds)

n - Number of points in the solution (typically around 1-10 times the number of seconds being simulated gives suitable resolution without taking too long)

The solution for (z, t) is then stored in two arrays:

    t.positions
    t.time
    
which can be visualised using matplotlib:

    import matplotlib.pyplot as plt
    plt.plot(t.times, t.positions)
    plt.show()
      
The power spectrum is also calculated and stored in:

    t.p_spectrum
    t.frequencies


## Further use

'equilibrium change for reduced laser power.py' performs these simulatons automatically for particles of various radii, defined in the 'various radii' directory. For each particle, a simulation is performed with a constant 300mW laser, so the equilibrium position of the particle is found. A second simulation is then run, where the laser power is halved to 150Mw after 1 second. The position of the particle is plotted against time for this second simulation. The code also automatically calculates the change in equilibrium position after the change in laser power, and plots this against particle radius.

'free fall technique.py' tests the possibility of estimating trapped particle size based on how far it falls when the laser is briefly switched off. This code simulates motion for particles of various radii, trapped in a 300mW laser which is switched off for 0.1 seconds after 5 seconds. The particle positions with time are plotted. 

'frequency sweep modelling.py' simulates the motion of differently sized particles when the laser power is modulated at increasing frequencies. The amplitude of motion for each particle, at each frequency, is calculated automatically and plotted. The frequency at which the amplitude of oscillation start to decrease depends on particle size, as well as the amplitude at low frequencies. This code can take a long time to produce results (30+ mins).
