# 1D-Optical-Levitation-Trap

## Basic Operation

First initialise a 'trapped particle' object, defined in trap_model.py:

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
