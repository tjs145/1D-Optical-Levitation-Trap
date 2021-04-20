# 1D-Optical-Levitation-Trap

## Basic Operation

First initialise a 'trapped particle' object:

    t = Trapped_particle(1e5, 0.4, parameter_filename='default_parameters.txt') 
    
Arguments:

Pressure in Pa (Can be constant or function of time)

Laser power in W (Can be constant or function of time)

Parameter_filename - The name of the file containing various model parameters (must be stored in current working directory), see 'default_parameters.txt' for an example.


Run a simulation with the trapped particle using 'run_ivp':

    t.run_ivp(1000,0,3000,10000)
