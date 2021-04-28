# -*- coding: utf-8 -*-
"""
Created on Tue Dec 08 15:25:46 2020

@author: tjsmi
"""

from scipy.integrate import odeint, solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from numerical_rad_force import numerical_pressure
from effective_reflectivity import calculate_reflectivity
from photophoretic_force import photophoretic_force_calculator

class Trapped_particle:
    def __init__(self, pressure, laser_power, parameter_filename='default_parameters.txt', length_scale=1e-6, time_scale=1e-3):
        '''
        Initialises the trapped particle object by reading parameters from the parameter file and assigns various variables
        Pressure in Pa (Can be constant or function of time)
        Laser power in W (Can be constant or function of time)
        parameter_filename - The name of the file containing various model parameters (must be stored in current working directory), see 'default_parameters.txt' for an example.
        length_scale and time_scale are the characteristic scales of the system, used to normalize the units to reduce floating point error. Floating point error is negligible and the code has only been tested using the default values so we wouldn't recommend changing these.
        '''
        self.parameters = dict(np.loadtxt(parameter_filename, dtype=str, delimiter=': '))
        viscosity = 1.85e-5
        molecular_mass = 4.8e-26
        temperature = 293
        radius = float(self.parameters['Particle radius'])
        self.normalized_radius = radius/length_scale
        self.volume = (4/3)*np.pi*radius**3
        wavelength = float(self.parameters['Wavelength'])
        self.normalized_wavelength = wavelength/length_scale
        self.m = self.volume* float(self.parameters['Particle density'])
        divergence = float(self.parameters['Beam divergence'])
        self.norm_beam_waist = self.normalized_wavelength/(np.pi*divergence)
        self.n = float(self.parameters['n_sphere'])
        self.drag_option = bool(int(self.parameters['Drag']))
        self.brownian_option = bool(int(self.parameters['Brownian motion']))
        self.numerical_pressure_option = bool(int(self.parameters['radiation pressure']))
        self.normalized_stokes = (6*np.pi*viscosity*radius*time_scale)/self.m
        self.g = 9.81*time_scale**2/length_scale
        self.MFP_factor = (viscosity/length_scale)*np.sqrt((np.pi*1.38e-23*temperature)/(2*molecular_mass))
        self.B = np.sqrt((2*time_scale*1.38e-23*temperature)/(6*np.pi*viscosity*radius))/length_scale
        self.p = pressure
        self.power = laser_power
        self.time_scale = time_scale
        self.file = parameter_filename

        #the parameters below are for the photphoretic force caluclator function
        self.PF_option = bool(int(self.parameters['Photophoretic']))
        self.temperature=float(self.parameters['Gas temperature'])
        self.n_i= float(self.parameters['n_sphere_imag'])
        self.k_p= float(self.parameters['Sphere thermal conductivity'])
        self.photophoretic_norm_factor=(time_scale**2)/(self.m*length_scale)

        if not self.numerical_pressure_option: #calculated the effcetive reflctivity if using the faster mirror model
            self.reflectivity = calculate_reflectivity(self.normalized_wavelength, 1, self.n, self.norm_beam_waist, self.normalized_radius)    
            self.power_factor = (2*self.reflectivity*time_scale**2)/(self.m*3e8*length_scale) #Converts units of power into normalised units of optical force

        else:
            self.power_factor = (time_scale**2)/(self.m*length_scale) #Normalises the units of optical force

        
    def pressure(self, t):
        '''
        It pressure is given as a 'callable' function, it's value is returned at the current time (t)
        It pressure is constant, then its value is returned
        '''
        if callable(self.p):
            return self.p(t)
        else:
            return self.p
        
    def normalized_MFP(self, t):
        '''
        Returns the mean-free-path (MFP) in normalised units. 
        MFP_factor is deifned in the __init__ function based on a standard definition for MFP.
        '''
        return self.MFP_factor/self.pressure(t)

    def kn(self, t):
        '''
        Returns the Knudsen number of the system.
        '''
        return self.normalized_MFP(t)/self.normalized_radius
        
    def drag_correction(self, t):
        '''
        A correction factor for Stokes law in the case of a rarefied gas.
        (Millikan, R.A., 1923. Stokes' Law of Fall Completely Corrected. Proceedings of the National Academy of Sciences of the United States of America, 9(3), p.67.)
        '''
        return 1/(1+self.kn(t)*(0.864+0.29*np.exp(-1.25/self.kn(t))))

    def drag(self, t):
        '''
        Returns the drag coefficient calculated using Stokes law and Millikans correction, if drag is switched on.
        Drag can be switched off in the parameter file, in which case this functions returns 0.
        '''
        if self.drag_option:
            return self.normalized_stokes*self.drag_correction(t)
        else:
            return 0

    def brownian(self, t):
        '''
        Can be switched on/off in the parameter file.
        If on, calculates random perturbations to the velocity of the particle, based on Einsteins theory of brownian motion.
        (Renn, J. (2005). Einstein's invention of Brownian motion. Annalen der Physik, 14(S1), 23-37.)
        This expresses the RMS change in position after a time delta_t. 
        We converted this into an RMS velocity perturbation by dividing by delta_t.
        delta_t refers to the average time between steps of the DE solver.
        If off, returns  0.
        '''
        if self.brownian_option:
            return np.random.normal(0, self.B/np.sqrt(self.delta_t*self.drag_correction(t)))
        else:
            return 0

    def laser_power(self, t):
        '''
        It power is given as a 'callable' function, it's value is returned at the current time (t)
        It power is constant, then its value is returned
        '''
        if callable(self.power):
            return self.power(t)
        else:
            return self.power

    def norm_beam_radius(self, z):
        '''
        Returns the radius of a gaussian beam at a distance z from the focus.
        '''
        return self.norm_beam_waist*np.sqrt(1+((self.normalized_wavelength*z)/(np.pi*self.norm_beam_waist**2))**2)

    def radiation_pressure(self, z, t):
        '''
        Returns the acceleration on the particle due to the scattering force.
        If the numerical method has been selected in the parameter file, this calculated from the numerical_rad_force file.
        If the numerical method for calculating optical force is not selected, the simpler mirror model is used.
        '''
        if self.numerical_pressure_option:
            return self.power_factor*numerical_pressure(self.laser_power(t), self.normalized_wavelength, 1, self.n, self.norm_beam_waist, self.normalized_radius, z)
        else:
            return self.laser_power(t)*self.power_factor*(1-np.exp(-(2*self.normalized_radius**2)/(self.norm_beam_radius(z)**2)))
        
    def photophoretic_force(self, z, t): #will not work correctly if imaginary part of RI is too large (greather than 1e-4)
        '''
        Caculates the acceleration due to photophoretic force using the module in 'photophoretic_force.py'
        Returns 0 if photophoretic forces are switched off in the parameter file.
        '''
        if self.PF_option:
            return self.photophoretic_norm_factor*photophoretic_force_calculator(self.pressure(t), self.temperature, self.laser_power(z, t), self.normalized_wavelength*1e-6, self.n, self.n_i, self.k_p, self.norm_beam_waist*1e-6, self.normalized_radius*1e-6, z*1e-6)
        else:
            return 0

    def gravity(self):
        return self.g

    def air_density(self, t):
        '''
        Returns the air density based on the assumption that density is proportional to pressure
        '''
        return 1.2e-5*self.pressure(t)

    def buoyancy(self, t):
        return self.g*self.volume*self.air_density(t)*(time_scale**2)/(self.m*length_scale)

    def ivp_func(self, t, y): 
        '''
        Defines the function that returns the first derivative of vector y (i.e. the differential equation that describes the system).
        Gets called by solve_ivp
        '''
        matrix = np.array([[0,1],[0, -self.drag(t)]])
        return np.matmul(matrix, y)+np.array([self.brownian(t), -1*self.photophoretic_force(y[0], t)+self.radiation_pressure(y[0], t)-self.gravity()+self.buoyancy(t)])

    def run_ivp(self, z0, v0, t_end, n, method='RK45'):
        '''
        Calculate numerical soultion for y (=[z,v]) in time range [0 - t_end]
        n is the number of data points to calculate and plot
        z0 and v0 are the initial values for normalized position and normalized velocity
        
        method argument describes which DE solver to use. We recommed:
        'RK45' - Good all round
        'Radau' - Good for systems with no random effects, so no brownian motion or random variations in laser power.
        'DOP853' - More accurate but slower. 
        
        Solution stored in 'self.times' and 'self.positions'.
        Power spectrum of solution stored in 'self.frequencies' and 'self.p_spectrum'
        '''
        y0 = np.array([z0, v0])
        t=np.linspace(0, t_end, n)
        self.delta_t = t_end/n
        sol=solve_ivp(self.ivp_func, [0, t_end], [z0, v0], t_eval=t, method=method)
        self.times=sol['t']
        y=sol['y']
        self.positions = y[0]
        self.velocities = y[1]

        self.p_spectrum = np.absolute(np.fft.fft(self.positions-np.mean(self.positions)))
        self.frequencies = np.fft.fftfreq(n, (t_end*self.time_scale)/n)

        print(self.file+' run complete')
