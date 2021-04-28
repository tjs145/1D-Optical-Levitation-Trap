import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
import time
import numpy as np

def calculate_reflectivity(wl, n_0, n_s, W_0, R):
    '''
    Calculates the effective reflectivity of the particle when using mirror approximation, for a system defined by the arguments:
    
    wl - wavelength of incident light
    n_0 - refractive index of surrounding medium
    n_s - refractiev index of the particle/sphere
    W_0 - size of the beam waist
    R - radius of the particle
    
    Compares radiation pressure exerted on a perfectly reflective mirror to radiation pressure predicted by Kim and Kim, 2017
    (Kim, S. B., & Kim, S. S. (2006). Radiation forces on spheres in loosely focused Gaussian beam: ray-optics regime. JOSA B, 23(5), 897-903.)
    '''
    P=1

    z_0=math.pi*(W_0**2)/wl # rayleigh range

        #angles are in radians
    def th2(th1): #Angle of refraction for a given angle of incidence, calculated using snells law
        return math.asin((n_0/n_s)*math.sin(th1))
        
    def p(th1): #The radial position for a given angle of incidence.
        return R*math.sin(th1)
    
    def r(th1):#power reflectance coefficient calculated using Fresnels equation
        th_2 = th2(th1)
        return 0.5*(((math.sin(th1-th_2)**2)/math.sin(th1+th_2)**2)+((math.tan(th1-th_2)**2)/math.tan(th1+th_2)**2))

    def t(th1):#transmittance coefficient
        return 1-r(th1) 

    def W(z):#beam radius
        return W_0*((1+(z/z_0)**2)**0.5)
    
    def I(th1,z):#intensity of a gaussian beam
        return ((2*P)/(math.pi*(W(z)**2)))*math.exp((-2*p(th1)**2)/(W(z)**2))
    
    def kernel_1(th1,z): #The expression to be integrated
        th_2 = th2(th1)
        re = r(th1)
        return (math.pi/3e8)*n_0*I(th1,z)*(1+re*math.cos(2*th1)-(t(th1)**2)*((math.cos(2*(th1-th_2))+re*math.cos(2*th1))/(1+(re**2)+2*re*math.cos(2*th_2))))*(R**2)*(math.sin(2*th1))


    numerical_forces = []   #Stores numerically calculated radiation pressures over a range of axial positions (z)
    z_values = np.arange(1000, 3000, 10)

    for z in z_values:
        numerical_forces.append(integrate.quad(kernel_1,0,math.pi/2,args=(z))[0])  #Calculates forces using the numerical method for the various z positions, by integrating over the pi/2 possible angles of incidence


    mirror_forces = []  #Stores the values for the radiation pressure on a perfect sphere for the same range of z positions
    
    def p_incident(z):
        return P*(1-np.exp(-(2*R**2)/(W(z)**2)))

    def force(z):
        return 2*p_incident(z)/3e8

    for z in z_values:
        mirror_forces.append(force(z))


    return np.mean(numerical_forces)/np.mean(mirror_forces)   #Returns the 'effective reflectivity' the factor of difference between the perfect mirror approximation and the numerical technique.
