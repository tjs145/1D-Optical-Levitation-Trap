# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 20:03:59 2021

@author: timkl

This code is used to calculate the axial component of the radiation force exerted on a transparent homogenous 
sphere by a gaussian beam for when the particle is located on the beam axis. The theoretical approach behind it is based on the paper:
"Radiation forces on spheres in loosely focused Gaussian beam: ray-optics regime"-Sang Bok Kim and Sang Soo Kim.

"""

import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
import time

def numerical_pressure(P, wl, n_0, n_s, W_0, R, z):
    '''
    calculate the axial component of the radiation pressure in newtons for a transparent sphere 
    exerted by a gaussian laser beam when the particle is located on the beam axis. P is the beam power, wl is the beam wavelength, 
    n_0 is the refractive index of air, n_s is the refractive index of the material that makes up
    the target, W_0 is the beam waist, and z is the axial position with origin being the beam waist.
    '''
    z_0=math.pi*(W_0**2)/wl # rayleigh range

        #angles are in radians
    def th2(th1):
        return math.asin((n_0/n_s)*math.sin(th1))
        
    def p(th1):
        return R*math.sin(th1)
    def r(th1):#power reflectance coefficient
        th_2 = th2(th1)
        return 0.5*(((math.sin(th1-th_2)**2)/math.sin(th1+th_2)**2)+((math.tan(th1-th_2)**2)/math.tan(th1+th_2)**2))

    def t(th1):#transmittance coefficient
        return 1-r(th1) 


    def W(z):#beam radius
        return W_0*((1+(z/z_0)**2)**0.5)
    def I(th1,z):#intensity of a gaussian beam
        return ((2*P)/(math.pi*(W(z)**2)))*math.exp((-2*p(th1)**2)/(W(z)**2))
    def kernel_1(th1,z):
        th_2 = th2(th1)
        re = r(th1)
        return (math.pi/3e8)*n_0*I(th1,z)*(1+re*math.cos(2*th1)-(t(th1)**2)*((math.cos(2*(th1-th_2))+re*math.cos(2*th1))/(1+(re**2)+2*re*math.cos(2*th_2))))*(R**2)*(math.sin(2*th1))
      
    return integrate.quad(kernel_1,0,math.pi/2,args=(z))[0]
