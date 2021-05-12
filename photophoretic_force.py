# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 19:48:17 2021

@author: timkl
"""
import scipy.integrate as integrate
import math

def photophoretic_force_calculator(p, T, P, wl,n_s, n_s_i, k_p, W_0, R, z): #k_p is the spheres thermal conduction coefficent
    
    pr=p
    


    n_r=n_s #real part of the spheres refractive index
    n_i=n_s_i # imag part of the spheres RI

    r=R#sphere radius
    R=8.3145 #gas constant


    k=1.14 #thermal creep  coefficient
    M= 28.55e-3# M needs to be in kg per mol of substance, 

    vis=1.85e-5 #gas viscosity
    alpha=0.7 #accomodation coefficient 
    



    x=(2*math.pi*r)/wl #size parameter
    J=2*n_r*n_i*x*(((3*(n_r-1))/(8*(n_r**2)))-(2/5)*n_r*n_i*x) # asymmetry factor for when x is much more than one and x*n_i is much less than one
    #J=0.5
    c=((8*R*T)/(M*math.pi))**0.5 #mean speed of the gas molecules
    D=(math.pi/2)*(((math.pi/3)*k)**0.5)*((c*vis)/T) #the D coefficient
    pr_h=((2/alpha)**0.5)*(3/math.pi)*D*(T/r) #the pressure at which the PF force reaches its maximum value

    def W(z):#beam radius
        return W_0*((1+((wl*z)/(math.pi*W_0**2))**2)**0.5)
    def Int(z,p):#intensity of a gaussian beam
        return ((2*P)/(math.pi*(W(z)**2)))*math.exp((-2*p**2)/(W(z)**2))
    def kernel_1(p,z):#kernel for integrating the beam intensity over the sphere to find power incident on it
        return ((2*P)/(math.pi*(W(z)**2)))*math.exp((-2*p**2)/(W(z)**2))*2*math.pi*p
    def I_av(z):#returns the average intensity of light incident on the sphere when its center is at position z in the beam
        P_sphere=integrate.quad(kernel_1,0,r,args=(z))
        return (P_sphere[0])/(math.pi*r**2)

    def f_hat(z):
        #return D*((alpha/2)**0.5)*(((r**2)*J)/k_p)*I_av(z)
        return D*((alpha/2)**0.5)*(((r**2)*J)/k_p)*I_av(z)

    #evaluates the photpohretic force on the sphere, for the sphere with its center at position z along the z axis ()
    
    return f_hat(z)*(2/((pr/pr_h)+(pr_h/pr)))
