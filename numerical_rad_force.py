import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
import time

def numerical_pressure(P, wl, n_0, n_s, W_0, R, z):
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
