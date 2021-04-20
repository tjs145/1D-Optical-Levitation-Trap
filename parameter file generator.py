import numpy as np
import os

os.chdir('varying_radii')

radii = np.arange(2e-6, 10e-6, 1e-6)
density = 963  #For silicon oil

for r in radii:
    mass = density*(4/3)*np.pi*r**3
    if r>=0e-6:
        parameters = np.array([['Particle radius', str(r)], ['Particle mass', str(mass)], ['Wavelength', '532e-9'], ['Beam divergence', '0.08309'], ['n_sphere', '1.52'], ['Drag', '1'], ['Brownian motion', '0'], ['radiation pressure', '1']])
    else:
        parameters = np.array([['Particle radius', str(r)], ['Particle mass', str(mass)], ['Wavelength', '532e-9'], ['Beam divergence', '0.08309'], ['n_sphere', '1.52'], ['Drag', '1'], ['Brownian motion', '0'], ['radiation pressure', '0']])
    np.savetxt('Standard radius'+str(int(round(r*1e6))).zfill(3)+'.csv', parameters, fmt='%s', delimiter=': ')





##os.chdir('varying_n')
##
##indices = np.arange(1.2, 2.6, 0.1)
##density = 963  #For silicon oil
##r=5e-6
##
##for n in indices:
##    mass = density*(4/3)*np.pi*r**3
##    if r>=20e-6:
##        parameters = np.array([['Particle radius', str(r)], ['Particle mass', str(mass)], ['Wavelength', '532e-9'], ['Beam divergence', '0.037'], ['n_sphere', str(n)], ['Drag', '1'], ['Brownian motion', '0'], ['radiation pressure', '1']])
##    else:
##        parameters = np.array([['Particle radius', str(r)], ['Particle mass', str(mass)], ['Wavelength', '532e-9'], ['Beam divergence', '0.037'], ['n_sphere', str(n)], ['Drag', '1'], ['Brownian motion', '0'], ['radiation pressure', '0']])
##    np.savetxt('RI '+str(int(round(10*n)))+'.csv', parameters, fmt='%s', delimiter=': ')

##os.chdir('varying_density')
##
##densities = np.arange(250, 7000, 250)
####density = 963  #For silicon oil
##r=5e-6
##n=1.5
##
##for density in densities:
##    mass = density*(4/3)*np.pi*r**3
##    if r>=20e-6:
##        parameters = np.array([['Particle radius', str(r)], ['Particle mass', str(mass)], ['Wavelength', '532e-9'], ['Beam divergence', '0.037'], ['n_sphere', str(n)], ['Drag', '1'], ['Brownian motion', '0'], ['radiation pressure', '1']])
##    else:
##        parameters = np.array([['Particle radius', str(r)], ['Particle mass', str(mass)], ['Wavelength', '532e-9'], ['Beam divergence', '0.037'], ['n_sphere', str(n)], ['Drag', '1'], ['Brownian motion', '0'], ['radiation pressure', '0']])
##    np.savetxt(str(int(round(density))).zfill(4)+'.csv', parameters, fmt='%s', delimiter=': ')


