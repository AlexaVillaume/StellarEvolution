import sys
import math
import numpy as np
import matplotlib.pyplot as plt

k = 1.38e-16
m_h = 1.67e-24
mu = 0.35

def calc_t(pressure, density):
    temperature = []
    for p, d in zip(pressure, density):
        temperature.append((p*mu*m_h)/(d*k))
    return temperature

radius_j = 6.9e4    # km
radius = np.linspace(0,1e5,1e6)
blargh = [value/radius_j for value in radius]
blargh = np.asarray(blargh)

jupiter_data = np.loadtxt('Jup_Guillot99.dat')

j_radius = [value for value in jupiter_data[:,1]]
j_pressure = [10**value for value in jupiter_data[:,2]]
j_temp = [10**value for value in jupiter_data[:,3]]
j_rho = [10**value for value in jupiter_data[:,4]]

density = 4.35*(np.sin(blargh*np.pi)/(blargh*np.pi))
pressure = 2.082e12*(density)**2
temperature = calc_t(pressure, density)

plt.plot(j_pressure, j_temp)
plt.plot(pressure, ((pressure*mu*m_h)/(density*k)), lw=5)
plt.xlabel('Pressure')
plt.ylabel('Temperature')
plt.xscale('log')
plt.yscale('log')
plt.show()
