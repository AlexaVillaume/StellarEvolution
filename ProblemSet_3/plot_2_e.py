import sys
import math
import numpy as np
import matplotlib.pyplot as plt

radius_j = 6.9e4    # km
radius = np.linspace(0,1e5,1e6)
blargh = [value/radius_j for value in radius]
blargh = np.asarray(blargh)

jupiter_data = np.loadtxt('Jup_Guillot99.dat')

j_radius = [value for value in jupiter_data[:,1]]
j_pressure = [10**value for value in jupiter_data[:,2]]
j_rho = [10**value for value in jupiter_data[:,4]]

plt.plot(j_radius, j_rho)
plt.plot(blargh, 4.35*(np.sin(blargh*np.pi)/(blargh*np.pi)))
plt.xlabel('Radius')
plt.ylabel('Density')
plt.show()


plt.plot(j_pressure, j_rho)
plt.plot(2.082e12*(4.35*(np.sin(blargh*np.pi)/(blargh*np.pi)))**2, 4.35*(np.sin(blargh*np.pi)/(blargh*np.pi)))
plt.xlabel('Pressure')
plt.ylabel('Density')
plt.show()
