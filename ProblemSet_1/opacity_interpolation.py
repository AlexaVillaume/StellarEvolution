'''
Develop a program, by interpolation in a table, for caculating the opacity of a
gas with a given density and temperature, and with the chemical compoistion that you
will use for the stellar model.

What do you obtain for the opacity in (cm^2 g^-1) for:
    (a) log T = 6.3 K, log rho = 0.3 (cgs)
    (b) log T = 5.0 K, log rho = -4.0 (cgs)

The temperature range you will need for the stellar models runs from 5x10^3 to 3x10^7 K,
and the dnesity range from 10^-9 to 10^3 g cm^-3, depending on T.
'''

import sys
import math
import numpy as np
from scipy import interpolate

def density_r(logT, logRho):
    # R=density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
    rho = 10**logRho
    T = 10**logT
    return math.log10(rho/(T*1.e-6)**3)

# Read in the table and sort the values into the appropriate arrays
logTs = []
opacities = []
with open('test_table.dat', 'r') as f:
    for i, row in enumerate(f):
        cols = row.split()
        if i == 0:
            logRs = np.asarray([float(value) for value in cols[1:len(cols)]])
        else:
            logTs.append(float(cols[0]))
            if len(logRs) == len(cols[1:len(cols)]):
                opacities.append([float(value) for value in cols[1:len(cols)]])
            else:
                tmp_opacity = [float(value) for value in cols[1:len(cols)]]
                while True:
                    if len(tmp_opacity) == len(logRs):
                        break
                    else:
                        tmp_opacity.append(9.999)
                opacities.append(tmp_opacity)

logTs = np.asarray(logTs)
blargh = np.empty((70,19))
for i, blah in enumerate(opacities):
    blargh[i] = blah

function = interpolate.RectBivariateSpline(logTs, logRs, blargh)

# Test
print "For part (a): "
print "The R value: ", density_r(6.3, 0.3)
print "The interpolated opacity: ", function(6.3, density_r(6.3, 0.3))

print "For part (b): "
print "The R value: ", density_r(5.0, -4.0)
print "The interpolated opacity: ",function(5.0, density_r(5.0, -4.0))
