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
import numpy as np
from scipy import interpolate

# Read in the table and sort the values into the appropriate arrays
logT = []
opacity = []
with open('test_table.dat', 'r') as f:
    for i, row in enumerate(f):
        cols = row.split()
        if i == 0:
            logR = np.asarray([float(value) for value in cols[1:len(cols)]])
        else:
            logT.append(float(cols[0]))
            if len(logR) == len(cols[1:len(cols)]):
                opacity.append([float(value) for value in cols[1:len(cols)]])
            else:
                tmp_opacity = [float(value) for value in cols[1:len(cols)]]
                while True:
                    if len(tmp_opacity) == len(logR):
                        break
                    else:
                        tmp_opacity.append(9.999)
                opacity.append(tmp_opacity)

logT = np.asarray(logT)
blargh = np.empty((70,19))
for i, blah in enumerate(opacity):
    blargh[i] = blah


# Nee
function = interpolate.RectBivariateSpline(logT, logR, blargh)
print function

# Input
# Correct the "R" values to be density values
# R=density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
