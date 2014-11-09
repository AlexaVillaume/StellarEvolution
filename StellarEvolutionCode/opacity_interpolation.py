'''
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

