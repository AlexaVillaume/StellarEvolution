'''

Develop a program, by interpolation in a table, for caculating the opacity of a
gas with a given density and temperature, and with the chemical compoistion that you
will use for the stellar model.

What do you obtain for the opacity in (cm^2 g^-1) for:
    (a) log T = 6.3 K, log rho = 0.3 (cgs)
    (a) log T = 5.0 K, log rho = -4.0 (cgs)

The temperature range you will need for the stellar models runs from 5x10^3 to 3x10^7 K,
and the dnesity range from 10^-9 to 10^3 g cm^-3, depending on T.

'''
import numpy as np
from astropy.table import Table
from scipy.interpolate import griddata
from scipy import interpolate

opacities = Table.read('test_table.dat', format='ascii.no_header', guess=False)

# Okay, the trick is getting the values read in correctly
# so I can access them all as arrays

# Need to get the temp values, r values into arrays
# Need to get the opacity values into an array-like structure
function = interpolate.interp2d(temp, r, opacity)

# Input
# Correct the input density values to "R" used in the table
