'''
'''
import math
import opacity_interpolation


"""
Constants
"""
G =
c =
rho_c =
pressure_c =
temperature_c =
mass_initial =

def load1(mass, radius, luminosity):
    """
    Computing the initial guesses at the core.

    At the core radius=0 and luminosity=0.

    We'll start at some tiny value of m, with an
    assumed pressure_c and temperature_c.
    """
    term1 = -(3*G)/(8*math.pi)
    term2 = (rho_c*4*math.pi/3)**(4/3)
    term3 = mass**(2/3)
    pressure = term1*term2*term3 + pressure_c
    # Test whether the gas is raidative
    # or convective
    if convective:
        temperature =
    else:
        term1 = -1/(2*a*c)
        term2 = (3/4*math.pi)**(2/3)
        term3 =
        temperature = term1*term2*term3*(rho_c**(4/3))*(mass**(2/3)) + temperature_c**4

    return [pressure, temperature]

"""
Computing the initial guesses at the surface.

guess outer radius and ltot

return radius and l.
"""
def load2():


#L_tot = 4*pi*R^2*sigma*Teff

# Use eddington approx. to find T, look up k
#kP=2g/3

