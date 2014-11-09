'''
'''
import math
import ../StellarEvolutionCode/opacity_interpolation
import ../StellarEvolutionCode/calc_density


"""
Constants
"""
G =
c =
# chatper 2.3
pressure_c =
temperature_c =
total_lum =
total_radius =
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
    # or convective by using the
    density_c = calc_density.density_is(math.log10(temperature_c, pressure_c, X, Y))
    opacity_c = opacity_interpolation.opacity_s(temperature_c, density_c)
    # Need to check on on the conditions
    if opacity_c < (2./5.):
        term1 = -(math.pi/6)**(1/3)
        term2 =
        temperature = math.exp(term1*G*term2*mass**(2/3) - math.log(temperature_c))
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

