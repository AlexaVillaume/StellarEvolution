'''
'''
import math
import ../StellarEvolutionCode/opacity_interpolation
import ../StellarEvolutionCode/calc_density

"""
Constants
"""
G = 6.673e-8           # km^3 kg^-1 s^-2
c = 3e5                # km s^-1
# Should be taken from chatper 2.3
# but using known values for now
# use density_c = 162.2 g cm^-3 to find pressure_c
pressure_c =
temperature_c = 1.57e7 # K
# taken for the sun
total_lum = 4e26       # Watts
total_radius = 7e5     # km
# just a random small number
mass_initial = 1e-6

def load1(mass, radius, luminosity):
    """
    Computing the initial guesses for pressure and
    temperture at the core.

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
    opacity_c = opacity_interpolation.opacity_is(temperature_c, density_c)
    # Need to check on on the conditions
    # If convective
    if opacity_c < (2./5.):
        term1 = -(math.pi/6)**(1/3)
        term2 =  (del_ad_c*density_c**(4/3))/pressure_c
        temperature = math.exp(term1*G*term2*mass**(2/3) - math.log(temperature_c))
    # if not convective, then raditive
    else:
        term1 = -1/(2*a*c)
        term2 = (3/4*math.pi)**(2/3)
        term3 = densitity_c*(e_n - e_v + e_g)
        temperature = term1*term2*term3*(rho_c**(4/3))*(mass**(2/3)) + temperature_c**4

    return [pressure, temperature]

"""
Computing the initial guesses at the surface.

guess outer radius and ltot

return radius and l.
"""
def load2():
    # define surface based on tau=2/3
    # use equation 11.13, solved for R, the ideal gas law to get P?
    #kP=2g/3

    radius =
    # Use eddington approx. to find T, using opacity_interpolation to find tau
    opacity = opacity_interpolation.opacity_is()
    t_eff = ((3/4)*(total_lum/(4*math.pi*(radius**2)*sigma))*(opacity + (2/3)))**(1/4)


def derivs():
