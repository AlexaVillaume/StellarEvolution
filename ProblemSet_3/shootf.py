'''
'''
import math
import ../StellarEvolutionCode/opacity_interpolation
import ../StellarEvolutionCode/calc_density

"""
Constants.

Characterstics for the star are all from the sun.
"""
G = 6.673e-8           # km^3 kg^-1 s^-2
c = 3e5                # km s^-1
pressure_c = 2.526e14  # g cm^-2
temperature_c = 1.57e7 # K
# taken for the sun
total_lum = 3.846e33   # erg s^-1
total_radius = 7e5     # km
total_mass =
# just a random small number
a = 7.56464e-15        # ergs cm^-2 K^-4
mass_h = 1.38e-16   #cgs
del_ad = 0.4
sigma  = 5.67e-5 #erg cm^-2 s^-1 K^-1
mass_initial = 1e-6
density_surface = 10e-7
X = 0.70
Y = 0.27

def percent_difference(value1, value2):
    average = (math.fabs(value1) + math.fabs(value2))/2
    return math.fabs(value1 - value2) / average

def calc_luminosity(radius, temperature)
    return 4*math.pi*(radius**2)*sigma*temperature**4

def calc_del_rad(density, pressure, temperature, luminosity)
    return (3/16*math.pi)*a*c*((G*density*luminosity*pressure)/(mass*temperature**4))

# Dummy function for now
def calc_e_n(density, temperature):
    return 1

def load1(mass, radius, luminosity):
    """
    Computing the initial guesses for pressure and
    temperture at the core.

    At the core radius=0 and luminosity=0.

    We'll start at some tiny value of m, with an
    assumed pressure_c and temperature_c.
    """
    density_c = calc_density.density_is(math.log10(temperature_c), math.log10(pressure_c), X, Y))
    e_n = calc_e_n(density_c, temperature_c)

    term1 = -(3*G)/(8*math.pi)
    term2 = (density_c*4*math.pi/3)**(4/3)
    term3 = mass**(2/3)
    pressure = term1*term2*term3 + pressure_c

    luminosity_c = (e_n)*mass

    opacity_c = opacity_interpolation.opacity_is(temperature_c, density_c)
    del_rad = calc_del_rad(density_c, pressure_c, temperature_c, luminosity_c)
    if del_rad >= del_ad:
        # if temperature gradiant is greater than del_ad, then convective
        term1 = -(math.pi/6)**(1/3)
        term2 =  (del_ad_c*density_c**(4/3))/pressure_c
        temperature = math.exp(term1*G*term2*mass**(2/3) - math.log(temperature_c))
    else:
        # if not convective, then radiative
        term1 = -1/(2*a*c)
        term2 = (3/4*math.pi)**(2/3)
        term3 = density_c*(e_n)
        temperature = term1*term2*term3*(rho_c**(4/3))*(mass**(2/3)) + temperature_c**4

    radius = (3/(4*math.pi*density_c))**(1/3) * mass**(1/3)

    return [pressure, temperature, radius, luminosity]

"""
Computing the initial guesses at the surface.

guess outer radius and ltot

return radius and l.
"""
def load2():
    # define surface based on tau=2/3
    mu = calc_density.mu_is(X, Y)
    kappa = opacity_interpolation.opacity_is(t_eff, density)

    pressure1 = (density/(mu*mass_h))*k*t_eff
    pressure2 = ((G*total_mass)/total_radius)*(2./3.)*(1./kappa)

    while True:
        if percent_difference(pressure1, pressure2) > 0.01:

        else:
            break
    #radius =

    return [pressure, temperature, radius, luminosity]

def derivs(pressure, temperature, radius, luminosity):
    d_pressure = -((G)/(4*math.pi))*((total_mass)/(radius**4))
    d_radius = (1/(4*math.pi))*(1/(density*radius**2))
    d_luminoisty = e_n

    density = calc_density.density_is(math.log10(temperature), math.log10(pressure), X, Y))
    del_rad = calc_del_rad(density, pressure, temperature)
    if del_ad < del_rad:
        d_temperature = *del_ad
    else:
        d_temperature = *del_rad

    return [d_pressure, d_temperature, d_luminosity, d_radius]
