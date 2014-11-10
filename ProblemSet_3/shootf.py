'''
'''
import sys
import math

import opacity_interpolation
import calc_density

"""
Constants:
"""
G = 6.673e-8           # cm^3 g^-1 s^-2
c = 3e10               # cm s^-1
sigma  = 5.67e-5       # erg cm^-2 s^-1 K^-1
a = 7.56464e-15        # ergs cm^-2 K^-4
mass_h = 1.673e-24     # cgs
del_ad = 0.4
k = 1.38e-16           # ergs K^-1

# For the Sun
pressure_c = 2.526e14  # g cm^-2
temperature_c = 1.57e7 # K
total_lum = 3.846e33   # erg s^-1
total_radius = 7e10     # cm
total_mass = 1.98e33    #g

mass_initial = 1e-6
density_surface = 10e-7
X = 0.70
Y = 0.27

def percent_difference(value1, value2):
    average = (math.fabs(value1) + math.fabs(value2))/2.
    return math.fabs(value1 - value2) / average

def calc_luminosity(radius, temperature):
    return 4.*math.pi*(radius**2)*sigma*temperature**4.

def calc_teff(radius, luminosity):
    return (luminosity/(4*math.pi*sigma*radius**2))**(1./4.)

def calc_del_rad(density, pressure, temperature, opacity, luminosity, mass):
    term1 = 3/(16*math.pi*a*c*G)
    term2 = (opacity*luminosity*pressure)/(mass*temperature**4)
    return term1*term2

# Dummy function for now
def calc_e_n(density, temperature):
    return 1

def load1(mass):
    """
    Computing the initial guesses for pressure and
    temperture at the core.

    At the core radius=0 and luminosity=0.

    We'll start at some tiny value of m, with an
    assumed pressure_c and temperature_c.
    """

    density_c = calc_density.density_is(math.log10(temperature_c), math.log10(pressure_c), X, Y)
    e_n = calc_e_n(density_c, temperature_c)

    # Calculate initial, tiny radius
    radius_initial = 3./(4*math.pi*density_c)**(1./3.)*mass**(1./3.)

    term1 = -(3.*G)/(8.*math.pi)
    term2 = (density_c*4.*math.pi/3.)**(4./3.)
    term3 = mass**(2./3.)
    pressure = term1*term2*term3 + pressure_c
    luminosity_c = (e_n)*mass

    opacity_c = 10**(opacity_interpolation.opacity_is(math.log10(temperature_c), math.log10(density_c)))
    del_rad = calc_del_rad(density_c, pressure_c, temperature_c, opacity_c, luminosity_c, mass)
    if del_rad >= del_ad:
        # if temperature gradiant is greater than del_ad, then convective
        term1 = -(math.pi/6.)**(1./3.)
        term2 =  (del_ad_c*density_c**(4./3.))/pressure_c
        temperature = math.exp(term1*G*term2*mass**(2./3.) - math.log(temperature_c))
    else:
        # if not convective, then radiative
        term1 = -1./(2.*a*c)
        term2 = (3/(4*math.pi))**(2./3.)
        term3 = opacity_c*(e_n)*density_c**(4./3.)
        term4 = mass**(2./3.)
        temperature = (temperature_c**4 - term1*term2*term3*term4)**(1./4.)

    radius = (3./(4.*math.pi*density_c))**(1./3.) * mass**(1./3.)

    return [pressure, temperature, radius, luminosity_c]

"""
Computing the initial guesses at the surface.

guess outer radius and ltot

return radius and l.
"""
def load2():
    t_eff = calc_teff(total_radius, total_lum)
    mu = calc_density.mu_is(X, Y)
    opacity_s = 10**opacity_interpolation.opacity_is(math.log10(t_eff), math.log10(density_surface))

    # These are more different than I would expect
    pressure1 = (density_surface*k*t_eff)/(mu*mass_h)
    pressure2 = (2*G*total_mass)/(3*opacity_s*total_radius**2)

    #while percent_difference(pressure1, pressure2) > 0.01:

    #    else:
    #        break
    #radius =

    #return [pressure, temperature, radius, luminosity]

def derivs(pressure, temperature, radius, luminosity):
    density = calc_density.density_is(math.log10(temperature), math.log10(pressure), X, Y)
    opacity = 10**opacity_interpolation.opacity_is(math.log10(temperature), math.log10(density))

    d_pressure = -((G)/(4*math.pi))*((total_mass)/(radius**4))
    d_radius = (1./(4.*math.pi))*(1./(density*radius**2))
    d_luminoisty = calc_e_n(density, temperature)

    # What mass is this supposed to be?
    del_rad = calc_del_rad(density, pressure, temperature, opacity, luminosity, total_mass)
    if del_ad < del_rad:
        d_temperature = -((G*total_mass*temperature)/(4*math.pi*pressure))*del_ad
    else:
        d_temperature = -((G*total_mass*temperature)/(4*math.pi*pressure))*del_rad

    return [d_pressure, d_temperature, d_radius, d_luminoisty]

test =  load1(mass_initial)
#print load2()
print derivs(test[0], test[1], test[2], test[3])
