'''
'''
import sys
import math
import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
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

def pressure_from_ideal(density, t_eff, mu):
    return (density*k*t_eff)/(mu*mass_h)

def calc_other_pressure(opacity):
    return (2*G*total_mass)/(3*opacity*total_radius**2)

def find_intersection(function1, function2, xvalues):
    return fsolve(lambda x: function1(x) - function2(x), xvalues)

# Dummy function for now
def calc_e_n(density, temperature):
    return 1

def outward_start(mass):
    """
    Computing the initial guesses for pressure and
    temperture at the core.

    At the core radius=0 and luminosity=0.

    We'll start at some tiny value of m, with an
    assumed pressure_c and temperature_c.
    """
    density_c = calc_density.density_is(math.log10(temperature_c), math.log10(pressure_c), X, Y)
    e_n = calc_e_n(density_c, temperature_c)

    term1 = -(3.*G)/(8.*math.pi)
    term2 = (density_c*4.*math.pi/3.)**(4./3.)
    term3 = mass**(2./3.)
    pressure = term1*term2*term3 + pressure_c
    luminosity_c = (e_n)*mass

    opacity_c = 10**(opacity_interpolation.opacity_is(math.log10(temperature_c), math.log10(density_c)))
    del_rad = calc_del_rad(density_c, pressure_c, temperature_c, opacity_c, luminosity_c, mass)
    if del_rad >= del_ad:
        # if temperature gradiant is greater than del_ad, then convective (because blobs are unstable)
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
def inward_start(test=True):
    t_eff = calc_teff(total_radius, total_lum)
    mu = calc_density.mu_is(X, Y)

    # Make an array of density values
    densities = np.linspace(1e-9,density_surface,1e5)
    # Calculate corresponding pressure arrays
    pressures1 = np.asarray([pressure_from_ideal(density, t_eff, mu) for density in densities])
    # Calculate an array of opacity values for computing pressure2
    opacities = [10**opacity_interpolation.opacity_is(math.log10(t_eff), math.log10(density)) for density in densities]
    pressures2 = np.asarray([calc_other_pressure(opacity) for opacity in opacities])

    intersection_index = np.abs(pressures1 - pressures2).argmin(0)
    # Need to address what to do if more then one index is returned
    surface_pressure = (pressures1[intersection_index] + pressures2[intersection_index]) / 2
    surface_density = densities[intersection_index]

    if test:
        plt.plot(densities, pressures2, lw=2)
        plt.plot(densities, pressures1, lw=2)
        plt.plot(surface_density, surface_pressure, ls='none', marker='o')
        plt.xlabel('Density')
        plt.ylabel('Pressure')
        plt.xscale('log')
        plt.yscale('log')
        plt.show()

    return [surface_pressure, t_eff, total_radius, total_lum]

def derivs(pressure, temperature, radius, luminosity, mass):
    """
    The mass given should be the enclosed mass, deal with that
    outside the function.
    """
    density = calc_density.density_is(math.log10(temperature), math.log10(pressure), X, Y)
    opacity = 10**opacity_interpolation.opacity_is(math.log10(temperature), math.log10(density))

    dpressure_dm = -((G)/(4*math.pi))*((mass)/(radius**4))
    dradius_dm = (1./(4.*math.pi))*(1./(density*radius**2))
    dluminoisty_dm = calc_e_n(density, temperature)

    del_rad = calc_del_rad(density, pressure, temperature, opacity, luminosity, mass)
    if del_rad >= del_ad:
        dtemperature_dm = -((G*mass*temperature)/(4*math.pi*pressure))*del_ad
    else:
        dtemperature_dm = -((G*mass*temperature)/(4*math.pi*pressure))*del_rad

    return [dpressure_dm, dtemperature_dm, dradius_dm, dluminoisty_dm]

core =  outward_start(mass_initial)
surface =  inward_start()
print derivs(core[0], core[1], core[2], core[3], mass_initial)
print derivs(surface[0], surface[1], surface[2], surface[3], total_mass)
