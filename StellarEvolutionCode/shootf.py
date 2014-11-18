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

mass_initial = 1e-6

def percent_difference(value1, value2):
    average = (math.fabs(value1) + math.fabs(value2))/2.
    return math.fabs(value1 - value2) / average

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
    density_surface = 10e-7
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

# I should make an input object
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

"""
Making testing suite.
"""
if __file__ == sys.argv[0]:
    core =  outward_start(mass_initial)
    surface =  inward_start()
    print derivs(core[0], core[1], core[2], core[3], mass_initial)
    print derivs(surface[0], surface[1], surface[2], surface[3], total_mass)

    # test on on multiple sets of stars (once I make star object)...abuse these functions
