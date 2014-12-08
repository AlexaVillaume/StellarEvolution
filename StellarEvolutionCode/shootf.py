'''
'''
import sys
import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import newton_krylov
import matplotlib.pyplot as plt
import opacity_interpolation
import calc_density
import utilities

def outward_start(star, mass, test=True):
    """
    Computing the initial guesses for pressure and
    temperture at the core.

    At the core radius=0 and luminosity=0.

    We'll start at some tiny value of m, with an
    assumed pressure_c and temperature_c.
    """
    core_density = calc_density.density_is(math.log10(star.core_temp), math.log10(star.core_pressure), \
                star.hydrogen_mass, star.helium_mass)

    e_n = utilities.calc_e_n(star, core_density, star.core_temp)

    term1 = -(3.*utilities.gravitational_constant)/(8.*math.pi)
    term2 = (core_density*4.*math.pi/3.)**(4./3.)
    term3 = mass**(2./3.)
    pressure = term1*term2*term3 + star.core_pressure
    core_luminosity = e_n*mass

    core_opacity = 10**(opacity_interpolation.opacity_is(math.log10(star.core_temp), \
                math.log10(core_density)))
    del_rad = star.calc_del_rad(core_density, star.core_pressure, star.core_temp, \
                core_opacity, core_luminosity, mass)
    if del_rad >= utilities.del_adiabatic:
        if test:
            print "Convective"
        term1 = -(math.pi/6)**(1./3)
        term2 = (utilities.del_adiabatic*core_density**(4/3.))/star.core_pressure
        temperature = math.exp(term1*utilities.gravitational_constant*term2*mass**(2./3.)) + star.core_temp
    else:
        if test:
            print "Radiative"
        term1 = -1/(2*utilities.radiation_density_constant*utilities.speed_of_light)
        term2 = (3/(4*math.pi))**(2./3.)
        term3 = (core_opacity*e_n*core_density)**(4./3.)
        temperature = (term1*term2*term3*mass*(2./3.))**(1./4.) + star.core_temp

    radius = (3./(4.*math.pi*core_density))**(1./3.) * mass**(1./3.)

    if test:
        print "Core Density: ", core_density
        print "Core luminosity: ", core_luminosity
        print "Core opacity: ", core_opacity
        print [pressure, temperature, radius, core_luminosity]


    return [pressure, temperature, radius, core_luminosity]

def inward_start(star, test=True):
    """
    Computing the initial guesses at the surface.

    guess outer radius and ltot
    """
    surface_density = 10e-7

    # Make an array of density values
    densities = np.linspace(1e-9,surface_density,1e5)
    # Calculate corresponding pressure arrays
    pressures1 = np.asarray([star.pressure_from_ideal(density) for density in densities])
    # Calculate an array of opacity values for computing pressure2
    opacities = [10**opacity_interpolation.opacity_is(math.log10(star.teff), math.log10(density)) for density in densities]
    pressures2 = np.asarray([star.calc_other_pressure(opacity) for opacity in opacities])

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

    return [surface_pressure, star.teff, star.total_radius, star.total_lum]

def derivatives(layer, enclosed_mass, star, test=True):
    """
    The enclosed_mass given should be the enclosed enclosed_mass, deal with that
    outside the function.
    """
    density = calc_density.density_is(math.log10(layer[1]), math.log10(layer[0]), star.hydrogen_mass, star.helium_mass)
    opacity = 10**opacity_interpolation.opacity_is(math.log10(layer[1]), math.log10(density))

    dpressure = -((utilities.gravitational_constant)/(4*math.pi))*((enclosed_mass)/(layer[2]**4))
    dradius = (1./(4.*math.pi))*(1./(density*layer[2]**2))
    dluminosity = utilities.calc_e_n(star, density, layer[1])

    del_rad = star.calc_del_rad(density, layer[0], layer[1], opacity, layer[3], enclosed_mass)
    if del_rad >= utilities.del_adiabatic:
        if test:
            print "Convective"
        dtemperature = -((utilities.gravitational_constant*enclosed_mass*layer[1])/(4*math.pi*layer[0]*layer[2]**4))*utilities.del_adiabatic
    else:
        if test:
            print "Radiative"
        dtemperature = -((utilities.gravitational_constant*enclosed_mass*layer[1])/(4*math.pi*layer[0]*layer[2]**4))*del_rad

    if test:
        print "Enclosed Mass: ", enclosed_mass
        print "Density: ", density
        print "Opacity: ", opacity
        print [dpressure, dtemperature, dradius, dluminosity], '\n'
    return [dpressure, dtemperature, dradius, dluminosity]

"""
Making testing suite. This now longer makes sense.
"""
if __file__ == sys.argv[0]:
    core =  outward_start(star, mass_initial)
    surface =  inward_start(star)
    print derivatives(star, core, mass_initial)
    print derivatives(star, surface, total_mass)
