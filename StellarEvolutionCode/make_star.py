"""
This will be the control program for running the stellar evolution code.
"""
import sys
import math
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
import calc_density
import shootf
import utilities



class star(object):
    def __init__(self, core_pressure, core_temp, total_lum, total_radius, total_mass, hydrogen_mass, helium_mass):
        """
        Things we know a priori about the star.

        Units:
        Core Pressure:    g cm^-2
        Core Temperature: K
        Total Luminosity: erg s^-1
        Total Radius:     cm
        Total Mass:       g
        """
        self.core_pressure = core_pressure
        self.core_temp = core_temp
        self.total_lum = total_lum
        self.total_radius = total_radius
        self.total_mass = total_mass
        self.hydrogen_mass = hydrogen_mass
        self.helium_mass = helium_mass
        self.teff = None

    def calc_teff(self):
        return  (self.total_lum/(4*math.pi*utilities.stefan_boltzmann_constant*self.total_radius**2))**(1./4.)

    # Need a more descriptive name for this function
    def calc_other_pressure(self, opacity):
        return (2*utilities.gravitational_constant*self.total_mass)/(3*opacity*self.total_radius**2)

    def pressure_from_ideal(self, density):
        mu = calc_density.mu_is(self.hydrogen_mass, self.helium_mass)
        return (density*utilities.boltzmann_constant*self.teff)/(mu*utilities.mass_of_hydrogen)

    def calc_luminosity(self, radius, temperature):
        return 4.*math.pi*(radius**2)*utilities.stefan_boltzmann_constant*temperature**4.

    def calc_del_rad(self, density, pressure, temperature, opacity, luminosity, mass):
        term1 = 3/(16*math.pi*utilities.radiation_density_constant* \
                utilities.speed_of_light*utilities.gravitational_constant)
        term2 = (opacity*luminosity*pressure)/(mass*temperature**4)
        return term1*term2

def generate_inputs():
    return 0

def get_structure(star, inner_masses, outer_masses, mass_step):
    """
    Need to run the re-do the integration until convergence is reached
    i.e. when the difference is zero.

    After the first integration, use Newton's Method to generate new
    initial guesses to give to shootf.
    """
    # Get outward and inward initial conditions
    outward_initial =  shootf.outward_start(star, mass_initial)
    inward_initial =  shootf.inward_start(star)

    difference = 1
    while difference > 0:
        out_inv_jac =  utilities.compute_jacobian(outward_initial, mass_initial, star.core_pressure)
        in_inv_jac =  utilities.compute_jacobian(inward_initial, mass_initial, star.core_pressure)
        outward_initial = outward_initial - np.dot(out_inv_jac, shootf.derivatives(outward_initial, mass_initial, star))
        inward_initial = inward_initial - np.dot(in_inv_jac, shootf.derivatives(inward_initial, star.total_mass, star))
        difference = shootf.integrate(star, inner_masses, outer_masses, mass_step, outward_initial, inward_initial)
        print difference
        sys.exit()

# All values for a two solar mass star
solar_2x = star(1.6032636e17, 20.47409576e6, 15.51844053*(3.846e33), 1.66086519*(7e10),  2*(1.98e33), 0.70, 0.28)
solar_2x.teff = solar_2x.calc_teff()
mass_step = 1e-5 * solar_2x.total_mass

fitting_point = solar_2x.total_mass*0.5
shootf.inner_masses = np.linspace(mass_step, fitting_point, num=100)

#99-100% of mass very fine steps, 1e-8
blah1 = solar_2x.total_mass*0.99
blah2 = solar_2x.total_mass*0.80
outer_masses_1 = np.logspace(math.log10(solar_2x.total_mass),  math.log10(blah1), 1000)
outer_masses_2 = np.logspace(math.log10(blah1),  math.log10(blah2), 1000)
outer_masses_3 = np.logspace(math.log10(blah2),  math.log10(fitting_point), 1000)
outer_masses = np.concatenate((outer_masses_1, outer_masses_2, outer_masses_3), axis=0)

get_structure(solar_2x, inner_masses, outer_masses, mass_step)


# All values for the Sun
solar = star(2.526e14, 1.57e7, 3.846e33, 7e10, 1.98e33, 0.70, 0.28)
solar.teff = solar.calc_teff()
mass_step = 1e-5 * solar.total_mass

fitting_point = solar.total_mass*0.5
inner_masses = np.linspace(mass_step, fitting_point, num=100)

#99-100% of mass very fine steps, 1e-8
blah1 = solar.total_mass*0.99
blah2 = solar.total_mass*0.80
outer_masses_1 = np.logspace(math.log10(solar.total_mass),  math.log10(blah1), 1000)
outer_masses_2 = np.logspace(math.log10(blah1),  math.log10(blah2), 1000)
outer_masses_3 = np.logspace(math.log10(blah2),  math.log10(fitting_point), 1000)
outer_masses = np.concatenate((outer_masses_1, outer_masses_2, outer_masses_3), axis=0)

#get_structure(solar, inner_masses, outer_masses, mass_step)

