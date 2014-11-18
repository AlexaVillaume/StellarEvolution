"""
This will be the control program for running the stellar evolution code.
"""
import math
import utilities

class star(object):
    def __init__(self, core_pressure, core_temp, total_lum, total_radius, total_mass, hydrogen_mass, helium_mass):
        """
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

    def calc_teff(self):
        return  (self.total_lum/(4*math.pi*utilities.stefan_boltzmann_constant*self.total_radius**2))**(1./4.)

    # Need a more descriptive name for this function
    def calc_other_pressure(self, opacity):
        return (2*G*self.total_mass)/(3*opacity*self.total_radius**2)

    def pressure_from_ideal(self, density, teff):
        mu = calc_density.mu_is(self.hydrogen_mass, self.helium_mass)
        return (density*utilities.boltzmann_constant*teff)/(mu*utitlities.mass_h)

    def calc_luminosity(radius, temperature):
        return 4.*math.pi*(radius**2)*utilities.stefan_boltzmann_constant*temperature**4.

    def calc_del_rad(density, pressure, temperature, opacity, luminosity, mass):
        term1 = 3/(16*math.pi*utilities.radiation_density_constant* \
                utilities.speed_of_light*utilities.gravitiational_constnat)
        term2 = (opacity*luminosity*pressure)/(mass*temperature**4)
        return term1*term2

   def calc_e_n(density, temperature):
        return 1

# All values for the Sun
solar = star(2.526e14, 1.57e7, 3.846e33, 7e10, 1.98e33, 0.70, 0.27)
teff = solar.calc_teff()
