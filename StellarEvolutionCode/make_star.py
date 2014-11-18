"""
This will be the control program for running the stellar evolution code.
"""


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
        return  (self.total_luminosity/(4*math.pi*sigma*self.total_radius**2))**(1./4.)

# All values for the Sun
solar = star(2.526e14, 1.57e7, 3.846e33, 7e10, 1.98e33, 0.70, 0.27)
