import math

"""
Constants:
"""
gravitational_constant     = 6.673e-8      # cm^3 g^-1 s^-2
speed_of_light             = 3e10          # cm s^-1
stefan_boltzmann_constant  = 5.67e-5       # erg cm^-2 s^-1 K^-1
radiation_density_constant = 7.56464e-15   # ergs cm^-2 K^-4
mass_of_hydrogen           = 1.673e-24     # cgs
del_adiabatic              = 0.4
boltzmann_constant         = 1.38e-16      # ergs K^-1
gas_constant               = 8.32e7        # K^-1 mol^-1

"""
Functions:
"""
def nond_temperature(temperature, power):
    """
    A function to non-dimesionlize temperature to be
    consistent with Kippenhahn, Weigert, and Weiss.
    """
    return temperature/(10**power)

def calc_e_n(self, density, temperature):
    t_9 = nond_temperature(temperature, 9)
    if nond_temperature(temperature, 7) >= math.fabs(2 - 0.1):
        psi =  2.0
    else:
        psi = 1.5
    f_11 = 1
    g_11 = 1. + 3.82*t_9 + 1.151*t_9**2 + 0.144*t_9**3 - 0.0114*t_9**4
    pp_e_n = (2.57e4)*psi*f_11*g_11*density*(self.hydrogen_mass**2)*(t_9**(-2./3.))*math.exp(-3.381/(t_9**(1./3.)))

    g_14_1 = 1- 2.0*t_9 + 3.41*t_9**2 - 2.43*t_9**3
    X_cno = 0.7*(1 - (self.hydrogen_mass + self.helium_mass))
    cno_e_n = (8.24e25)*g_14_1*X_cno*self.hydrogen_mass*density*(t_9**(-2./3.))*math.exp(-15.231*(t_9**(-1./3.)) - (t_9/0.8)**2)

    return pp_e_n + cno_e_n


