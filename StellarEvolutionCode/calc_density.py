'''
PURPOSE:
    Functions to calculate the density of gas given a chemical composition, temperature,
    and pressure.

CALLING SEQUENCE:
    rho = density_is(temperature, pressure, mass_hydrogen, mass_helium)
    Where,
    mass_hydrogen = mass fraction of hydrogen
    mass_helium = mass fraction of helium
'''
import utilities

def mu_is(mass_hydrogen, mass_helium):
    return 2./(1. + 3.*mass_hydrogen + 0.5*mass_helium)

def gas_pressure_is(rho, logT, mass_hydrogen, mass_helium):
    mu = mu_is(mass_hydrogen, mass_helium)
    return r*rho*(10**logT) / mu

def rad_pressure_is(logT):
    return 0.33*utilities.radiation_density_constant*(10*logT)**4

def density_is(logT, logP, mass_hydrogen, mass_helium):
    mu = mu_is(mass_hydrogen, mass_helium)
    top = 10**logP - ((0.33)*utilities.radiation_density_constant*(10**logT)**4)
    bottom = utilities.gas_constant*(10**logT)
    return mu*(top/bottom)

