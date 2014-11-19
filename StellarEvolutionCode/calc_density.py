'''
PURPOSE:
    Functions to calculate the density of gas given a chemical composition, temperature,
    and pressure.

CALLING SEQUENCE:
    rho = density_is(temperature, pressure, X, Y)
    Where,
    X = mass fraction of hydrogen
    Y = mass fraction of helium
'''
a = 7.56e-15 # radiation density constant [erg cm^-3 k^-4]
k = 1.38e-16 # Boltzmann consant [erg K^-1]
r = 8.32e7   # Gas constant [erg K^-1 mol^-1]

def mu_is(X, Y):
    return 2./(1. + 3.*X + 0.5*Y)

def gas_pressure_is(rho, logT, X, Y):
    mu = mu_is(X, Y)
    return r*rho*(10**logT) / mu

def rad_pressure_is(logT):
    return 0.33*a*(10*logT)**4

def density_is(logT, logP, X, Y):
    mu = mu_is(X, Y)
    top = 10**logP - ((0.33)*a*(10**logT)**4)
    bottom = r*(10**logT)
    return mu*(top/bottom)

