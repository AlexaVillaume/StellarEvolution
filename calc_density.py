'''
PURPOSE:
    Functions to calculate the density of gas given a chemical composition, temperature,
    and pressure.

CALLING SEQUENCE:
    rho = density_is(temperature, pressure, X, Y)
    Where,
    X = abundance of hydrogen
    Y = abundance of helium
'''
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

