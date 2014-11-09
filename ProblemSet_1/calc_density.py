'''

Problem Set 1, Problem 4.

Develop a simple numerical program to calculate the density, given the pressure,
temperature, and X, Y, and Z. Include radiation pressue and assume the gas is
fully ionized. Evalute the density for:

    (a) X=0.00, Y=0.98, log T = 7.55, log P = 16.85 (cgs)
    (b) X=0.70, Y=0.28, log T = 6.91, log P = 16.87 (cgs)

What is Beta, the ratio of gas pressure to total pressure in each case?

'''

a = 7.56e-15 # radiation density constant [erg cm^-3 k^-4]
k = 1.38e-16 # Boltzmann consant [erg K^-1]
r = 8.32e7   # Gas constant [erg K^-1 mol^-1]

a_X = 0
a_Y = 0.98
a_T = 7.55
a_P = 16.85

b_X = 0.70
b_Y = 0.28
b_T = 6.91
b_P = 16.87

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


print "For part a: "
rho = density_is(a_T, a_P, a_X, a_Y)
print "The density is: ", rho
print "Beta is: ", gas_pressure_is(rho, a_T, a_X, a_Y)/10**a_P
print "For part b: ",
rho = density_is(b_T, b_P, b_X, b_Y)
print "The density is: ", rho
print "Beta is: ", gas_pressure_is(rho, b_T, b_X, b_Y)/10**b_P
