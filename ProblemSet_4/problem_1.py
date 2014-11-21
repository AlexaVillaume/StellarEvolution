"""
Calculating the power-law t_6erature dependence and Gamow peak energies for
thermonuclear reactions.
"""
import sys
import math

def w_is(a_j, a_k, z_j, z_k):
    term = (a_j*a_k)/(a_j + a_k)
    return (z_j**2)*(z_k**2)*term

# t is the dimensionless t_6erature and n is the power it's in ref. to
def t_6_is(t, n):
    return t*10**n

# t_6 in K
def t_is(t_6):
    return t_6/10**7

def gamow_peak_is(w, t):
    return 5.665*(w**(1./3.))*(t**(2./3.))

# nu is the power in equation 18.36
def nu_is(w, t):
    return (6.574*(w**(1./3.))*(t**(-1./3))) - 2./3.

hydrogen_a = 1.
hydrogen_z = 1.
be_7_a = 7.
be_7_z = 4.
n_14_a = 14.
n_14_z = 7.

# for H + H
w_a = w_is(hydrogen_a, hydrogen_a, hydrogen_z, hydrogen_z)
t_6 = 10.
t_a = t_is(t_6_is(t_6, 6))
print "Part a: "
print "The Gamow peak is: ", gamow_peak_is(w_a, t_a)
print "The power-law temperature dependence is: ", nu_is(w_a, t_a)
print '\n'

# for 7Be + H
w_b = w_is(be_7_a, hydrogen_a, be_7_z, hydrogen_z)
t_6 = 15.
t_b = t_is(t_6_is(t_6, 6))
print "Part b: "
print "The Gamow peak is: ", gamow_peak_is(w_b, t_a)
print "The power-law temperature dependence is: ", nu_is(w_b, t_b)
print '\n'

# for 14N + H
w_c = w_is(n_14_a, hydrogen_a, n_14_z, hydrogen_z)
t_6 = 15.
t_c = t_is(t_6_is(t_6, 6))
print "Part c: "
print "The Gamow peak is: ", gamow_peak_is(w_c, t_a)
print "The power-law temperature dependence is: ", nu_is(w_c, t_c)
print '\n'


