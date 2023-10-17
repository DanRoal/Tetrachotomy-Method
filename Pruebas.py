"""
Primer acercacimiento al intentar hacer un programa que haga el método de tetracotomía para funciones de
variable compleja.
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad

def x(t, p, x_0):
    return p*(abs(sp.cos(t))*sp.cos(t) + abs(sp.sin(t))*sp.sin(t) + x_0)

def y(t, q, y_0):
    return q*(abs(sp.cos(t))*sp.cos(t) - abs(sp.sin(t))*sp.sin(t) + y_0)

t = sp.Symbol('t')
p = 1
q = 1
x_0 = 0
y_0 = 0

f = lambda t: sp.diff(x(t,p,x_0), t)

print(quad(f, 0, sp.pi/2))