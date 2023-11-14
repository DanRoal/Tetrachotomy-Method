"""
Primer acercacimiento al intentar hacer un programa que haga el método de tetracotomía para funciones de
variable compleja.
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time
import random
#definimos el cambio de variable
def P(x):
    return x**3 *(6*x**2 - 15*x + 10)

def delta_P(x):
    return 30*(x-1)**2 * x**2

#definimos la función que queremos integrar


########### integramos normalmente ############
N=4000

pi = 3.1415926535897932384626433

h = (1)/N

########### integramos con metodo ############


print(pi)

x_0, x_1 = 0, 1
y_0, y_1 = 0, 1

t = sp.symbols('t')
x = sp.symbols('x')

##Definimos las curvas de integración
def gamma_down(t):
    return (4) * (x_1 - x_0) * t + x_0 + y_0*1j
gamma_Delta_down = 4 * (x_1 - x_0)

def gamma_right(t):
    return x_1+ ((4) *(y_1 - y_0) * (t-1/4) + y_0)*1j
gamma_Delta_right = 4 * (y_1 - y_0)*1j

def gamma_up(t):
    return (4) *(x_1 - x_0) * (3/4 - t) + x_0+ y_1*1j
gamma_Delta_up = -4 * (x_1 - x_0)

def gamma_left(t):
    return x_0+ ((4) * (y_1 - y_0) * (1 - t) + y_0)*1j
gamma_Delta_left = -4 * (y_1 - y_0)*1j

#definimos la función que queremos integrar
def f(z: complex):
    return z/(z - (0.5 + 0.5j))
def f_1(z: complex):
    return 1/(z - (0.5 + 0.5j))

#nuestro intento con la aproximación
def func(array):
    t=len(array)*[0]
    i=0
    for e in array:
        if e:
           t[i]=e
           i+=1
    return t

n = 200000
arr = [random.randint(0,9) for i in range(n)]
tic = time.time()
x = func(arr)
toc = time.time()
print("Elapsed time in seconds: ", toc-tic)