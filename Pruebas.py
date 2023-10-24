"""
Primer acercacimiento al intentar hacer un programa que haga el método de tetracotomía para funciones de
variable compleja.
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time
#definimos el cambio de variable
def P(x):
    return x**3 *(6*x**2 - 15*x + 10)

def delta_P(x):
    return 30*(x-1)**2 * x**2

#definimos la función que queremos integrar
def f(x):
    return np.sqrt(1 - x**2)

########### integramos normalmente ############
N=4000000

pi = 3.1415926535897932384626433

h = (1)/N

########### integramos con metodo ############
tic = time.time()
integral_met = sum([delta_P(i/N) * f(P(i/N)) for i in range(N)]) * h
toc = time.time()
print(integral_met*4 , toc-tic)### La convergencia es rapidísima 
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

print()

#Integral precisa
int_abajo = sp.integrate(f(gamma_down(t))*gamma_Delta_down, (t, 0, 1/4))
int_derecha = sp.integrate(f(gamma_right(t))*gamma_Delta_right, (t, 1/4, 1/2))
int_arriba = sp.integrate(f(gamma_up(t))*gamma_Delta_up, (t, 1/2, 3/4))
int_izquierda = sp.integrate(f(gamma_left(t))*gamma_Delta_left, (t, 3/4, 1))

#nuestro intento con la aproximación

expr = sp.lambdify(t, (t)**2 * (1-t)**2)
producto = (1/sp.integrate(expr(t), (t, 0, 1))) * h

int_0_abajo = sum([f(gamma_down(P(i/N)))*delta_P(i/N) for i in range(int(N/4))]) *h* gamma_Delta_down
int_0_derecha = sum([f(gamma_right(P(i/N)))*delta_P(i/N) for i in range(int(N/4),int(N/2))]) *h* gamma_Delta_right
int_0_arriba = sum([f(gamma_up(P(i/N)))*delta_P(i/N) for i in range(int(N/2),int(3*N/4))]) *h* gamma_Delta_up
int_0_izquierda = sum([f(gamma_left(P(i/N)))*delta_P(i/N) for i in range(int(3*N/4),N)]) *h* gamma_Delta_left

resultado = sp.simplify(int_abajo + int_derecha + int_arriba + int_izquierda)

resultadoaprox = sp.simplify(int_0_abajo + int_0_derecha + int_0_arriba + int_0_izquierda)

print(resultado)
print(resultadoaprox)
