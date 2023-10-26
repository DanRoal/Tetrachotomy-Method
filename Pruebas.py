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


########### integramos normalmente ############
N=40000

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

tic = time.time()

print()

#Integral precisa
int_abajo = sp.integrate(f(gamma_down(t))*gamma_Delta_down, (t, 0, 1/4))
int_derecha = sp.integrate(f(gamma_right(t))*gamma_Delta_right, (t, 1/4, 1/2))
int_arriba = sp.integrate(f(gamma_up(t))*gamma_Delta_up, (t, 1/2, 3/4))
int_izquierda = sp.integrate(f(gamma_left(t))*gamma_Delta_left, (t, 3/4, 1))

#nuestro intento con la aproximación
toc = time.time()

tic_1 = time.time()
lista_P_i = [P(i/N) for i in range(N+1)]

int_0_abajo = sum([f(gamma_down(lista_P_i[i]*0.25))*delta_P(i/N) for i in range(N)])* gamma_Delta_down
int_0_derecha = sum([f(gamma_right(lista_P_i[i]*0.25 + 0.25))*delta_P(i/N) for i in range(N)])* gamma_Delta_right
int_0_arriba = sum([f(gamma_up(lista_P_i[i]*0.25 + 0.5))*delta_P(i/N) for i in range(N)])* gamma_Delta_up
int_0_izquierda = sum([f(gamma_left(lista_P_i[i]*.25 + 0.75))*delta_P(i/N) for i in range(N)])* gamma_Delta_left

toc_1 = time.time()
resultado = sp.simplify(int_abajo + int_derecha + int_arriba + int_izquierda)

resultadoaprox = sp.simplify((int_0_abajo + int_0_derecha + int_0_arriba + int_0_izquierda)*h*.25)


print(f"Resultado exacto: {resultado}")
print(f"Resultado aproximado: {resultadoaprox}")

print(f"tiempo ejecución exacto: {toc-tic}")
print(f"tiempo ejecución aproximado: {toc_1-tic_1}")


