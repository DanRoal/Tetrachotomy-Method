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


tic_1 = time.time()
lista_P_i = [P(i/N) for i in range(N+1)]
l_DelP = [delta_P(i/N) for i in range(N+1)]

int_0_abajo = sum([f_1(gamma_down(lista_P_i[i]*0.25))*l_DelP[i] for i in range(N)])* gamma_Delta_down
int_0_derecha = sum([f_1(gamma_right(lista_P_i[i]*0.25 + 0.25))*l_DelP[i] for i in range(N)])* gamma_Delta_right
int_0_arriba = sum([f_1(gamma_up(lista_P_i[i]*0.25 + 0.5))*l_DelP[i] for i in range(N)])* gamma_Delta_up
int_0_izquierda = sum([f_1(gamma_left(lista_P_i[i]*.25 + 0.75))*l_DelP[i] for i in range(N)])* gamma_Delta_left

int_1_abajo = sum([f(gamma_down(lista_P_i[i]*0.25))*l_DelP[i] for i in range(N)])* gamma_Delta_down
int_1_derecha = sum([f(gamma_right(lista_P_i[i]*0.25 + 0.25))*l_DelP[i] for i in range(N)])* gamma_Delta_right
int_1_arriba = sum([f(gamma_up(lista_P_i[i]*0.25 + 0.5))*l_DelP[i] for i in range(N)])* gamma_Delta_up
int_1_izquierda = sum([f(gamma_left(lista_P_i[i]*.25 + 0.75))*l_DelP[i] for i in range(N)])* gamma_Delta_left

int_2_abajo = sum([f(gamma_down(lista_P_i[i]*0.25))*l_DelP[i]*gamma_down(lista_P_i[i]*0.25) for i in range(N)])* gamma_Delta_down
int_2_derecha = sum([f(gamma_right(lista_P_i[i]*0.25 + 0.25))*l_DelP[i]*gamma_right(lista_P_i[i]*0.25 + 0.25) for i in range(N)])* gamma_Delta_right
int_2_arriba = sum([f(gamma_up(lista_P_i[i]*0.25 + 0.5))*l_DelP[i]*gamma_up(lista_P_i[i]*0.25 + 0.5) for i in range(N)])* gamma_Delta_up
int_2_izquierda = sum([f(gamma_left(lista_P_i[i]*.25 + 0.75))*l_DelP[i]*gamma_left(lista_P_i[i]*.25 + 0.75) for i in range(N)])* gamma_Delta_left


resultadoaprox = sp.simplify((int_0_abajo + int_0_derecha + int_0_arriba + int_0_izquierda)*h*.25)

toc_1 = time.time()



print(f"Resultado aproximado: {resultadoaprox}")

print(f"tiempo ejecución aproximado: {toc_1-tic_1}")
