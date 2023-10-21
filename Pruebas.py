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
N= 1000

h = (1)/N

#obtenemos los límites de cada subintervalo
lista_x = [i*h for i in range(N+1)]

#hacemos una lista con los valores de la función en los límites de cada subintervalo
lista_y = [f(i) for i in lista_x]

#calculamos la suma de los valores de la función en los límites de cada subintervalo
integral = (lista_y[0]+lista_y[-1])/2 + sum(lista_y[1:-1])

#multiplicamos la suma por el ancho de cada subintervalo
integral = integral * h
print(integral*4)

########### integramos con metodo ############
tic = time.time()
integral_met = sum([delta_P(i/N) * f(P(i/N)) for i in range(N)]) * h
toc = time.time()
print(integral_met*4 , toc-tic)### La convergencia es rapidísima 