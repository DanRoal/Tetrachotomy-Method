"""
Implementación de código

21/10/2023
Autor: DanRoal
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

#definimos el cambio de variable
x_0, x_1 = 0, 1
y_0, y_1 = 0, 1
def P(x):
    return x**3 *(6*x**2 - 15*x + 10)

def delta_P(x):
    return 30*(x-1)**2 * x**2

def gamma_down(t):
    return complex((x_1 - x_0) * t + x_0, y_0)

def gamma_right(t):
    return complex(x_1, (y_1 - y_0) * t + y_0)

def gamma_up(t):
    return complex((x_1 - x_0) * (1 - t) + x_0, y_1)

def gamma_left(t):
    return complex(x_0, (y_1 - y_0) * (1 - t) + y_0)

#definimos la función que queremos integrar
def f(z):
    return 1/((0.5 + 0.5j)-z)

N=50
h = (1)/N

"""
int_down_0 = sum([delta_P(i/N) * f(gamma_down(P(i/N)))* 1/(P(i/N)) for i in range(N)]) * h
int_right_0 = sum([delta_P(i/N) * f(gamma_right(P(i/N)))* 1/(P(i/N)) for i in range(N)]) * h
int_up_0 = sum([delta_P(i/N) * f(gamma_up(P(i/N)))* 1/(P(i/N)) for i in range(N)]) * h
int_left_0 = sum([delta_P(i/N) * f(gamma_left(P(i/N)))* 1/(P(i/N)) for i in range(N)]) * h

inttegral_0 = int_down_0 + int_right_0 + int_up_0 + int_left_0
"""
int_down_1 = sum([delta_P(i/N) * f(gamma_down(P(i/N))) for i in range(N)]) * h
int_right_1 = sum([delta_P(i/N) * f(gamma_right(P(i/N))) for i in range(N)]) * h
int_up_1 = sum([delta_P(i/N) * f(gamma_up(P(i/N))) for i in range(N)]) * h
int_left_1 = sum([delta_P(i/N) * f(gamma_left(P(i/N))) for i in range(N)]) * h

integral_1 = int_down_1 + int_right_1 + int_up_1 + int_left_1

int_down_2 = sum([delta_P(i/N) * f(gamma_down(P(i/N)))* P(i/N) for i in range(N)]) * h
int_right_2 = sum([delta_P(i/N) * f(gamma_right(P(i/N)))* P(i/N) for i in range(N)]) * h
int_up_2 = sum([delta_P(i/N) * f(gamma_up(P(i/N)))* P(i/N) for i in range(N)]) * h
int_left_2 = sum([delta_P(i/N) * f(gamma_left(P(i/N)))* P(i/N) for i in range(N)]) * h

integral_2 = int_down_2 + int_right_2 + int_up_2 + int_left_2

print(integral_1, integral_2)