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

t = sp.symbols('t')
x = sp.symbols('x')
def P(x, a, b):
    numerador = sp.integrate((t-a)**2 * (b-t)**2, (t, a, x))
    denominador = sp.integrate((t-a)**2 * (b-t)**2, (t, a, b))
    return numerador/denominador

def delta_P(argumento, a, b):
    delta = sp.diff(P(x, a, b), x)
    funcion =  sp.lambdify(x, delta)
    return funcion(argumento)

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
    return z/(z - (0.1 + 0.4j))
def f_1(z):
    return 1/(z - (0.1 + 0.4j))

N=100
h = (1)/N



int_down_0 = sum([delta_P(i/N, 0, sp.pi/2) * f_1(gamma_down(P(i/N, 0, sp.pi/2))) for i in range(N)]) * h
int_right_0 = sum([delta_P(i/N, 0, sp.pi/2) * f_1(gamma_right(P(i/N, 0, sp.pi/2))) for i in range(N)]) * h
int_up_0 = sum([delta_P(i/N, 0, sp.pi/2) * f_1(gamma_up(P(i/N, 0, sp.pi/2))) for i in range(N)]) * h
int_left_0 = sum([delta_P(i/N, 0, sp.pi/2) * f_1(gamma_left(P(i/N, 0, sp.pi/2))) for i in range(N)]) * h

integral_0 = int_down_0 + int_right_0 + int_up_0 + int_left_0

int_down_1 = sum([delta_P(i/N, sp.pi/2,sp.pi) * f(gamma_down(P(i/N,sp.pi/2 , sp.pi))) for i in range(N)]) * h
int_right_1 = sum([delta_P(i/N, sp.pi/2,sp.pi) * f(gamma_right(P(i/N,sp.pi/2 , sp.pi))) for i in range(N)]) * h
int_up_1 = sum([delta_P(i/N, sp.pi/2,sp.pi) * f(gamma_up(P(i/N,sp.pi/2 , sp.pi))) for i in range(N)]) * h
int_left_1 = sum([delta_P(i/N, sp.pi/2,sp.pi) * f(gamma_left(P(i/N,sp.pi/2 , sp.pi))) for i in range(N)]) * h

integral_1 = int_down_1 + int_right_1 + int_up_1 + int_left_1
"""
int_down_2 = sum([delta_P(i/N) * f(gamma_down(P(i/N)))* P(i/N) for i in range(N)]) * h
int_right_2 = sum([delta_P(i/N) * f(gamma_right(P(i/N)))* P(i/N) for i in range(N)]) * h
int_up_2 = sum([delta_P(i/N) * f(gamma_up(P(i/N)))* P(i/N) for i in range(N)]) * h
int_left_2 = sum([delta_P(i/N) * f(gamma_left(P(i/N)))* P(i/N) for i in range(N)]) * h

integral_2 = int_down_2 + int_right_2 + int_up_2 + int_left_2
"""

print(integral_0)
print(integral_1)
print(integral_1/integral_0)
