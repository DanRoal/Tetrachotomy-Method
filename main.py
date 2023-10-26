"""
Implementación de código

21/10/2023
Autor: DanRoal
"""

import numpy as np
import sympy as sp

#definimos el cambio de variable
x_0, x_1 = 0, 1
y_0, y_1 = 0, 1

t = sp.symbols('t')
x = sp.symbols('x')


def P(x, m=2):
    numerador = sp.integrate((t)**2 * (1-t)**2, (t, 0, x))
    denominador = sp.integrate((t)**2 * (1-t)**2, (t, 0, 1))
    return numerador/denominador

def delta_P(argumento, m=2):
    delta = sp.diff(P(x, m), x)
    funcion =  sp.lambdify(x, delta)
    return funcion(argumento)

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

N=25
h = (1)/N


def integral_abajo(a,b):
    Expr = (t-a)**2 * (b-t)**2
    P_denominador = 1/sp.integrate(Expr, (t, a, b))
    P_denominador = float(P_denominador)
    Expr = sp.lambdify(t, Expr)

    list_gamma = [gamma_down(P(i/N,a,b)) for i in range(N)]
    producto = h*gamma_Delta_down*P_denominador
    
    int_down_0 = sum([f_1(list_gamma[i])* Expr(i/N) for i in range(N)]) *producto
    int_down_1 = sum([f(list_gamma[i])* Expr(i/N) for i in range(N)]) *producto
    int_down_2 = sum([f(list_gamma[i])* Expr(i/N) * P(i/N,a ,b) for i in range(N)]) *producto
    return [int_down_0, int_down_1, int_down_2]

def integral_derecha(a,b):
    Expr = (t-a)**2 * (b-t)**2
    P_denominador = 1/sp.integrate(Expr, (t, a, b))
    P_denominador = float(P_denominador)
    Expr = sp.lambdify(t, Expr)

    list_gamma = [gamma_right(P(i/N,a,b)) for i in range(N)]
    producto = h*gamma_Delta_right*P_denominador

    int_right_0 = sum([f_1(list_gamma[i])* Expr(i/N) for i in range(N)]) *producto
    int_right_1 = sum([f(list_gamma[i])* Expr(i/N) for i in range(N)]) *producto
    int_right_2 = sum([f(list_gamma[i])* Expr(i/N) * P(i/N,a ,b) for i in range(N)]) *producto
    return [int_right_0, int_right_1, int_right_2]

def integral_arriba(a,b):
    Expr = (t-a)**2 * (b-t)**2
    P_denominador = 1/sp.integrate(Expr, (t, a, b))
    P_denominador = float(P_denominador)
    Expr = sp.lambdify(t, Expr)

    list_gamma = [gamma_up(P(i/N,a,b)) for i in range(N)]
    producto = h*gamma_Delta_up*P_denominador

    int_up_0 = sum([f_1(list_gamma[i])* Expr(i/N) for i in range(N)]) * producto
    int_up_1 = sum([f(list_gamma[i])* Expr(i/N) for i in range(N)]) * producto
    int_up_2 = sum([f(list_gamma[i])* Expr(i/N) * P(i/N,a ,b) for i in range(N)]) * producto
    return [int_up_0, int_up_1, int_up_2]

def integral_izquierda(a,b):
    Expr = (t-a)**2 * (b-t)**2
    P_denominador = 1/sp.integrate(Expr, (t, a, b))
    P_denominador = float(P_denominador)
    Expr = sp.lambdify(t, Expr)

    list_gamma = [gamma_left(P(i/N,a,b)) for i in range(N)]
    producto = h*gamma_Delta_left*P_denominador

    int_left_0 = sum([f_1(list_gamma[i])* Expr(i/N) for i in range(N)]) *producto
    int_left_1 = sum([f(list_gamma[i])* Expr(i/N) for i in range(N)]) *producto
    int_left_2 = sum([f(list_gamma[i])* Expr(i/N) * P(i/N,a ,b) for i in range(N)]) *producto
    return [int_left_0, int_left_1, int_left_2]

lista_abajo = integral_abajo(0, np.pi/2)
lista_derecha = integral_derecha(np.pi/2, np.pi)
lista_arriba = integral_arriba(np.pi, 3*np.pi/2)
lista_izquierda = integral_izquierda(3*np.pi/2, 2*np.pi)

integral_0 = lista_abajo[0] + lista_derecha[0] + lista_arriba[0] + lista_izquierda[0]
integral_1 = lista_abajo[1] + lista_derecha[1] + lista_arriba[1] + lista_izquierda[1]
integral_2 = lista_abajo[2] + lista_derecha[2] + lista_arriba[2] + lista_izquierda[2]

print(integral_0)
print(integral_1)
print(integral_2)
print(integral_1/integral_0)
print(integral_2/integral_1)
