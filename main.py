"""
Implementación de código

21/10/2023
Autor: DanRoal
"""

import numpy as np
import sympy as sp
import time
from numba import njit

#definimos el cambio de variable
x_0, x_1 = 0, 1
y_0, y_1 = 0, 1

t, x, z= sp.symbols('t x z')
a, b = sp.symbols('a b')

def P(x, m=2):
    numerador = sp.integrate((t)**(2*m -2) * (1-t)**(2*m -2), (t, 0, x))
    denominador = sp.integrate((t)**(2*m -2) * (1-t)**(2*m -2), (t, 0, 1))
    return numerador/denominador

def delta_P(argumento, m=2):
    delta = sp.diff(P(x, m), x)
    funcion =  sp.lambdify(x, delta)
    return funcion(argumento)

def Pdefault(x):
    return x**3 *(6*x**2 - 15*x + 10)

def delta_Pdefault(x):
    return 30*(x-1)**2 * x**2

##Definimos las curvas de integración
def gamma_down(t, x_0 = 0, x_1=1, y_0=0, y_1=1):
    return (4) * (x_1 - x_0) * t + x_0 + y_0*1j


def gamma_right(t, x_0 = 0, x_1=1, y_0=0, y_1=1):
    return x_1+ ((4) *(y_1 - y_0) * (t-1/4) + y_0)*1j

def gamma_up(t, x_0 = 0, x_1=1, y_0=0, y_1=1):
    return (4) *(x_1 - x_0) * (3/4 - t) + x_0+ y_1*1j

def gamma_left(t, x_0 = 0, x_1=1, y_0=0, y_1=1):
    return x_0+ ((4) * (y_1 - y_0) * (1 - t) + y_0)*1j

#definimos la función que queremos integrar
def f(z):
    return z/((z - (0.2 + 0.2j))*(z - (0.8 + 0.8j)))
def f_1(z):
    return 1/((z - (0.2 + 0.2j))*(z - (0.8 + 0.8j)))

N=4000
h = (1)/N

lista_P = [Pdefault(i/N) for i in range(N)]
list_delP = [delta_Pdefault(i/N) for i in range(N)]

def integral_abajo(f = f,x_0=0, x_1=1, y_0=0, y_1=1):

    gamma_Delta_down = 4 * (x_1 - x_0)

    L_gamm = [gamma_down(lista_P[i]*0.25, x_0, x_1,y_0,y_1) for i in range(N)]
    lista_f = [f(L_gamm[i]) * list_delP[i] for i in range(N)]
    producto = h*gamma_Delta_down*0.25
    
    int_down_0 = sum([lista_f[i] * 1/L_gamm[i] for i in range(1,N)]) *producto
    int_down_1 = sum([lista_f[i] for i in range(N)]) *producto
    int_down_2 = sum([lista_f[i] * L_gamm[i] for i in range(N)]) *producto
    return [int_down_0, int_down_1, int_down_2]

def integral_derecha(f = f,x_0=0, x_1=1, y_0=0, y_1=1):

    gamma_Delta_right = 4 * (y_1 - y_0)*1j

    L_gamm = [gamma_right(lista_P[i]*0.25 + 0.25, x_0, x_1,y_0,y_1) for i in range(N)]
    lista_f = [f(L_gamm[i]) * list_delP[i] for i in range(N)]
    producto = h*gamma_Delta_right*0.25

    int_right_0 = sum([lista_f[i] * 1/L_gamm[i] for i in range(1,N)]) *producto
    int_right_1 = sum([lista_f[i] for i in range(N)]) *producto
    int_right_2 = sum([lista_f[i] * L_gamm[i] for i in range(N)]) *producto
    return [int_right_0, int_right_1, int_right_2]

def integral_arriba(f= f, x_0=0, x_1=1, y_0=0, y_1=1):

    gamma_Delta_up = -4 * (x_1 - x_0)

    L_gamm = [gamma_up(lista_P[i]*0.25 + 0.5, x_0, x_1,y_0,y_1) for i in range(N)]
    lista_f = [f(L_gamm[i]) * list_delP[i] for i in range(N)]
    producto = h*gamma_Delta_up*0.25

    int_up_0 = sum([lista_f[i] * 1/L_gamm[i]for i in range(1,N)]) * producto
    int_up_1 = sum([lista_f[i] for i in range(N)]) * producto
    int_up_2 = sum([lista_f[i] * L_gamm[i] for i in range(N)]) * producto
    return [int_up_0, int_up_1, int_up_2]

def integral_izquierda(f=f,x_0=0, x_1=1, y_0=0, y_1=1):

    gamma_Delta_left = -4 * (y_1 - y_0)*1j

    L_gamm = [gamma_left(lista_P[i]*0.25 + 0.75, x_0, x_1,y_0,y_1) for i in range(N)]
    lista_f = [f(L_gamm[i]) * list_delP[i] for i in range(N)]
    producto = h*gamma_Delta_left*0.25
    
    int_left_0 = sum([lista_f[i] * 1/L_gamm[i] for i in range(1,N)]) *producto
    int_left_1 = sum([lista_f[i] for i in range(N)]) *producto
    int_left_2 = sum([lista_f[i] * L_gamm[i] for i in range(N)]) *producto
    return [int_left_0, int_left_1, int_left_2]

def integrales(funcion, x_0, x_1, y_0, y_1):
    """
    función que obtiene las tres integrales de la función dada en el rectángulo dado.
    """
    lista_abajo = integral_abajo(funcion, x_0, x_1, y_0, y_1)
    lista_derecha = integral_derecha(funcion, x_0, x_1, y_0, y_1)
    lista_arriba = integral_arriba(funcion, x_0, x_1, y_0, y_1)
    lista_izquierda = integral_izquierda(funcion, x_0, x_1, y_0, y_1)

    int_0 = lista_abajo[0] + lista_derecha[0] + lista_arriba[0] + lista_izquierda[0]
    int_1 = lista_abajo[1] + lista_derecha[1] + lista_arriba[1] + lista_izquierda[1]
    int_2 = lista_abajo[2] + lista_derecha[2] + lista_arriba[2] + lista_izquierda[2]

    return[int_0, int_1, int_2]

def Encontrar_polos(funcion, x_0, x_1, y_0, y_1, tol = 1e-10):
    """
    Función que encuentra polos de una función compleja en un rectángulo dado.
    
    El parámetro tol es la diferencia que deben tener 2 cocientes de 
    integrales para que se considere el mísmo número. En otras palabra es la forma de determinar
    si existe un solo polo en el rectángulo dado.
    """
    ints = integrales(funcion, x_0, x_1, y_0, y_1)

    if abs(ints[0].real) and abs(ints[0].imag) and abs(ints[1].real) and abs(ints[1].imag) <= tol:
        return
    polo_1 = ints[1]/ints[0]
    polo_2 = ints[2]/ints[1]


    if abs(polo_1.real - polo_2.real) <= tol and abs(polo_1.imag - polo_2.imag) <= tol:
        resultado = (ints[1]/ints[0] + ints[2]/ints[1])/2
        return resultado
        
    else:
        cuadro_izq = Encontrar_polos(funcion, x_0, (x_0 + x_1)/2, (y_0 + y_1)/2, y_1, tol)
        cuadro_der = Encontrar_polos(funcion, (x_0 + x_1)/2, x_1, (y_0 + y_1)/2, y_1, tol)
        cuadro_arr = Encontrar_polos(funcion, (x_0 + x_1)/2, x_1, y_0, (y_0 + y_1)/2, tol)
        cuadro_abj = Encontrar_polos(funcion, x_0, (x_0 + x_1)/2, y_0, (y_0 + y_1)/2, tol)

        lista = [cuadro_izq, cuadro_der, cuadro_arr, cuadro_abj]
        filtrada = [i for i in lista if i != None]
        
        if len(filtrada) == 0:
            return
        elif len(filtrada) == 1:
            return filtrada[0]
        elif len(filtrada) == 2:
            return filtrada[0], filtrada[1]
        elif len(filtrada) == 3:
            return filtrada[0], filtrada[1], filtrada[2]
        elif len(filtrada) == 4:
            return filtrada[0], filtrada[1], filtrada[2], filtrada[3]

def amortiguado(t):
    return 1/(t*(-1-2j*(5+2j)) + 8**2)


tic = time.time()
polos = Encontrar_polos(amortiguado, -100, 100, -100, 100)
toc = time.time()

print(polos)
print(f"tiempo de ejecución: {toc-tic}")