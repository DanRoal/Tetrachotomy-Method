"""
construcción de función logaritmo a la cual se le pueda modificar el corte de rama

"""
import cmath
from cmath import log, polar, sqrt, cos, sin, pi

def custom_log(z, corte = pi, flip = False):
    """
    Esta función calcula el logaritmo de un número complejo z, con corte de rama en corte
    """
    r, fase = polar(z)
    diferencia = pi - corte
    im = z.imag
    re = z.real
    if flip:
        return log(complex(re,-im)) + (fase - diferencia)*1j
    return log(complex(re,im)) - (fase - diferencia)*1j

def custom_sqrt(z, corte = pi, flip= False):
    """
    Calcula la raíz cuadrada de un número complejo 
    """
    r, fase = polar(z)
    diferencia = pi - corte
    im = z.imag
    re = z.real
    if flip:
        return sqrt(complex(re,-im)) - diferencia
    return sqrt(r)*(cos((fase -diferencia)/2) + sin((fase - diferencia)/2)*1j)

numero = cmath.rect(-1,(0)*pi)

ima = numero.imag
print(sqrt(-1.0-0.0j))
print(log(complex(-1.0,0.0)))
print(custom_log(1j, corte = pi/2, flip = False))
