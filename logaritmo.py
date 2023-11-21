"""
construcción de función logaritmo a la cual se le pueda modificar el corte de rama

"""
import cmath
from cmath import log, polar, sqrt, cos, sin, pi

def custom_log(z, corte = pi, flip = True):
    """
    Esta función calcula el logaritmo de un número complejo z, con corte de rama en corte
    """
    r, fase = polar(z)

    real = log(r)
    
    if flip:
        if fase >= (corte-pi) and fase < (corte+pi):
            imag = fase
        else:
            imag = 2*pi - fase
        return complex(real, imag)
    else:
        if fase > (-corte-pi) and fase <= (-corte+pi):
            imag = fase
        else:
            imag = -2*pi - fase
        return complex(real, imag)

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

print(log(complex(-1.0,0.0)))
print(custom_log(1j, corte = pi/2, flip = False))

