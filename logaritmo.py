"""
construcción de función logaritmo a la cual se le pueda modificar el corte de rama

"""

from cmath import log, phase, polar, pi

def logaritmo(z, corte = 0):
    """
    Esta función calcula el logaritmo de un número complejo z, con corte de rama en corte
    """
    r, theta = polar(z)
    return log(r) + (theta + corte*2*pi)*1j

print(log(-1-0.000001j))


