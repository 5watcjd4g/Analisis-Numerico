import numpy as np


def f(x):
    return np.e**(np.sin(x)-np.cos(x**2))


#Donde f es la funcion, x es donde se evalúa, grado equivale a grado de la función
def taylor(f, x, grado):
    grd = 0
    res = 0

    for grd in range (grado + 1):
        res += (x**grd)/factorial(grd)
        res1 = round(res,4)
        print(str(grd) + " Aproximación: " + str(res1))
    
    return res

def factorial(num):
    aux = num
    result = 1
    while aux != 1 and aux != 0:
        result *= aux
        aux -= 1

    return result

taylor(f, 2, 2)