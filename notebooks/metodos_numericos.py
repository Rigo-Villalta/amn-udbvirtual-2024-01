# Módulo para guardar los algoritmos de métodos numéricos

from math import sqrt

def biseccion(a, b, func, tolerancia=0.00001, iteracion=1, resultado=[]):
    """
    El método de bisección toma un intervalo [a, b] y retorna una raíz aproximada
    'a' y 'b', func(a)* func(b) tiene que ser menor a uno, si no no se puede desarrollar el método.
    El ultimo valor retornado en 'c' es la raíz aproximada.
    """
    if func(a) * func(b) >= 0:
        raise ValueError("Los argumentos no tienen signos diferentes o no son enteros")
    if not callable(func):
        raise TypeError("No ha suministrado una función")
    if tolerancia > 0.00001:
        raise ValueError("El toleracia no puede ser mayor a 0.00001")

    # Para tabular vamos a ingresar los datos en una lista multidimensional
    # la primera columna es para los encabezados de las tablas.

    c = (a + b) / 2
    error = abs(func(c)) if iteracion == 1 else abs(round((c - b), 15))
    if iteracion == 1:
        resultado.append(["# de iteración", "a", "b", "c", "error"])

    if abs(error) <= (0 + tolerancia):
        resultado.append(
            [
                str(iteracion),
                f"{a:,.15f}",
                f"{b:,.15f}",
                f"{c:,.15f}",
                f"{error:,.5E}"
            ]
        )
        return resultado
    elif func(a) * func(c) < 0:
        resultado.append(
            [str(iteracion), f"{a:,.15f}", f"{b:,.15f}", f"{c:,.15f}", f"{error:,.5E}"]
        )
        return biseccion(a, c, func, tolerancia, iteracion + 1, resultado)
    elif func(b) * func(c) < 0:
        resultado.append(
            [str(iteracion), f"{a:,.15f}", f"{b:,.15f}", f"{c:,.15f}", f"{error:,.5E}"]
        )
        return biseccion(b, c, func, tolerancia, iteracion + 1, resultado)
    else:
        raise ValueError("El valor procesado no es correcto")


def newton_raphson(
    x0, funcion, derivada, tolerancia=0.00001, iteracion=1, resultado=[]
):
    """
    Método de newton raphson, retorna una lista de listas para generar una tabla
    de rich con la utilidad indicada.
    Se asume que se ha probado que xo es suficientemente cercano a la raíz para
    poder encontrarla, se establece un máximo de iteraciones.
    La derivada debe ser suministrada.
    El último x1 es la raíz y por tanto el valor de aproximación buscado
    """
    if iteracion == 1:
        resultado.append(["# de iteración", "x0", "x1", "error"])
    elif iteracion == 100:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", "Máximas iteraciones", "posibles"]
        )
        return resultado
    x1 = x0 - (funcion(x0) / derivada(x0))
    error = abs(x1 - x0)
    if error < tolerancia:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{error:,.5E}"])
        return resultado
    else:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{error:,.5E}"]
        )
        return newton_raphson(
            x0=x1,
            funcion=funcion,
            derivada=derivada,
            tolerancia=tolerancia,
            iteracion=iteracion + 1,
            resultado=resultado,
        )


def punto_fijo(x0, g_de_x, tolerancia=0.00001, iteracion=1, resultado=[]):
    """
    Método de newton raphson, retorna una lista de listas para generar una tabla
    de rich con la utilidad indicada.
    Se asume que se tiene un x0 cercano a la raíz y que se proporciona la función
    g(x) la cual es una función de la form x = f(x).
    El último x1 es la raiz aproximada.
    """
    if iteracion == 1:
        resultado.append(["# de iteración", "x0", "x1", "error"])
    elif iteracion == 100:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", "Máximas iteraciones", "posibles"]
        )
        return resultado
    x1 = g_de_x(x0)
    error = abs(x1 - x0)
    if error < tolerancia:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{error:,.5E}"])
        return resultado
    else:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{error:,.5E}"]
        )
        return punto_fijo(
            x0=x1,
            g_de_x=g_de_x,
            tolerancia=tolerancia,
            iteracion=iteracion + 1,
            resultado=resultado,
        )


def metodo_de_secante(x0, x1, f, tolerancia=0.00001, iteracion=1, resultado=[]):
    """
    Método de secante. Se solicitan dos puntos x0 y x1 y que la raíz esté entre ellos y
    la función f. El útimo x1 mostrado es la raíz aproximada encontrada.
    """
    if iteracion == 1:
        resultado.append(["# de iteración", "x0", "x1", "x2", "error"])
    elif iteracion == 100:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", "Máximas iteraciones", "posibles"]
        )
        return resultado

    x2 = x1 - ((f(x1) * (x1 - x0)) / (f(x1) - f(x0)))
    error = abs(x2 - x1)
    if error < tolerancia:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{error:,.5E}"])
        return resultado
    else:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{error:,.5E}"]
        )
        return metodo_de_secante(
            x0=x1,
            x1=x2,
            f=f,
            tolerancia=tolerancia,
            iteracion=iteracion + 1,
            resultado=resultado,
        )


def posicion_falsa(x0, x1, f, tolerancia=0.00001, iteracion=1, resultado=[]):
    """
    Método de posición falsaSe solicitan dos puntos x0 y x1 y que la raíz esté entre ellos y
    la función f. El último x1 es la raíz aproximada.
    """

    x2 = x1 - ( (f(x1) * (x1 - x0)) / (f(x1) - f(x0)))

    if iteracion == 1:
        resultado.append(["# de iteración", "x0", "x1", "x2", "error"])
        error= abs(x2 -x1)
    elif iteracion == 100:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", "Máximas iteraciones", "posibles"]
        )
        return resultado
    else:
        error = abs(x2 - x0)

    if error < tolerancia:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", "<-- solución"])
        return resultado
    else:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{error:,.5E}"])
        if f(x0)*f(x2) < 0:
            return posicion_falsa(
               x2, x0, f, tolerancia, iteracion+1, resultado   
            )
        else:
            return posicion_falsa(
               x2, x1, f, tolerancia, iteracion+1, resultado
            )


def steffensen(x0, g, toleracia=0.00001, iteracion=1, resultado=[]):
    """
    la función g suministrada debe ser en forma g(x) = x desde la función original
    """
    if iteracion == 1:
        resultado.append(["# de iteración", "x0", "x1", "x2", "x3", "error"])
    elif iteracion == 100:
        resultado.append([str(iteracion), "máximas", " iteraciones",  "posibles", "--", "--"])
        return resultado
        
    
    x1 = g(x0)
    x2 = g(x1)
    x3 = x0 - (x1 - x0)**2 / (x2 - 2*x1 + x0)
    error = abs(x3 - x0)
    if error < toleracia:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{x3:,.15f}", f"{error:,.5E}"])
        return resultado
    else:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{x3:,.15f}", f"{error:,.5E}"])
        return steffensen(
            x3,
            g,
            toleracia,
            iteracion + 1,
            resultado
        )


def muller(x0, x1, x2, funcion, toleracia=0.00001, iteracion=1, resultado=[]):
    """
    Método de muller: Para este método tomo 3 estimaciones iniciales y una función
    f(x) = 0.
    """
    if iteracion == 1:
        resultado.append(["# de iteración", "x0", "x1", "x2", "x3", "error"])
    elif iteracion == 100:
        resultado.append([str(iteracion), "máximas", " iteraciones",  "posibles", "--", "--"])
        return resultado

    f0 = funcion(x0)
    f1 = funcion(x1)
    f2 = funcion(x2)
    a = (((f0-f2)*(x1-x2)) - ((f1 - f2)*(x0 - x2))) / ((x0-x2)*(x1-x2)*(x0-x1))
    b = (((f1-f2)*((x0-x2)**2)) - ((f0 - f2)*((x1 - x2)**2))) / ((x0-x2)*(x1-x2)*(x0-x1))
    c=f2
    if b > 0:
        x3 = x2 - ((2*c) / (b + sqrt(b**2 - 4*a*c)))
    else:
        x3 = x2 - ((2*c) / (b - sqrt(b**2 - 4*a*c)))
    error = abs(x3 - x2)
    if error < toleracia:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{x3:,.15f}", f"{error:,.5E}"])
        return resultado
    else:
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{x3:,.15f}", f"{error:,.5E}"])
        return muller(x1, x2, x3, funcion, toleracia, iteracion + 1, resultado)
