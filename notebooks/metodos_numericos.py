# Módulo para guardar los algoritmos de métodos numéricos


def biseccion(a, b, func, tolerancia=0.00001, iteracion=1, resultado=[]):
    """
    'a' y 'b', func(a)* func(b) tiene que ser menor a uno,
    si no no se puede desarrollar el método.
    La 'tolerancia' tiene que ser menor a 0.00001 (el predeterminado)
    Si el valor absoluto de a y b dividido entre 2 es menor a 0 + 'tolerancia'
    , se retorna ese dato si no se vuelve a hacer la operación
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
                "<-- raíz aproximada.",
            ]
        )
        return resultado
    elif func(a) * func(c) < 0:
        resultado.append(
            [str(iteracion), f"{a:,.15f}", f"{b:,.15f}", f"{c:,.15f}", f"{error:,.15f}"]
        )
        return biseccion(a, c, func, tolerancia, iteracion + 1, resultado)
    elif func(b) * func(c) < 0:
        resultado.append(
            [str(iteracion), f"{a:,.15f}", f"{b:,.15f}", f"{c:,.15f}", f"{error:,.15f}"]
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
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", "<-- solución"])
        return resultado
    else:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{error:,.15f}"]
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
    g(x)
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
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", "<-- solución"])
        return resultado
    else:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{error:,.15f}"]
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
    Método de secante acorde a algorimo de clase:
    https://www.udbvirtual.edu.sv/materiales_didacticos/AMN941/clase3.html
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
        resultado.append([str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", "<-- solución"])
        return resultado
    else:
        resultado.append(
            [str(iteracion), f"{x0:,.15f}", f"{x1:,.15f}", f"{x2:,.15f}", f"{error:,.15f}"]
        )
        return metodo_de_secante(
            x0=x1,
            x1=x2,
            f=f,
            tolerancia=tolerancia,
            iteracion=iteracion + 1,
            resultado=resultado,
        )
