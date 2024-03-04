from sympy import Number, symbols


def lagrange(datos_x, datos_y, decimales):
    """
    dada 2 listas de puntos retorna el polinomio de lagrange en forma algebráica
    y con un redondeo dado por la cantidad de 'decimales'
    """
    x = symbols("x")
    if len(datos_x) != len(datos_y):
        print("ERROR")
        return 1
    polinomio = float(0.0)
    for k in range(len(datos_x)):
        t = float(1.0)
        for j in range(len(datos_x)):
            if j != k:
                t = t * ((x - datos_x[j]) / float(datos_x[k] - datos_x[j]))
        polinomio += t * datos_y[k]
    polinomio = polinomio.xreplace(
        {n.evalf(): round(n, decimales) for n in polinomio.atoms(Number)}
    )
    return polinomio


def neville(datos_x, datos_y, x):
    """
    Finds an interpolated value using Neville's algorithm.

    Input
      datax: input x's in a list of size n
      datay: input y's in a list of size n
      x: the x value used for interpolation

    Output
      p[0]: the polynomial of degree n
    """
    n = len(datos_x)
    p = n * [0]
    for k in range(n):
        for i in range(n - k):
            if k == 0:
                p[i] = datos_y[i]
            else:
                p[i] = ((x - datos_x[i + k]) * p[i] + (datos_x[i] - x) * p[i + 1]) / (
                    datos_x[i] - datos_x[i + k]
                )
    return p[0]
