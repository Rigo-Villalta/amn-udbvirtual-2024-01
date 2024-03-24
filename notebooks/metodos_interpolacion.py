from sympy import symbols, nsimplify
from numpy import zeros, float64, diag, linalg, matmul, set_printoptions

set_printoptions(precision=16)


def lagrange(datos_x, datos_y, valor_a_aproximar):
    """
    método de lagrange completo que retorna una lista con:
    - Las funciones 'L' construidas
    - El polinomio de interpolación de lagrange
    - El valor aproximado de la función de el valor a aproximar
    """
    x = symbols("x")
    n = len(datos_x)
    resultado = [["Concepto", "Dato"]]
    polinomio = float(0.0)
    for i in range(n):
        numerador = float(1.0)
        denominador = float(1.0)
        for j in range(n):
            if j != i:
                numerador = numerador * (x - datos_x[j])
                denominador = denominador * (datos_x[i] - datos_x[j])
        L = numerador / denominador
        L_formateado = (
            f"{nsimplify(numerador, tolerance=0.00001, rational=True)}\n"
            f"{'-'* len(str(numerador))}\n"
            f"{' ' * ((len(str(numerador))- len(str(denominador))) // 2)}{float(denominador)}"
        )
        evaluacion_en_l_i = float(L.evalf(subs={x: valor_a_aproximar}))
        resultado.append(
            [
                f"La evaluación de {valor_a_aproximar} en L_{i} es: ",
                f"{str(evaluacion_en_l_i)}\n",
            ]
        )
        resultado.append([f"El L_{i} es: ", f"{str(L_formateado)}\n"])
        polinomio += L * datos_y[i]
    resultado.append(["El polinomio de Lagrange es: ", str(polinomio)])
    resultado.append([" ", " "])
    valor_aproximado = polinomio.evalf(subs={x: valor_a_aproximar})
    resultado.append(
        [f"{valor_a_aproximar} Approximado en  f(x) es: ", str(valor_aproximado)]
    )
    return resultado


def neville(datos_x, datos_y, x):
    """
    Se retorna una array de numpy con forma de matriz
    para cumplir los criterios de la interpolación de neville.
    El útlimo valor 'matriz[n][n]' es el valor aproximado de la función
    """
    n = len(datos_x)
    matriz = zeros((n, n), dtype=float64)
    matriz[:, 0] = datos_y
    for i in range(1, n):
        for j in range(1, i + 1):
            matriz[i][j] = (
                ((x - datos_x[i - j]) * matriz[i][j - 1])
                - ((x - datos_x[i]) * matriz[i - 1][j - 1])
            ) / (datos_x[i] - datos_x[i - j])
    return matriz


def diferencias_divididas(datos_x, datos_y, valor_a_aproximar):
    """
    método de diferencias divididas que retorna una lista con:
    - El valor aproximado de la función de el valor a aproximar
    """
    n = len(datos_x)
    matriz = zeros((n, n), dtype="f")
    matriz[:, 0] = datos_y
    for i in range(1, n):
        for j in range(1, i + 1):
            matriz[i][j] = (matriz[i][j - 1] - matriz[i - 1][j - 1]) / (
                datos_x[i] - datos_x[i - j]
            )
    coeficientes = diag(matriz)
    x = symbols("x")
    polinomio = coeficientes[0]
    for i in range(1, n):
        termino = coeficientes[i]
        for k in range(0, i):
            termino = termino * (x - datos_x[k])
        polinomio = polinomio + termino
    valor_aproximado = polinomio.evalf(subs={x: valor_a_aproximar})
    return (matriz, polinomio, valor_aproximado)


def hermite(datos_x, datos_y, datos_y_prima, valor_a_aproximar):
    """
    método de hermite, par este método necesitados suministrar como argumentos:
    - Una lista de valores datos_x que son los valores que se pasan a una función
    - Una lista de valores datos_y que son lo valores de datos_x evaluados en la función
    - una lista de valores datos_y_prima que son los valores de datos_x evaluados en la
      primera derivada de la función
    - El valor a aproximar

    Esta función retorna una tupla de tres elementos con:
    - La matriz que contienen los polinomios de interpolación de Hermite.
    - El polinomio interpolador
    - El valor aproximado que resulta de evaluar el valor a aproximar en el polinomio interpolador
    """
    n = len(datos_x)
    matriz = zeros((n * 2, n * 2), dtype=float64)
    z = []
    for dato in datos_x:
        z.append(dato)
        z.append(dato)
    primera_columna = []
    for dato in datos_y:
        primera_columna.append(dato)
        primera_columna.append(dato)
    matriz[:, 0] = primera_columna
    for i in range(1, n * 2):
        j = 1
        if i % 2 == 0:
            matriz[i][j] = (matriz[i][j - 1] - matriz[i - 1][j - 1]) / (
                z[i] - z[i - 1]
            )
        else:
            matriz[i][j] = datos_y_prima[i // 2]

    for i in range(2, n * 2):
        for j in range(2, i + 1):
            matriz[i][j] = round((matriz[i][j - 1] - matriz[i - 1][j - 1]) / (z[i] - z[i - j]), 16)
    x = symbols("x")
    polinomio = matriz[0][0]
    for i in range(1, n * 2 -1):
        termino = matriz[i][i]
        for k in range(0, i):
            termino = termino * (x - z[k])
        polinomio = polinomio + termino
    valor_aproximado = polinomio.evalf(subs={x: valor_a_aproximar})
    return (matriz, polinomio, valor_aproximado)


def trazador_cubico_natural(datos_x, datos_y):
    """
    Nos proporcionan un conjunto de datos datos_x, datos_y y a partir
    de estos se construye el trazador cúbico con n polinomios. Esta función
    retorna una tupla de 3 elementos con:
    - La matriz de trazadores
    - Los trazadores cúbicos para cada tramo.
    - 
    """
    n = len(datos_x)
    x = symbols("x")

    # calculo la lista de datos h, es decir la distancia entre nodos x
    h = []
    for i in range(0, n - 1):
        h.append(datos_x[i + 1] - datos_x[i])

    # Declaro la matriz de coeficientes, los lleno de ceros
    # y coloco los unos en los extremos de la diagonal
    coeficientes = zeros((n, n), dtype=float64)
    coeficientes[0][0] = 1
    coeficientes[n - 1][n - 1] = 1
    # LLenamos el resto de la matriz acorde a las fórmulas
    for i in range(1, n - 1):
        for j in range(i - 1, i + 2):
            if j == i - 1:
                coeficientes[i][j] = h[i - 1]
            elif j == i:
                coeficientes[i][j] = 2 * (h[i - 1] + h[i])
            elif j == i + 1:
                coeficientes[i][j] = h[i]

    # Ahora hago los términos independientes
    TI = zeros((n, 1), dtype=float64)
    for i in range(0, n):
        if i == 0 or i == n - 1:
            TI[i][0] = 0
        else:
            TI[i][0] = 3 * (datos_y[i + 1] - datos_y[i]) - 3 * (
                datos_y[i] - datos_y[i - 1]
            )

    # Calculo los c
    c = matmul(linalg.inv(coeficientes), TI)

    # TODO: Calcular los b y los d
    b = []
    d = []
    for i in range(0, n - 1):
        b.append(
            ((datos_y[i + 1] - datos_y[i]) / h[i])
            - ((h[i] * (2 * c[i] + c[i + 1])) / 3)
        )

    for i in range(0, n - 1):
        d.append((c[i + 1] - c[i]) / (3 * h[i]))

    # Generamos el trazador cúbico como expresión simbólica
    x = symbols("x")
    polinomio_trazador = []
    for i in range(0, n - 1):
        polinomio_i = (
            datos_y[i]
            + b[i] * (x - datos_x[i])
            + float(c[i][0]) * (x - datos_x[i]) ** 2
            + d[i] * (x - datos_x[i]) ** 3
        )
        polinomio_i = polinomio_i[0]
        polinomio_trazador.append(polinomio_i)
    return (coeficientes, TI, polinomio_trazador)
