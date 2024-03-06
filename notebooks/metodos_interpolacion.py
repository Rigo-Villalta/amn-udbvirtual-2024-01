from sympy import symbols, nsimplify, Number
from numpy import zeros, float128, diag
import numpy

numpy.set_printoptions(precision=16)


def lagrange(datos_x, datos_y, valor_a_aproximar, decimales_a_usar):
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
    polinomio_a_mostrar = polinomio.xreplace(
        {n.evalf: round(n, decimales_a_usar) for n in polinomio.atoms(Number)}
    )
    resultado.append(["El polinomio de Lagrange es: ", str(polinomio_a_mostrar)])
    resultado.append([" ", " "])
    valor_aproximado = round(
        polinomio.evalf(subs={x: valor_a_aproximar}), decimales_a_usar
    )
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
    matriz = zeros((n, n), dtype=float128)
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
    método de lagrange completo que retorna una lista con:
    - Las funciones 'L' construidas
    - El polinomio de interpolación de lagrange
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


'''
Esta es la prueba contra la clase del vídeo:
https://www.udbvirtual.edu.sv/materiales_didacticos/AMN941/clase6.html
matriz = diferencias_divididas([0, 1, 2, 3, 4], [0.0, 0.75, 2.25, 3.0, 2.25], 1.5)
'''
