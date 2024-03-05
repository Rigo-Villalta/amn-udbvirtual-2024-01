from sympy import symbols, factor, simplify

def lagrange_completo(datos_x, datos_y, valor_a_aproximar, decimales_a_usar):
    """
    método de lagrange completo que retorna una lista con:
    - Las funciones 'L' construidas
    - El polinomio de interpolación de lagrange
    - El valor aproximado de la funcio de el valor a aproximar
    """
    x = symbols("x", positive = True)
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
        numerador = (simplify(factor(numerador)))
        L =  numerador / float(denominador)
        L_formateado = f"{numerador}\n{'-'* len(str(numerador))}\n{' ' * ((len(str(numerador))- len(str(denominador))) // 2)}{denominador}"
        resultado.append([f"El L_{i} es: ", str(L_formateado)])
        resultado.append([" ", " "])
        polinomio += L * datos_y[i]
    resultado.append(["El polinomio de lagrange es: ", str(polinomio)])
    resultado.append([" ", " "])
    valor_aproximado = round(polinomio.evalf(18, subs={x: valor_a_aproximar}), decimales_a_usar)
    print(valor_aproximado)
    resultado.append([f"{valor_a_aproximar} Approx. en  f(x) es: ", str(valor_aproximado)])
    return resultado
            


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
