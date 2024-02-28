from sympy import Number, symbols

def lagrange(Lx, Ly, decimales):
    """
    dada 2 listas de puntos retorna el polinomio de lagrange en forma algebráica
    """
    x=symbols('x')
    if  len(Lx)!= len(Ly):
        print("ERROR")
        return 1
    polinomio=float(0.0)
    for k in range ( len(Lx) ):
        t=float(1.0)
        for j in range ( len(Lx) ):
            if j != k:
                t=t* ( (x-Lx[j]) / float(Lx[k]-Lx[j])) 
        polinomio+= t*Ly[k]
    polinomio = polinomio.xreplace({n.evalf() : round(n, decimales) for n in polinomio.atoms(Number)})
    return polinomio