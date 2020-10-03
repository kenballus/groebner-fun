# These are Python 2-style cmp functions. You can wrap them with functools.cmp_to_key to use them with Python 3's sorting.
def lex(m1, m2):
    for sym in sorted({**m1.variables, **m2.variables}, key=lambda var: var.symbol):
        if sym in m1.variables and sym in m2.variables:
            if m1.variables[sym] > m2.variables[sym]:
                return 1
            elif m1.variables[sym] < m2.variables[sym]:
                return -1
        elif sym in m1.variables:
            return 1
        elif sym in m2.variables:
            return -1

    return 0

def deglex(m1, m2):
    if m1.degree > m2.degree:
        return 1
    elif m1.degree < m2.degree:
        return -1
    else:
        return lex(m1, m2)

def grevlex(m1, m2):
    if m1.degree > m2.degree:
        return 1
    elif m1.degree < m2.degree:
        return -1
    else:
        for sym in sorted({**m1.variables, **m2.variables}, key=lambda var: var.symbol, reverse=True):
            if sym in m1.variables and sym not in m2.variables:
                return -1
            elif sym not in m1.variables and sym in m2.variables:
                return 1
            elif sym in m1.variables and sym in m2.variables:
                return m2.variables[sym] - m1.variables[sym]

    return 0