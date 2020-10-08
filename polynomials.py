from functools import cmp_to_key
from fractions import Fraction

from orderings import grevlex

COEFFICIENT_TYPES = (int, Fraction) # For now.
ORDERING = grevlex

class Variable:
    """ This class is unnecessary. I could have just used chars. """

    # Maps symbols to objects. If someone tries to construct a Variable with a used symbol, we'll
    # give them the previously-created one instead.
    other_vars = {}
    def __new__(cls, symbol):
        if symbol in Variable.other_vars:
            return Variable.other_vars[symbol]
        return super().__new__(cls)

    def __init__(self, symbol):
        assert symbol in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

        self.symbol = symbol
        Variable.other_vars[symbol] = self

    def __eq__(self, other):
        return isinstance(other, symbol) and self.symbol == other.symbol

    def __hash__(self):
        return hash(self.symbol)

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return self.symbol


class Monomial:
    def __init__(self, coefficient, variables={}): # Mutable default parameters are usually a big mistake in Python, but it's actually ok here.
        assert any(isinstance(coefficient, coefficient_type) for coefficient_type in COEFFICIENT_TYPES)
        assert all(isinstance(var, Variable) and isinstance(variables[var], int) for var in variables)

        self.coefficient = coefficient
        self.variables = variables if self.coefficient != 0 else {}
        # I would copy variables to avoid accidentally changing this object, but Python attributes are all public anyway

        # Clean up the variables dict in case they gave us 0 power variables
        to_delete = []
        for var in self.variables:
            if self.variables[var] == 0:
                to_delete.append(var)
        for var in to_delete:
            del self.variables[var]

        self.degree = sum(self.variables[var] for var in self.variables)

    def __add__(self, other):
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            other = Monomial(other)
        
        if isinstance(other, Monomial):
            if self.variables == other.variables:
                return Monomial(self.coefficient + other.coefficient, self.variables)
            else:
                return Polynomial(self, other)
        elif isinstance(other, Polynomial):
            return Polynomial(self) + other

    def __sub__(self, other):
        return self + other * -1

    def __mul__(self, other):
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            return Monomial(self.coefficient * other, self.variables)
        
        elif isinstance(other, Monomial):
            coefficient = self.coefficient * other.coefficient

            new_variables = {}
            for var in {**self.variables, **other.variables}:
                new_variables[var] = 0
                if var in self.variables:
                    new_variables[var] += self.variables[var]
                if var in other.variables:
                    new_variables[var] += other.variables[var]

            return Monomial(coefficient, new_variables)

        elif isinstance(other, Polynomial):
            return Polynomial(self) * other

    def __truediv__(self, other):
        coefficient = Fraction(self.coefficient, other.coefficient)

        new_variables = {}

        for var in {**other.variables, **self.variables}:
            # If we're dividing by something with a variable we don't have, then we can't do the division
            if var not in self.variables:
                return 0, self

            new_power = self.variables[var] - (other.variables[var] if var in other.variables else 0)

            if new_power < 0:
                return 0, self

            if new_power != 0:
                new_variables[var] = new_power

        return Monomial(coefficient, new_variables), 0

    def __pow__(self, exp):
        assert isinstance(exp, int)

        new_variables = self.variables.copy()
        for var in new_variables:
            new_variables[var] *= exp
        return Monomial(coefficient * exp, new_variables)

    def __eq__(self, other):
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            return self.coefficient == other and all(self.variables[var] == 0 for var in self.variables)
        elif isinstance(other, Monomial):
            return self.coefficient == other.coefficient and self.variables == other.variables
        elif isinstance(other, Polynomial):
            return Polynomial(self) == other
        else:
            raise TypeError("Cannot compare Monomial and", type(other))

    def __str__(self):
        return str(self.coefficient) + \
               "".join(str(var) + ('^' + str(self.variables[var]) if self.variables[var] != 1 else "") \
                    for var in sorted(self.variables, key=lambda var: var.symbol))

    def __repr__(self):
        return str(self.coefficient) + \
               "".join(str(var) + ('^' + str(self.variables[var]) if self.variables[var] != 1 else "") \
                    for var in sorted(self.variables, key=lambda var: var.symbol))

    def __hash__(self):
        return hash(repr(self))


def monomial_lcm(m1, m2):
    """ Note: The zero-checks here aren't necessary if the Monomial object is being used correctly. """
    new_variables = {}
    for var in {**m1.variables, **m2.variables}:
        new_variables[var] = 0
        if var in m1.variables and m1.variables[var] != 0:
            new_variables[var] = m1.variables[var]
        if var in m2.variables:
            new_variables[var] = max(new_variables[var], m2.variables[var])
        if new_variables[var] == 0:
            del new_variables[var]

    return Monomial(1, new_variables)


class Polynomial:
    def __init__(self, *monomials):
        assert monomials != () and all(isinstance(monomial, Monomial) for monomial in monomials)

        self.monomials = []

        monomials = list(monomials)
        monomials.sort(key=cmp_to_key(ORDERING), reverse=True)

        self.monomials.append(monomials[0])
        for monomial in monomials[1:]:
            if ORDERING(self.monomials[-1], monomial) == 0:
                self.monomials[-1] = monomial + self.monomials[-1]
            else:
                self.monomials.append(monomial)

        i = 0
        while i < len(self.monomials):
            if self.monomials[i].coefficient == 0:
                del self.monomials[i]
            else:
                i += 1

        if self.monomials == []:
            self.monomials = [Monomial(0)]

    def __mul__(self, other):
        """ Simple n^2 multiplication. Better algorithms exist, but this is easier. """
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            other = Polynomial(Monomial(other))
        elif isinstance(other, Monomial):
            other = Polynomial(other)
        
        if isinstance(other, Polynomial):
            products = []
            for m1 in self.monomials:
                for m2 in other.monomials:
                    products.append(m1 * m2)
            return Polynomial(*products)
        else:
            raise TypeError("Cannot multiply Polynomial with", type(other))

    def __truediv__(self, other):
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            other = Polynomial(Monomial(other))
        elif isinstance(other, Monomial):
            other = Polynomial(other)

        quotients = []
        remainders = []

        curr = self
        while curr != 0:
            curr_monomial = leading_term(curr)
            quot, rem = curr_monomial / leading_term(other)
            if rem != 0:
                remainders.append(curr_monomial)
                curr -= curr_monomial
            else:
                quotients.append(quot)
                curr -= other * quot

        quot = Polynomial(*quotients) if quotients != [] else Monomial(0)
        rem = Polynomial(*remainders) if remainders != [] else Monomial(0)
        return quot, rem

    def __add__(self, other):
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            other = Polynomial(Monomial(other))
        elif isinstance(other, Monomial):
            other = Polynomial(other)
        
        if isinstance(other, Polynomial):
            new_monomials = []
            i = j = 0
            while i < len(self.monomials) and j < len(other.monomials):
                if ORDERING(self.monomials[i], other.monomials[j]) < 0:
                    new_monomials.append(other.monomials[j])
                    j += 1
                elif ORDERING(self.monomials[i], other.monomials[j]) > 0:
                    new_monomials.append(self.monomials[i])
                    i += 1
                else:
                    new_monomials.append(self.monomials[i] + other.monomials[j])
                    i += 1
                    j += 1

            while i < len(self.monomials):
                new_monomials.append(self.monomials[i])
                i += 1
            while j < len(other.monomials):
                new_monomials.append(other.monomials[j])
                j += 1

            return Polynomial(*new_monomials) # I think this will already be in order, so maybe we can bypass the reordering in __init__?
        else:
            raise TypeError("Cannot add Polynomial with", type(other))

    def __eq__(self, other):
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            return len(self.monomials) == 1 and self.monomials[0] == other
        elif isinstance(other, Monomial):
            return self == Polynomial(other)
        elif isinstance(other, Polynomial):
            return self.monomials == other.monomials
        else:
            raise TypeError("Cannot compare Polynomial with", type(other))

    def __sub__(self, other):
        return self + other * -1

    def __str__(self):
        if self.monomials == []:
            return "1"
        return " + ".join(map(str, self.monomials))

    def __repr__(self):
        if self.monomials == []:
            return "1"
        return " + ".join(map(str, self.monomials))

    def __hash__(self):
        return hash(str(self))


def leading_term(p):
    return p.monomials[0]

def leading_monomial(p):
    return Monomial(1, p.monomials[0].variables)

def leading_coefficient(p):
    return leading_term(p).coefficient

def s_polynomial(p1, p2):
    lcm = monomial_lcm(leading_monomial(p1), leading_monomial(p2))
    q1, r1 = lcm / leading_term(p1)
    q2, r2 = lcm / leading_term(p2)
    
    assert r1 == r2 == 0

    return q1 * p1 - q2 * p2

def lead_reducible(p1, p2):
    _, r = leading_monomial(p1) / leading_monomial(p2)
    return r == 0

def divisible(p1, p2):
    _, r = p1 / p2
    return r == 0

def buchberger(*polynomials):
    """ Takes some number of polynomials, returns a Groebner basis for them. """
    polynomials = list(polynomials)
    tried = set()
    while True:
        to_add = None
        for i in range(len(polynomials)):
            for j in range(i + 1, len(polynomials)):
                p1 = polynomials[i]
                p2 = polynomials[j]
                if (p1, p2) in tried:
                    continue
                tried.add((p1, p2))

                s = s_orig = s_polynomial(p1, p2)

                # Divide s by our polynomials until r == 0 or we're out of polynomials.
                # If we run out of polynomials, that means our GB is incomplete!
                for p in polynomials:
                    _, r = s / p
                    if r == 0:
                        break
                    else:
                        s = r
                else:
                    to_add = s_orig
                    break # Double break. Makes me long for goto
            else:
                continue
            break # continuation of the double break

        if to_add is not None:
            polynomials.append(to_add)
        else:
            break

    return polynomials

def minimize_gb(*polynomials):
    """ This probably doesn't work """
    polynomials = list(polynomials)

    tried = set()
    while True:
        index_to_delete = None
        for i in range(len(polynomials)):
            for j in range(len(polynomials)):
                if i == j:
                    continue
                p1 = polynomials[i]
                p2 = polynomials[j]

                if (p1, p2) in tried:
                    continue
                tried.add((p1, p2))

                if lead_reducible(p1, p2):
                    index_to_delete = i
                    break
            if index_to_delete is not None:
                del polynomials[index_to_delete]
                break
        if index_to_delete is None:
            break

    # Then make all leading monomials monic
    for i in range(len(polynomials)):
        polynomials[i] = (polynomials[i] / leading_coefficient(polynomials[i]))[0]

    return polynomials

class Ideal:
    """ Definitely incomplete """
    def __init__(self, *generators):
        assert generators != () and all(isinstance(g, Polynomial) for g in generators)
        self.generators = generators

    def __str__(self):
        return f"ideal{self.generators}"

    def __contains__(self, p):
        assert isinstance(p, Polynomial)

        return any(lead_reducible(p, g) for g in buchberger(*self.generators))


# For testing
x = Variable("x")
y = Variable("y")
z = Variable("z")

m1 = Monomial(2, {x: 1, z: 1})
m2 = Monomial(3, {y: 2})
m3 = Monomial(1, {x: 1})
m4 = Monomial(1, {x: 3, y: 2, z: 2})

p1 = Polynomial(m1, m2)
p2 = Polynomial(m3)
p3 = Polynomial(m4, m3)