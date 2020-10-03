from string import ascii_uppercase, ascii_lowercase
from collections import defaultdict
from copy import copy
from functools import cmp_to_key
from fractions import Fraction

from orderings import grevlex

COEFFICIENT_TYPES = (int, Fraction) # For now.

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
        assert symbol in ascii_lowercase + ascii_uppercase

        self.symbol = symbol
        Variable.other_vars[symbol] = self

    def __eq__(self, other):
        return self.symbol == other.symbol

    def __hash__(self):
        return hash(self.symbol)

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return self.symbol

class Monomial:
    """ As it is, this class doesn't do well if you give it a zero-power variable. """
    def __init__(self, coefficient, variables={}): # Mutable default parameters are usually a big mistake in Python, but it's actually ok here.
        assert any(isinstance(coefficient, coefficient_type) for coefficient_type in COEFFICIENT_TYPES)
        assert all(isinstance(var, Variable) and isinstance(variables[var], int) for var in variables)

        self.coefficient = coefficient
        self.variables = copy(variables)

        self.degree = sum(self.variables[var] for var in self.variables)

    def __add__(self, other):
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            return self + Monomial(other)
        elif isinstance(other, Monomial):
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

            new_variables = defaultdict(lambda: 0)
            for var in {**self.variables, **other.variables}:
                if var in self.variables:
                    new_variables[var] += self.variables[var]
                if var in other.variables:
                    new_variables[var] += other.variables[var]

            return Monomial(coefficient, dict(new_variables))

        elif isinstance(other, Polynomial):
            return Polynomial(self) * other

    def divmod(self, other):
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

        new_variables = copy(self.variables)
        for var in new_variables:
            new_variables[var] *= exp
        return Monomial(coefficient * exp, new_variables)

    def __eq__(self, other):
        return self.coefficient == other.coefficient and self.variables == other.variables

    def __str__(self):
        return (str(self.coefficient) if self.coefficient != 1 else "") + \
               "".join(str(var) + ('^' + str(self.variables[var]) if self.variables[var] != 1 else "") \
                    for var in sorted(self.variables, key=lambda var: var.symbol))

    def __repr__(self):
        return (str(self.coefficient) if self.coefficient != 1 else "") + \
               "".join(str(var) + ('^' + str(self.variables[var]) if self.variables[var] != 1 else "") \
                    for var in sorted(self.variables, key=lambda var: var.symbol))

    def __hash__(self):
        return hash(repr(self))

class Polynomial:
    def __init__(self, *monomials):
        assert monomials != [] and all(isinstance(monomial, Monomial) for monomial in monomials)
        
        self.ordering = grevlex
        self.monomials = []

        self.leading_monomial = self.monomials[0]

        monomials = list(monomials)
        monomials.sort(key=cmp_to_key(self.ordering), reverse=True)

        self.monomials.append(monomials[0])
        for monomial in monomials[1:]:
            if self.ordering(self.monomials[-1], monomial) == 0:
                self.monomials[-1] = monomial + self.monomials[-1]
            else:
                self.monomials.append(monomial)

    def __mul__(self, other):
        """ Simple n^2 multiplication. Better algorithms exist, but this is easier. """
        if any(isinstance(other, coefficient_type) for coefficient_type in COEFFICIENT_TYPES):
            return self * Monomial(other)
        elif isinstance(other, Monomial):
            return self * Polynomial(other)
        elif isinstance(other, Polynomial):
            new_monomials = []
            for i in range(len(self.monomials)):
                for j in range(len(other.monomials)):
                    new_monomials.append(self.monomials[i] * other.monomials[j])

            return Polynomial(*new_monomials)

    def divmod(self, other):
        new_monomials = []
        remainders = []

        for monomial in self.monomials:
            quot, rem = divmod(other.leading_monomial, monomial)
            

    def __add__(self, other):
        new_monomials = []
        i = j = 0
        while i < len(self.monomials) and j < len(other.monomials):
            if self.ordering(self.monomials[i], other.monomials[j]) < 0:
                new_monomials.append(other.monomials[j])
                j += 1
            elif self.ordering(self.monomials[i], other.monomials[j]) > 0:
                new_monomials.append(self.monomials[i])
                i += 1
            else:
                new_monomials.append(self.monomials[i] + other.monomials[j])
                i += 1
                j += 1

        return Polynomial(*new_monomials) # I think this will already be in order, so maybe we can bypass the reordering in __init__?

    def __sub__(self, other):
        return self + other * -1

    def __str__(self):
        return " + ".join(map(str, self.monomials))

    def __repr__(self):
        return " + ".join(map(str, self.monomials))



# For testing
x = Variable("x")
y = Variable("y")
z = Variable("z")

m1 = Monomial(2, {x: 1, z: 1})
m2 = Monomial(3, {y: 2})
m3 = Monomial(1, {x: 1})
m4 = Monomial(9, {x: 10, y: 11, z: 100})

p1 = Polynomial(m1, m2)
p2 = Polynomial(m3)