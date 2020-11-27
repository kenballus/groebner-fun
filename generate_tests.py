from polynomials import *
from random import randint

of = open("test_cases.txt", "w")

for _ in range(10):
    gens = [random_polynomial() for _ in range(randint(2, 5))]
    query = random_polynomial()
    print(gens)
    print(query)
    print(query in Ideal(*gens))