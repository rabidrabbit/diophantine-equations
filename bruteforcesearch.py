"""
Implementation of a brute-force search for our Diophantine equation,
given that the bounds on the coefficients are sufficiently small.

@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""


class BruteForceSearch:
    def __init__(self, a, b, A, B, alpha, beta, num_terms, primes, debug_flag=False):
        self.a = a
        self.b = b
        self.A = A
        self.B = B
        self.alpha = alpha
        self.beta = beta
        self.num_terms = num_terms
        self.primes = primes

    def search(self, bound):
        """
        Parallel implementation of brute-force search.
        """
        memo = []
        solutions = []

        return solutions