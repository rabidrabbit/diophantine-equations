"""
Implementation of a brute-force search for our Diophantine equation,
given that the bounds on the coefficients are sufficiently small.

@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""
import math
import multiprocessing
import itertools
from itertools import combinations, islice


class BruteForceSearch:
    def __init__(self, u0, u1, a, b, A, B, num_terms, primes, chunk_size=256, debug_flag=False):
        self.u0 = u0
        self.u1 = u1
        self.a = a
        self.b = b
        self.A = A
        self.B = B
        self.num_terms = num_terms
        self.primes = primes

        self.chunk_size = 256

    def search(self, bound):
        memo = self.pre_compute_memoization(bound)
        solutions = []

        for t_vec in combinations(list(range(1, bound + 1)), self.num_terms + len(self.primes)):
            sum_term = sum([memo[t] for t in t_vec[:self.num_terms]])
            product_term = math.prod(self.primes[i] ** t_vec[i + self.num_terms - 1] for i in range(len(self.primes)))
            if sum_term == product_term:
                solutions.append(t_vec)
        return solutions

    def pre_compute_memoization(self, max_param):
        memo = [self.u0, self.u1]
        a = self.u0
        b = self.u1
        next_term = self.u0
        for i in range(max_param - 2):
            next_term = self.B * a + self.A * b
            a = b
            b = next_term
            memo.append(next_term)
        return memo

if __name__ == "__main__":
    brute_force_search = BruteForceSearch(
            u0 = 0,
            u1 = 1,
            a = 1,
            b = 1,
            A = 1,
            B = 1,
            num_terms = 2,
            primes = [2, 3, 5]
    )

    solutions = brute_force_search.search(100)
    print(solutions)
