"""
Computes the constants that are needed to 
find an upper-bound on our parameters.

@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""
import math
from fractions import Fraction
from sage.all import *
from sage.rings.number_field.S_unit_solver import log_p

def padic_log(value, ideal, prec=50):
    log_val = log_p(value, ideal, prec)
    return log_val

def padic_order(value, p):
    L = Qp(p)
    return L(value).ordp()

class Constants:
    def __init__(self, a, b, A, B, alpha, beta, num_terms, w, primes, debug_flag=False):
        self.a = a
        self.b = b
        self.A = A
        self.B = B
        self.alpha = alpha
        self.beta = beta
        self.num_terms = num_terms
        self.delta = abs(alpha - beta)
        self.sqrtdelta = math.sqrt(self.delta)
        self.w = w
        self.primes = primes

        self.constants = {}
    
    def calculate_constants(self):
        c1 = self.num_terms * (abs(self.a) + abs(self.b))/self.delta
        self.constants["c1"] = c1

        c2 = self.log_star(c1) / (math.log(self.alpha) * abs(self.w))
        self.constants["c2"] = c2

        golden_ratio = (1 + math.sqrt(5))/2
        c3_left_term = 2
        c3_right_term = 2
        c3 = max(c3_left_term, c3_right_term)
        self.constants["c3"] = c3

        c4 = abs(self.a) * (1 - 1 / golden_ratio) / (2 * self.sqrtdelta)
        self.constants["c4"] = c4

        c5 = (math.log(abs(self.w)) - math.log(c4)) / math.log(self.alpha)
        self.constants["c5"] = c5

        c6_list = self.calculate_c6_list()
        self.constants["c6_list"] = c6_list

        c8_list = self.calculate_c8_list()
        self.constants["c8_list"] = c8_list

        c9_list = self.calculate_c9_list(c6_list)
        self.constants["c9_list"] = c9_list

        c10_list = self.calculate_c10_list(c5, c8_list)
        self.constants["c10_list"] = c10_list

        d0 = self.calculate_d0()
        self.constants["d0"] = d0

        d1 = self.calculate_d1(d0)
        self.constants["d1"] = d1

        d2 = self.calculate_d2()
        self.constants["d2"] = d2

        N_list = self.calculate_Nt()

        dt_list = self.calculate_Dt(C5_list, N_list, d1)

        return self.constants

    def calculate_c6_list(self):
        """
        Calculates the constants c_{6, i} for 1 <= i <= s
        """
        c6_list = []
        for p in self.primes:
            term = math.log(2 * abs(self.a) * self.alpha) / math.log(p)
            c6_list.append(term)
        return c6_list
    
    def calculate_c8_list(self):
        c8_list = []
        c7 = max(math.log(abs(2 * self.a * self.alpha)), math.log(abs(2 * self.b * self.beta)))
        for p in self.primes:
            c8_list.append(max(c7, math.log(p)))
        return c8_list
    
    def calculate_c9_list(self, c6_list):
        var('x')
        K = NumberField(x ** 2 - self.A * x - self.B, name='root1')

        c9_list = []
        for i in range(len(c6_list)):
            """
            Uses the approach from X et al., as opposed to the original
            implementation of the p-adic logarithm in SAGE
            L(log_p(a/(1 - a), p1, 50).absolute_norm()).ordp()
            """
            p = self.primes[i]
            p_ideal = K.primes_above(p)[0]
            alpha_over_beta = K.gen(0) / (K.gen(0) - self.A)
            print(alpha_over_beta)
            p_order = padic_order(padic_log(alpha_over_beta, p_ideal, prec=50).norm(), p)
            term = p_order + (2 / math.log(p)) + c6_list[i]
            c9_list.append(term)
        return c9_list

    def calculate_c10_list(self, c5, c8_list):
        c10_list = []
        for i in range(len(c8_list)):
            sum_l = sum([c8_list[i] * math.log(p) for p in self.primes])
            term = sum_l / math.log(self.alpha) + c5
            c10_list.append(term)
        return c10_list

    def calculate_C1t(self, d0):
        return d0 * (num_terms - t + 1) + (t - 1) * abs(self.b) / self.sqrtdelta
    
    def calculate_C2t(self, t, d0):
        return (t - 1) * abs(self.b / self.a) + (num_terms - t + 1) * d0 * self.sqrtdelta / abs(a)
    
    def calculate_C3t(self, t, d0):
        return calculate_C1t(t, d0) * self.sqrtdelta / abs(a)

    def calculate_C4t(self, t, d0):
        return max(calculate_C3t(t, d0), calculate_C2t(t, d0))
    
    def calculate_C5t(self, t, d2):
        return 2 * (math.log(d2) + 0.5 * math.log(self.delta)) + t * math.log(4)
        
    def calculate_Nt(self, ):
        N_list = [None]
        return

    def calculate_Dt(self, C5_list, N_list, d1):
        Dt_list = [None, None]

        for t in range(2, self.num_terms):
            coefficient = 1.4 * (30 ** 5) * (2 ** 4.5) * 4
            numeric_constant = coefficient * 4 * (1 + math.log(2)) * (1 + math.log(d1)) * (2 * math.log(self.alpha))
            prime_constant = math.prod([2 * math.log(p) for p in self.primes])
            t_term = C5_list[t - 2]
            term = numeric_constant * prime_constant * t_term
            Dt_list.append(term)

        return Dt_list

    def calculate_d0(self):
        return (abs(self.a) + abs(self.b)) / self.sqrtdelta

    def calculate_d1(self, d0):
        return 1 + max([d0 * self.alpha / math.log(p) for p in self.primes])
    
    def calculate_d2(self):
        max_a = max(abs(a), 1/abs(a))
        max_b = max(abs(b), 1/abs(b))
        return max(max_a, max_b)

    def log_star(self, x):
        return math.log(max(x, 1))


if __name__ == "__main__":
    constants = Constants(
        a = 1,
        b = 1,
        A = 1,
        B = 1,
        alpha = (1 + sqrt(5))/2,
        beta = (1 - sqrt(5))/2,
        num_terms = 3,
        w = 1,
        primes = [2, 3]
    )

    c = constants.calculate_constants()
    print(c)