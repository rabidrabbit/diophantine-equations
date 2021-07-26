"""
Computes the constants that are needed to 
find an upper-bound on our parameters.

@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""
import math
from fractions import Fraction
from sage.all import *
from sage.rings.number_field.S_unit_solver import log_p

def get_primes_list(min_num, max_num):
    return list(primes(max_num))

def padic_log(value, ideal, prec=50):
    log_val = log_p(value, ideal, prec)
    return log_val

def padic_order(value, p):
    L = Qp(p)
    return L(value).ordp()

class Constants:
    def __init__(self, a, b, A, B, alpha, beta, delta, num_terms, w, primes, debug_flag=False):
        self.a = a
        self.b = b
        self.A = A
        self.B = B
        self.alpha = alpha
        self.beta = beta
        self.num_terms = num_terms
        self.delta = delta
        self.sqrtdelta = abs(alpha - beta)
        self.w = w
        self.primes = primes

        self.constants = {}

    def get_coefficients(self):
        return self.constants
    
    def calculate_constants(self):
        c1 = self.num_terms * (abs(self.a) + abs(self.b))/self.sqrtdelta
        self.constants["c1"] = c1

        c2 = self.log_star(c1 / abs(self.w)) / math.log(self.alpha)
        self.constants["c2"] = c2

        golden_ratio = (1 + math.sqrt(5))/2
        c3_left_term = self.log_star(2 * self.num_terms * abs(self.b) * golden_ratio / (abs(self.a) * (golden_ratio - 1))) / math.log(self.alpha)
        c3_right_term = self.log_star(2 * self.num_terms * abs(self.b) * golden_ratio / (abs(self.a) * (golden_ratio - 1))) / math.log(abs(self.alpha / self.beta))
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

        c11_list = self.calculate_c11_list(c6_list, c8_list, c9_list)
        self.constants["c11_list"] = c11_list

        c10_list = self.calculate_c10_list(c5, c11_list)
        self.constants["c10_list"] = c10_list

        d0 = self.calculate_d0()
        self.constants["d0"] = d0

        d1 = self.calculate_d1(d0)
        self.constants["d1"] = d1

        d2 = self.calculate_d2()
        self.constants["d2"] = d2

        C5t_list = self.calculate_C5t_list(d2)

        Nt_list = self.calculate_Nt(2, C5t_list, d0, d1)
        self.constants["Nt_list"] = Nt_list

        C6 = self.calculate_C6(c3)
        self.constants["C6"] = C6

        n1_bound = self.calculate_n1_bound(C6, c10_list, Nt_list)
        self.constants["n1_bound"] = n1_bound

        Z_bounds = self.calculate_Z_bounds(n1_bound)
        self.constants["Z_bounds"] = Z_bounds

        return self.constants

    def update_constants(self, diff_bound):
        """
        Given a new bound on n_1 - n_k, recalculate the bounds on n_1 and z_i.
        """
        self.constants["Nt_list"][-1] = diff_bound
        n1_bound = self.calculate_n1_bound(
                self.constants["C6"],
                self.constants["c10_list"],
                self.constants["Nt_list"]
        )
        self.constants["n1_bound"] = n1_bound

        Z_bounds = self.calculate_Z_bounds(n1_bound)
        self.constants["Z_bounds"] = Z_bounds
        return

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
            """
            p = self.primes[i]
            p_ideal = K.primes_above(p)[0]
            alpha_over_beta = K.gen(0) / (self.A - K.gen(0))
            p_order = padic_order(padic_log(alpha_over_beta, p_ideal, prec=50).norm(), p) / 2
            term = p_order + (2 / math.log(p)) + c6_list[i]
            c9_list.append(term)
        return c9_list

    def calculate_c10_list(self, c5, c11_list):
        c10_list = []
        for i in range(len(self.primes)):
            sum_l = sum([c11_list[i] * math.log(p) for p in self.primes])
            term = sum_l / math.log(self.alpha) + c5
            c10_list.append(term)
        return c10_list
    
    def calculate_c11_list(self, c6_list, c8_list, c9_list):
        c11_list = []
        for i in range(len(self.primes)):
            p = self.primes[i]
            max_term = max(math.log(self.alpha), math.log(p))
            left_term = self.calculate_bugeaud_constant(p) * max_term * c8_list[i] + c6_list[i]
            term = max(left_term, c9_list[i])
            c11_list.append(term)

        return c11_list 

    def calculate_C6(self, c3):
        P = max(self.primes)
        sum_term = max(math.log(abs(2 * self.a * self.alpha)), math.log(abs(2 * self.b * self.beta))) + 0.24 / (self.num_terms - 1)
        term = 17.5 * (self.num_terms - 1) * math.log(self.alpha) * sum_term
        C6 = max(c3, term, P ** 10)
        return C6

    def calculate_C1t(self, t, d0):
        if t < 2:
            raise ValueError("Invalid value of t.")
        return d0 * (self.num_terms - t + 1) + (t - 1) * abs(self.b) / self.sqrtdelta
    
    def calculate_C2t(self, t, d0):
        if t < 2:
            raise ValueError("Invalid value of t.")
        return (t - 1) * abs(self.b / self.a) + (self.num_terms - t + 1) * d0 * self.sqrtdelta / abs(self.a)
    
    def calculate_C3t(self, t, d0):
        if t < 2:
            raise ValueError("Invalid value of t.")
        return self.calculate_C1t(t, d0) * self.sqrtdelta / abs(self.a)

    def calculate_C4t(self, t, d0):
        if t < 2:
            raise ValueError("Invalid value of t.")
        return max(self.calculate_C3t(t, d0), self.calculate_C2t(t, d0))
    
    def calculate_C5t_list(self, d2):
        C5t_list = [None]
        for t in range(2, self.num_terms + 1):
            term = 2 * (math.log(d2) + math.log(self.sqrtdelta)) + t * math.log(4)
            C5t_list.append(term)
        return C5t_list
        
    def calculate_Nt(self, t, C5t_list, d0, d1):
        Nt_list = [None]

        for t in range(2, self.num_terms + 1):
            numerator = self.calculate_Dt(t, Nt_list, C5t_list, d0, d1) + math.log(self.calculate_C4t(t, d0))
            alpha_over_beta = self.alpha / abs(self.beta)
            denominator = math.log(min(alpha_over_beta, self.alpha))
            term = numerator / denominator
            Nt_list.append(term)

        return Nt_list

    def calculate_Dt(self, t, Nt_list, C5t_list, d0, d1):
        num_primes = len(self.primes)
        coefficient = 1.4 * (30 ** (num_primes + 5)) * ((num_primes + 2) ** 4.5) * 4
        numeric_constant = coefficient * (1 + math.log(2)) * (1 + math.log(d1)) * (2 * math.log(self.alpha))
        prime_constant = math.prod([2 * math.log(p) for p in self.primes])
        t_term = C5t_list[t - 1]
        sum_Nt = 0 if t == 2 else sum(Nt_list[1:(t - 1)])

        term = numeric_constant * prime_constant * (t_term + 2 * math.log(self.alpha) * sum_Nt)

        return term

    def calculate_d0(self):
        return (abs(self.a) + abs(self.b)) / self.sqrtdelta

    def calculate_d1(self, d0):
        return 1 + max([d0 * self.alpha / math.log(p) for p in self.primes])
    
    def calculate_d2(self):
        max_a = max(abs(self.a), 1/abs(self.a))
        max_b = max(abs(self.b), 1/abs(self.b))
        return max(max_a, max_b) + 1

    def calculate_Z_bounds(self, n1_bound):
        """
        Calculates bounds on z_i (i.e. the parameter on the exponents of the primes)
        """
        Z_bounds = []
        for p in self.primes:
            coefficient = 2 * math.log(self.alpha) / math.log(p)
            term = coefficient * n1_bound
            Z_bounds.append(term)
        return Z_bounds

    def calculate_n1_bound(self, C6, c10_list, Nt_list):
        """
        Calculates bounds on n_1 (i.e. the parameters of the recurrence sequence).
        """
        n1_bound = -1
        power_of_2 = 2 ** (self.num_terms + 1)
        Nk = Nt_list[self.num_terms - 1]
        for i in range(len(self.primes)):
            kth_root_term = math.pow(c10_list[i] * (self.num_terms - 1) * Nk, 1/(self.num_terms + 1))
            log_term = math.log( (self.num_terms + 1) ** (self.num_terms + 1) * c10_list[i] * (self.num_terms - 1) * Nk )
            first_term = power_of_2 * (kth_root_term * log_term) ** (self.num_terms + 1)
            euler_term = (2 * math.e) ** (2 * self.num_terms + 2)
            n1_bound = max(n1_bound, first_term, euler_term, C6)
        return n1_bound

    def calculate_bugeaud_constant(self, prime):
        """
        Calculates Bugeaud constant as defined in the paper and in their paper
        on bounding a linear form in two p-adic logarithms.
        """
        # First, indirectly calculates the residual degree of the extension by
        # calculating the ramification index
        square_free_part = squarefree_part(self.delta)
        residual_degree = 1 if square_free_part % prime == 0 else 2

        term = 947 * (prime ** residual_degree) / ((math.log(prime) ** 4))
        return term

    def log_star(self, x):
        return math.log(max(x, 1))


if __name__ == "__main__":
    constants = Constants(
        a = 1,
        b = 1,
        A = 1,
        B = 1,
        delta = 5,
        alpha = (1 + sqrt(5))/2,
        beta = (1 - sqrt(5))/2,
        num_terms = 2,
        w = 1,
        primes = get_primes_list(1, 200)
    )

    c = constants.calculate_constants()
    print(c)
