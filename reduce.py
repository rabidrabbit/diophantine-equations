import math
from constants import Constants, padic_log, padic_order
from itertools import combinations
from sage.all import *

class BoundReduce:
    def __init__(self, constants, tries=50):
        self.constants = constants
        self.coefficients = constants.calculate_constants()
        self.tries = tries

    def real_reduce(self, bound, large_constant):
        """
        Integral version of the LLL algorithm, where no rational operations are done.
        See de Weger's thesis for more information.
        """
        approximation_matrix = self.generate_approximation_matrix(large_constant)
        LLL_matrix = approximation_matrix.transpose().LLL().transpose()
        print(LLL_matrix)
        return

    def generate_approximation_matrix(self, large_constant):
        n = len(self.constants.primes)
        primes_row = [round(large_constant * log(p)) for p in self.constants.primes]
        approximation_matrix = []
        for i in range(n):
            zero_row = [0] * n
            zero_row[i] = 1
            approximation_matrix.append(zero_row)
        approximation_matrix[n - 1] = primes_row
        return Matrix(ZZ, approximation_matrix)

    def padic_reduce(self, bound):
        """
        Employs methodology from Pink and Zieglier (2016).
        """
        # Extremely unlikely edge case.
        if bound < self.constants.num_terms:
            raise ValueError("bounds are already sufficient.")

        Z_bounds = []
        for i in range(len(self.constants.primes)):
            p = self.constants.primes[i]
            r = math.ceil(math.log(self.coefficients["n1_bound"], p))
            prec = r + self.tries + 10
            L = Qp(p, prec)

            current_z_bound = -1
            for t_vec in combinations(list(range(1, bound + 1)), self.constants.num_terms):
                alpha, beta = self.calculate_alphabeta(p, prec)

                z0_left_term = self.constants.a * (1 + sum([alpha ** t for t in t_vec]))
                z0_right_term = alpha - beta
                z0 = L(z0_left_term.norm()).ordp() - L(z0_right_term.norm()).ordp()

                tau = self.calculate_tau(t_vec, alpha, beta)

                # Unlikely case, and we may assume the p-adic log is injective for our purposes.
                if tau == 1:
                    p_order = log(alpha / beta).norm().ordp() / 2
                    new_bound = math.log(self.coefficients["n1_bound"], p) + p_order + z0
                    Z_bounds.append(max(z0 + 3/2, new_bound))
                else:
                    # Remember that zeta is in Q_p, not in the field extension.
                    zeta = log(tau) / log(alpha / beta)
                    vzeta = zeta.norm().ordp()
                    zeta_list = list(zeta.expansion())

                    R_max = self.find_Rmax(50, vzeta, zeta_list)
                    if R_max == -1:
                        m0 = self.find_m0(p, zeta_list)
                        z_bound = math.log(self.coefficients["n1_bound"] - m0, p) + (alpha / beta).ordp() + z0
                    else:
                        zeta_inverse = log(alpha / beta) / log(tau)
                        vzeta_inverse = zeta_inverse.ordp()
                        z_bound = log(tau).ordp() + R_max + vzeta_inverse
                    current_z_bound = max(current_z_bound, z_bound)
            print(current_z_bound)
            Z_bounds.append(current_z_bound)
        return Z_bounds

    def cont_fraction_reduce(self):
        return

    def calculate_tau(self, t_vec, alpha, beta):
        numerator = 1 + sum([alpha ** t for t in t_vec])
        denominator = 1 + sum([beta ** t for t in t_vec])
        return numerator / denominator

    def calculate_alphabeta(self, p, prec):
        """
        Returns sqrt(delta) in p-adic representation
        """
        var('x')
        sqrtdelta = None
        try:
            M = Qp(p, prec).extension(x ** 2 - self.constants.delta, names="padicroot")
            sqrtdelta = K.gen(0)
        except NotImplementedError:
            try:
                M = Qp(p, prec)
                sqrtdelta = M(self.constants.delta).sqrt()
            except:
                # Exceptional case (i.e. p = 2).
                M = Qp(p, prec).extension(x ** 2 - self.constants.A * x - self.constants.B, names="padicroot")
                alpha = M.gen(0)
                beta = self.constants.A - alpha
                return (alpha, beta)

        alpha = (self.constants.A + sqrtdelta) / 2
        beta = self.constants.A - alpha
        return (alpha, beta)

    def find_Rmax(self, r, vzeta, zeta_list, tries=50):
        Rmax = -1
        for i in range(r, r + tries):
            # Checks if i is a valid index
            if zeta_list[i - vzeta] != []:
                if i > Rmax:
                    Rmax = i
            elif i == r + tries - 1:
                print("warning: no R found...")
        return Rmax

    def find_m0(self, p, zeta_list):
        # In this case, it must be the case that zeta is an integer.
        m0 = 0
        for i in range(len(zeta_list)):
            if zeta_list[i] != []:
                m0 += (p ** i) * zeta_list[i][0]
        return m0

if __name__ == "__main__":
    constants_gen = Constants(
        a = 1,
        b = 1,
        A = 1,
        B = 1,
        alpha = (1 + math.sqrt(5))/2,
        beta = (1 - math.sqrt(5))/2,
        delta = 5,
        num_terms = 3,
        w = 1,
        primes = [2, 3, 5]
    )

    br = BoundReduce(constants_gen)
    br.real_reduce(500, 10**100)
    #br.padic_reduce(3000)
