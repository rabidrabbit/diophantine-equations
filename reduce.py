import math
from constants import Constants, log_p
from itertools import combinations
from sage.all import *

class BoundReduce:
    def __init__(self, constants):
        print("initialized.")
        self.constants = constants
        self.coefficients = constants.calculate_constants()
        print("done calculating")

    def real_reduce(self):
        return

    def padic_reduce(self, bound):
        """
        Employs methodology from Pink and Zieglier (2016).
        """
        # Extremely unlikely edge case.
        if bound < 1:
            raise ValueError("bounds are already sufficient.")

        var('x')
        K = NumberField(x ** 2 - self.constants.A * x - self.constants.B, name='y')
        alpha = K.gen(0)
        beta = self.constants.A - alpha

        Z_bounds = []
        for i in range(len(self.constants.primes)):
            p = self.constants.primes[i]
            L = Qp(p)
            for t_vec in combinations(list(range(1, bound + 1)), self.constants.num_terms):
                tau = self.calculate_tau(t_vec, alpha, beta)

                print("tau: " + str(tau))
                z0_left_term = self.constants.a * (1 + sum([alpha ** t for t in t_vec]))
                z0_right_term = alpha - beta
                z0 = L(z0_left_term.norm()).ordp() - L(z0_right_term.norm()).ordp()

                # Unlikely case, and we may assume the p-adic log is injective for our purposes.
                if tau == 1:
                    p_ideal = K.primes_above(p)[0]
                    p_order = L(padic_log(alpha / beta, p_ideal, prec=50).norm()).ordp()
                    new_bound = math.log(self.coefficients.N_bounds[i], p) + p_order + z0
                    Z_bounds.append(max(z0 + 3/2, new_bound))
                else:
                    r = math.ceil(math.log(self.coefficients.N_bounds[i], p))
                    p_ideal = K.primes_above(p)[0]
                    zeta = padic_log(tau, p_ideal, prec=50) / padic_log(beta / alpha, p_ideal, prec=50)
                    zeta_list = zeta.list()
                    print(zeta_list) 
                    R_max = find_Rmax(50)
                    if R_max == -1:
                        raise ValueError("no max value of R found for " + str(t_vec))

        return Z_bounds

    def calculate_tau(self, t_vec, alpha, beta):
        numerator = 1 + sum([alpha ** t for t in t_vec])
        denominator = 1 + sum([beta ** t for t in t_vec])
        return numerator / denominator

    def find_Rmax(self, r, zeta_list, tries=50):
        Rmax = -1
        for i in range(r, r + tries):
            if zeta_list[i - vzeta] != None:
                break
        return Rmax

if __name__ == "__main__":
    constants_gen = Constants(
        a = 1,
        b = 1,
        A = 1,
        B = 1,
        alpha = (1 + math.sqrt(5))/2,
        beta = (1 - math.sqrt(5))/2,
        num_terms = 3,
        w = 1,
        primes = [2, 3, 5]
    )

    br = BoundReduce(constants_gen)
    br.padic_reduce(3000)