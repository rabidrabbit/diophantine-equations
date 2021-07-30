"""
Class for reducing upper-bounds on the parameters of our Diophantine equations.

@authors: Brian Ha, Lily MacBeath, Luisa Velasco
"""


import math
import logging
from constants import Constants, padic_order
from itertools import combinations_with_replacement
from sage.all import *
from timeit import default_timer as timer

class BoundReduce:
    def __init__(self, constants, threshold=0.05, tries=50, flags = {}):
        """
        Initialization of the BoundReduce class.
        Threshold defines the minimum percentage reduction for another iteration of the reduction
        process to occur.
        """
        self.constants = constants
        self.coefficients = constants.calculate_constants()
        self.threshold = threshold
        self.tries = tries
        self.flags = flags

        self.initialize_flags()

    def initialize_flags(self):
        if self.flags["DEBUG_FLAG"]:
            logging.basicConfig(level=logging.INFO)
            logging.info("Verbose mode on.")
        
        self.flags["REAL_REDUCE_FLAG"] = True
        self.flags["PADIC_REDUCE_FLAG"] = True

    def print_summary(self):
        return

    def reduce(self, threshold):
        """
        Completely reduces the upper-bound.
        """
        def percentage_change(old, new):
            return (old - new) / old

        real_reduction_iterations = 0
        padic_reduction_iterations = 0
        cont_reduction_iterations = 0
        factor = len(self.constants.primes) + 5

        # First, go through the real reduction loop.
        current_n1_bound = self.coefficients['n1_bound']
        current_diff_bound = None
        while True:
            real_reduction_iterations += 1
            logging.info("Real Reduction - Iteration %d" % real_reduction_iterations)

            large_constant = self.calculate_large_constant(current_n1_bound, factor)
            logging.info("Large constant contains %d digits " % large_constant.ndigits())

            # Find a new bound on n_1 - n_k
            new_diff_bound = self.real_reduce(current_n1_bound, large_constant)

            logging.info("Current bound on n1: " + str(current_n1_bound))
            self.update_real_constants(new_diff_bound)
            logging.info("New bound on n1 - nk: " + str(new_diff_bound))
            logging.info("New bound on n1: " + str(self.coefficients["n1_bound"]))
            
            if percentage_change(current_n1_bound, self.coefficients["n1_bound"]) < self.threshold:
                logging.info("New bound did not improve in the real step; real reduction process is done.")
                factor = factor + 5
                break

            current_n1_bound = self.coefficients['n1_bound']
            current_diff_bound = new_diff_bound

        # Second, go through the p-adic reduction loop.
        current_Z_bounds = self.coefficients['Z_bounds']
        while True:
            padic_reduction_iterations += 1
            logging.info("p-adic Reduction - Iteration %d" % padic_reduction_iterations)

            new_Z_bounds = self.padic_reduce(ceil(current_diff_bound))
            logging.info("Current bound on n1: " + str(current_n1_bound))
            new_n1_bound = self.update_padic_constants(new_Z_bounds)
            logging.info("New bound on n1: " + str(new_n1_bound))
            if percentage_change(current_n1_bound, new_n1_bound) < self.threshold:
                logging.info("New bound did not improve in the p-adic step; p-adic reduction process is done.")
                break

            current_n1_bound = new_n1_bound
        logging.info("Final reduced bound on n1: " + str(current_n1_bound))

        return self.constants

    def update_real_constants(self, diff_bound):
        """
        Updates the constants given new bounds on n_1 - n_k.
        """
        self.constants.update_constants(diff_bound)
        self.coefficients = self.constants.get_coefficients()
        return

    def update_padic_constants(self, Z_bounds):
        n1_bound = self.constants.update_padic_constants(Z_bounds)
        self.coefficients = self.constants.get_coefficients()
        return n1_bound

    def calculate_large_constant(self, bound, factor):
        """
        Calculates an appropriate large constant given a bound.
        """
        minimum_exponent = ceil(math.log(bound, 10) * factor)
        return ZZ(10 ** minimum_exponent)

    def real_reduce(self, bound, large_constant):
        """
        Finds a better upper-bound on our values as defined in the procedure within our paper.
        """
        approximation_matrix = self.generate_approximation_matrix(large_constant)
        LLL_matrix = self.generate_LLL_matrix(approximation_matrix)
        GS_matrix = self.generate_GS_matrix(LLL_matrix) 
        minimal_vector_bound = self.calculate_minimal_vector_bound(large_constant, LLL_matrix, GS_matrix)
        S, T = self.calculate_S_and_T(self.coefficients["Z_bounds"])
        new_bound = self.real_reduce_new_bound(large_constant, minimal_vector_bound, S, T)
        return new_bound

    def generate_approximation_matrix(self, large_constant):
        """
        Generates approximation matrix needed for LLL computations.
        Note that we prepare the matrix assuming its columns generate the lattice.
        """
        n = len(self.constants.primes) + 1
        primes_row = [round(large_constant * log(p)) for p in self.constants.primes]
        primes_row.append(round(RR(large_constant) * log(self.constants.alpha)))
        approximation_matrix = []
        for i in range(n):
            zero_row = [0] * n
            zero_row[i] = 1
            approximation_matrix.append(zero_row)
        approximation_matrix[n - 1] = primes_row
        return Matrix(ZZ, approximation_matrix)

    def generate_LLL_matrix(self, matrix):
        """
        Assumes the input is formatted such that the column generates the lattice,
        and returns a matrix whose columns still generate the lattice.
        """
        LLL_matrix = matrix.transpose().LLL().transpose()
        return LLL_matrix

    def generate_GS_matrix(self, matrix):
        """
        Assumes the input is formatted such that the column generates the lattice,
        and returns a matrix whose columns still generate the lattice.
        """
        GS_matrix = Matrix(QQ, matrix.transpose().gram_schmidt()[0]).transpose()
        return GS_matrix

    def calculate_minimal_vector_bound(self, large_constant, LLL_matrix, GS_matrix, prec=100):
        """
        Calculates a lower-bound for the minimal vector.
        Names of the constants are named according to the paper.
        """
        R = RealField(prec)
        n = len(self.constants.primes) + 1

        # First, setup the vector vy.
        vy = [0 for _ in range(n)]
        eta_0 = R(self.constants.w * sqrt(self.constants.delta) / self.constants.a)
        vy[-1] = -R(large_constant) * R(math.log(eta_0))
        vy[-1] = vy[-1].floor()
        vy = vector(ZZ, vy)

        # Second, calculate the constants needed.
        sigma = self.calculate_sigma(LLL_matrix, vy, prec)
        c2 = max([LLL_matrix.column(0).norm()**2 / GS_matrix.column(i).norm()**2 for i in range(n)]) 

        # Lastly, calculate the lower-bound.
        minimal_vector_bound = (1 / c2) * sigma * LLL_matrix.column(0).norm()**2
        return minimal_vector_bound

    def calculate_sigma(self, LLL_matrix, vy, prec=100):
        """
        Calculates the value of sigma, which is the distance from the last non-zero entry
        in the passed-in vector to the nearest integer.

        A modified implementation of minimal_vector() from SAGE's S_unit_solver.
        """
        R = RealField(prec)
        z = (LLL_matrix.inverse()) * vy
        z_diff = []
        for elem in z:
            current_term = abs(R(elem - elem.round()))
            if current_term != 0:
                z_diff.append(current_term)
        return 1 if len(z_diff) == 0 else z_diff[-1]

    def calculate_S_and_T(self, z_bounds):
        """
        Calculate the values S and T to find a better bound on n_1 - n_k.
        Note that X_i = Z_i, where Z_i represents the upper-bound on z_i.
        """
        S = sum([RR(zb) ** 2 for zb in z_bounds])
        T = (self.coefficients['n1_bound'] + sum(z_bounds)) / 2 
        return (S, T)

    def real_reduce_new_bound(self, large_constant, minimal_vector_bound, S, T):
        """
        Calculates a new bound on (n_1 - n_k), given appropriate values.
        Names of constants correspond to the names defined in the paper.
        """
        c3 = RR(2 + 4 * self.constants.num_terms * abs(self.constants.b) / abs(self.constants.a))
        c4 = RR(math.log(min(abs(self.constants.alpha / self.constants.beta), self.constants.alpha))) 
        new_bound = (1 / c4) * (log(RR(large_constant) * c3) - log(sqrt(RR(minimal_vector_bound) - RR(S)) - RR(T)))
        return new_bound

    def padic_reduce(self, diff_bound):
        """
        Employs methodology from Pink and Zieglier (2016).
        """
        Z_bounds = []
        for i in range(len(self.constants.primes)):
            p = self.constants.primes[i]
            logging.info("Examining the prime %d for the p-adic reduction." % p)
            r = ceil(math.log(self.coefficients["n1_bound"], p))
            prec = r + self.tries + 10
            L = Qp(p, prec)

            current_z_bound = -1

            # Pre-computation step (i.e. memoization)
            alpha, beta = self.constants.calculate_alphabeta(p, prec)
            alpha_over_beta_ordp = (alpha / beta).ordp()
            alpha_t = self.precompute_t_power(alpha, diff_bound)
            beta_t = self.precompute_t_power(beta, diff_bound)
            z0_right_term = L((alpha - beta).norm()).ordp()
            log_alpha_over_beta = log(alpha / beta)

            # Pre-computations for tau == 1 
            p_order = log_alpha_over_beta.norm().ordp() / 2
            iterations = 0

            for t_vec in combinations_with_replacement(list(range(0, diff_bound + 1)), self.constants.num_terms - 1):
                iterations += 1
                if iterations % 1000 == 0:
                    print(iterations)

                z0_left_term = self.constants.a * (1 + sum([alpha_t[t] for t in t_vec]))
                z0 = L(z0_left_term.norm()).ordp() - z0_right_term
                tau = self.calculate_tau(t_vec, alpha_t, beta_t)

                # Unlikely case, and we may assume the p-adic log is injective for our purposes.
                if tau == 1:
                    new_bound = math.log(self.coefficients["n1_bound"], p) + p_order + z0
                    current_z_bound = max(z0 + 3/2, new_bound, current_z_bound)
                else:
                    # Remember that zeta is in Q_p, not in the field extension.
                    log_tau = tau.log(p_branch=0)
                    zeta = log_tau / log_alpha_over_beta
                    vzeta = zeta.ordp()
                    zeta_list = self.get_expansion(prec, zeta)

                    R_max = self.find_Rmax(r, vzeta, zeta_list)
                    if R_max == -1:
                        m0 = self.find_m0(p, zeta_list)
                        z_bound = math.log(self.coefficients["n1_bound"] - m0, p) + alpha_over_beta_ordp + z0
                    else:
                        zeta_inverse = log_alpha_over_beta / log_tau
                        vzeta_inverse = zeta_inverse.ordp()
                        z_bound = log_tau.ordp() + R_max + vzeta_inverse
                    current_z_bound = max(current_z_bound, z_bound)
            Z_bounds.append(current_z_bound)
        return Z_bounds

    def precompute_t_power(self, term, diff_bound):
        return [term ** i for i in range(0, diff_bound + 1)]

    def calculate_tau(self, t_vec, alpha_t, beta_t):
        numerator = 1 + sum([beta_t[t] for t in t_vec])
        denominator = 1 + sum([alpha_t[t] for t in t_vec])
        return (self.constants.b * numerator) / (self.constants.a * denominator)

    def get_expansion(self, prec, padic_num):
        """
        Obtains the expansion of a p-adic number.
        Necessary as SAGE's p-adic representations are different for Eisenstein and unramified extensions.
        """
        padic_expansion = list(padic_num.expansion())
        if isinstance(padic_expansion[0], list):
            return padic_expansion
        else:
            # Eistenstein extension case.
            padic_list = []
            for i in range(0, len(padic_expansion), 2):
                term = [padic_expansion[i]]
                padic_list.append(term)

            # Fill the rest of the list to the sufficient precision.
            for i in range(prec - len(padic_list)):
                padic_list.append([])    
            return padic_list

    def find_Rmax(self, r, vzeta, zeta_list, tries=50):
        Rmax = -1
        for i in range(r, r + tries):
            # Checks if i is a valid index
            if zeta_list[i - vzeta] != []:
                if i > Rmax:
                    Rmax = i
            elif i == r + tries - 1:
                logging.debug("No R found for Rmax...")
        return Rmax

    def find_m0(self, p, zeta_list):
        # In this case, it must be the case that zeta is an integer.
        m0 = 0
        for i in range(len(zeta_list)):
            if zeta_list[i] != []:
                m0 += (p ** i) * zeta_list[i][0]
        return m0


if __name__ == "__main__":
    start = timer()
    constants_gen = Constants(
        a = 1,
        b = 1,
        A = 2,
        B = 1,
        alpha = (2 + math.sqrt(8))/2,
        beta = (2 - math.sqrt(8))/2,
        delta = 8,
        num_terms = 4,
        w = 1,
        primes = [2, 3, 5]
    )

    br = BoundReduce(constants_gen, flags={"DEBUG_FLAG": True})
    br.reduce(threshold=0.01)
    end = timer()
    print(end - start)
