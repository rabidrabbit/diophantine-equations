"""
Class for reducing upper-bounds on the parameters of our Diophantine equations.

@authors: Brian Ha, Lily MacBeath, Luisa Velasco
"""


import math
import logging
from constants import Constants, padic_log, padic_order
from itertools import combinations
from sage.all import *

class BoundReduce:
    def __init__(self, constants, threshold = 0.05, tries=50, flags = {}):
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

    def print_summary(self):
        return

    def reduce(self, threshold):
        """
        Completely reduces the upper-bound.
        """
        def percentage_change(old, new):
            return (new - old) / old
        real_reduction_iterations = 0
        padic_reduction_iterations = 0
        cont_reduction_iterations = 0

        # First, go through the real reduction loop.
        current_n1_bound = self.coefficients['n1_bound']
        while True:
            real_reduction_iterations += 1
            logging.info("Real Reduction - Iteration %d" % real_reduction_iterations)

            large_constant = self.calculate_large_constant(current_n1_bound, len(self.constants.primes) + 1)
            logging.info("Large constant contains %d digits " % large_constant.ndigits())

            # Find a new bound on n_1 - n_k
            new_diff_bound = self.real_reduce(current_n1_bound, large_constant)
            logging.info("Current bound on n1: " + str(current_n1_bound))
            self.update_constants(new_diff_bound)
            logging.info("New bound on n1: " + str(self.coefficients["n1_bound"]))
            
            if (current_n1_bound - self.coefficients['n1_bound']) / self.coefficients['n1_bound'] < self.threshold:
                logging.info("New bound did not improve; real reduction process is done.")
                break

            current_n1_bound = self.coefficients['n1_bound']

        # Second, go through the p-adic reduction loop.

        return self.constants

    def update_constants(self, bound):
        """
        Updates the constants given new bounds on n_1 - n_k.
        """
        self.constants.update_constants(bound)
        self.coefficients = self.constants.get_coefficients()
        return

    def calculate_large_constant(self, bound, factor):
        """
        Calculates an appropriate large constant given a bound.
        """
        minimum_exponent = math.ceil(math.log(bound, 10) * factor)
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
        n = len(self.constants.primes)
        primes_row = [round(large_constant * log(p)) for p in self.constants.primes]
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
        return matrix.transpose().LLL().transpose()

    def generate_GS_matrix(self, matrix):
        """
        Assumes the input is formatted such that the column generates the lattice,
        and returns a matrix whose columns still generate the lattice.
        """
        return Matrix(QQ, matrix.transpose().gram_schmidt()[0]).transpose()

    def calculate_minimal_vector_bound(self, large_constant, LLL_matrix, GS_matrix, prec=100):
        """
        Calculates a lower-bound for the minimal vector.
        Names of the constants are named according to the paper.
        """
        R = RealField(prec)
        n = len(self.constants.primes)

        # First, setup the vector vy.
        vy = [0 for _ in range(n)]
        eta_0 = R(self.constants.w * sqrt(self.constants.delta) / self.constants.a)
        vy[-1] = -R(large_constant * eta_0).floor() 
        vy = vector(ZZ, vy)

        # Second, calculate the constants needed.
        sigma = self.calculate_sigma(LLL_matrix, vy, prec)
        c2 = max([LLL_matrix.column(i).norm()**2 / GS_matrix.column(i).norm()**2 for i in range(n)]) 

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
        S = sum([zb ** 2 for zb in z_bounds])
        T = (self.coefficients['n1_bound'] + sum(z_bounds)) / 2 
        return (S, T)

    def real_reduce_new_bound(self, large_constant, minimal_vector_bound, S, T):
        """
        Calculates a new bound on (n_1 - n_k), given appropriate values.
        Names of constants correspond to the names defined in the paper.
        """
        c3 = RR(2 + 2 * self.constants.num_terms * abs(self.constants.b) / abs(self.constants.a))
        c4 = RR(math.log(min(abs(self.constants.alpha / self.constants.beta), self.constants.alpha))) 
        new_bound = (1 / c4) * (log(RR(large_constant) * c3) - log(sqrt(RR(minimal_vector_bound)**2 - RR(S)) - RR(T)))
        return new_bound

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
        Returns (alpha, beta) as a tuple in p-adic representation.
        Attempts to get around NotImplementedError from SAGE by extending the
        field Qp in different ways (i.e. using some ideas from Hensel's lemma).
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
            except NotImplementedError:
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

    br = BoundReduce(constants_gen, flags={"DEBUG_FLAG": True})
    br.reduce(threshold=0.05)
    #br.padic_reduce(3000)
