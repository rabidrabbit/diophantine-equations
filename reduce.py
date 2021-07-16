from constants import Constant, log_p
from itertools import product


class BoundReduce:
    def __init__(self, constants):
        self.constants = constants
    
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
        beta = 2 * alpha - self.constants.A

        Z_bounds = []
        for i in range(len(self.constants.primes)):
            p = self.constants.primes[i]
            L = Qp(p)
            for t_vec in product(range(1, bound + 1), repeat = 3):
                tau = self.calculate_tau(t_vec)

                z0_left_term = self.constants.a * (1 + sum([alpha ** t for t in t_vec]))
                z0_right_term = alpha - beta
                z0 = L(z0_left_term).ordp() - L(z0_right_term).ordp()

                # Unlikely case, but we may assume the p-adic log is injective for our purposes.
                if tau == 1:
                    p_ideal = K.primes_above(p)[0]
                    p_order = L(padic_log(alpha / beta, p_ideal, prec=50).norm()).ordp()
                    new_bound = math.log(self.constants.N_bounds[i], p) + p_order + z0
                    Z_bounds.append(max(z0 + 3/2, new_bound))
                else:
                    r = math.ceil(math.log(self.constants.N_bounds[i], p))
                    
        return Z_bounds

    def calculate_tau(self, t_vec):
        print(t_vec)

if __name__ == "__main__":
    br = BoundReduce([1])
    br.padic_reduce(3)