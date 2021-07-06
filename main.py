"""
Solves the family of Diophantine equations of the form
u_{n_1} + ... + u_{n_k} = wp_1^{z_1} ... p_s^{z_s},
under some mild technical restrictions.

@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""
import math
from fpylll import IntegerMatrix, LLL

def LLL_real_reduce(bound, x_list):
    adjusted_bound = int(math.log10(bound)) + 1
    approximation_matrix = generate_lattice_approximation_matrix(adjusted_bound, x_list)

    print(approximation_matrix)

    L = LLL.reduction(IntegerMatrix.from_matrix(approximation_matrix))
    print(L.nrows)

    new_bound = 2
    return new_bound

def LLL_padic_reduce(test):
    return

def calculate_lattice_min_distance(col_vectors, vy):
    """
    Calculates l(L, y), where y is a vector.
    """
    lattice_min_distance = math.inf

    # First, verify if v_y is in the lattice
    for c in col_vectors:
        print(c)

    print("lattice min distance: " % lattice_min_distance)

def calculate_euclidean_distance(vx, vy):
    """
    Calculates the Euclidean distance
    """

def get_column_vectors(matrix):
    """
    Extract the column vectors from a matrix (in list representation)
    """
    column_vectors = []
    for i in range(len(matrix)):
        col = []
        for j in range(len(matrix)):
            col.append(matrix[j][i])
        column_vectors.append(col)
    return column_vectors

def calculate_new_bound():
    return

def calculate_d1():
    return

def calculate_constants(alpha, beta, a, b):
    """
    Calculates constants that are defined within our paper.
    """
    constants = {}
    constants["c1"] = 2
    constants["c2"] = 2
    constants["c3"] = 2
    constants["c4"] = 2
    constants["c5"] = 2

    constants["d1"] = 2

    return constants

def generate_lattice_approximation_matrix(adjusted_bound, x_list):
    """
    Generates the lattice approximation matrix.
    """
    size = len(x_list)
    approximation_matrix = [[0 for _ in range(size - 1)] for _ in range(size - 1)]
    IntegerMatrix.identity(size - 1).to_matrix(approximation_matrix)

    # Append an empty column to the identity matrix
    for row in approximation_matrix:
        row.append(0)

    # Create the row that represents the approximation of the linear form
    row_approx = []
    for i in range(len(x_list)):
        row_approx.append(round(adjusted_bound * x_list[i]))
    approximation_matrix.append(row_approx)
    return approximation_matrix

if __name__ == "__main__":
    alpha = (1 + math.sqrt(5))/2
    beta = (1 - math.sqrt(5))/2
    a = 1
    b = 1

    constants = calculate_constants(alpha, beta, a, b)
    print(constants)

    LLL_real_reduce(100, [1, 2, 3])