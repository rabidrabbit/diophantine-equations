"""
Solves the family of Diophantine equations of the form
u_{n_1} + ... + u_{n_k} = wp_1^{z_1} ... p_s^{z_s},
under some mild technical restrictions.

@authors: Brian Ha, Lily McBeath, Luisa Velasco
"""
import math
import numpy as np
from sage.all import *
from fpylll import IntegerMatrix, GSO, LLL

def LLL_real_reduce(bound, x_list):
    large_constant = 2
    # STEP 1: Find the LLL-reduced basis
    adjusted_bound = int(math.log10(bound)) + 1
    approximation_matrix = generate_lattice_approximation_matrix(adjusted_bound, x_list)

    L = LLL.reduction(IntegerMatrix.from_matrix(approximation_matrix))
    L_matrix = [[0 for _ in range(L.nrows)] for _ in range(L.ncols)]
    L.to_matrix(L_matrix)
    
    L_cols = get_column_vectors(L_matrix)

    # STEP 2: Perform the Gram-Schmidt process
    GS = GSO.Mat(L, flags=GSO.INT_GRAM); _ = GS.update_gso()
    GSM = GS.G
    GS_matrix = [[0 for _ in range(GSM.nrows)] for _ in range(GSM.ncols)]
    GSM.to_matrix(GS_matrix)

    GS_cols = get_column_vectors(GS_matrix)
    print(GS_cols)

    # STEP 3: Find c_2 and sigma
    # TODO: Fixed "adjusted_bound"
    vy = [0 for _ in range(adjusted_bound)]
    vy[-1] = -math.floor(large_constant * x_list[0])
    c2 = calculate_max_frac_norm_squared(L_cols, GS_cols)
    lat_min_dist_squared = calculate_lattice_min_distance_squared(L_cols, vy)
    sigma = calculate_sigma(L_matrix, vy)

    c1 = sigma * calculate_norm_squared(L_cols[0]) / c2
    print("value of c1: %f" % c1)

    if lat_min_dist_squared < c1:
        raise ValueError("lattice min distance squared inequality is invalid.")

    # STEP 4: Calculate values needed for Lemma 6
    S, T = calculate_S_and_T()

    # STEP 5: Find a new bound
    new_bound = (math.log(large_constant * c3) - math.log(math.sqrt((c_1 ** 2) - S) - T))
    if new_bound < 0:
        raise ValueError("new bound is a negative number.")
    return new_bound

def LLL_padic_reduce(test):
    return

def calculate_distance_to_nearest_int(x):
    """
    Used in calculating sigma.
    """
    y = abs(x)
    frac_part = y - int(y)
    if frac_part == 0:
        return ValueError("frac_part is 0.")
    return min(frac_part, 1 - frac_part)

def find_last_nonzero_index(v):
    """
    Used in calculating sigma.
    """
    last_index = math.inf
    for i in range(len(v)):
        if v[-(i + 1)] != 0:
            return len(v) - (i + 1)
    raise ValueError("Zero vector passed in.")

def calculate_sigma(matrix, v):
    """
    Calculate sigma as defined in the paper.
    """
    # STEP 1: Calculate the vector z
    np_inverted_matrix = np.linalg.inv(np.array(matrix))
    z = np.dot(np_inverted_matrix, v)

    # STEP 2: Find the largest index such that the entry is non-zero
    last_index = find_last_nonzero_index(z)

    # STEP 3: Calculate the distance to the nearest integer
    # TODO: Fix floating point errors.
    return calculate_distance_to_nearest_int(z[last_index])

def calculate_S_and_T(x_list):
    """
    Calculates S and T and returns it as a tuple (S, T)
    """
    # STEP 1: Calculate S
    S = sum([x ** 2 for x in x_list])
    
    # STEP 2: Calculate T
    T = ((sum(x_list) + c20)/2) ** 2
    
    return (S, T)

def calculate_norm_squared(v):
    return sum([x ** 2 for x in v])

def calculate_max_frac_norm_squared(lll_reduced_basis, gs_basis):
    return max([calculate_norm_squared(lll_reduced_basis[0])/calculate_norm_squared(v) for v in gs_basis])

def calculate_lattice_min_distance_squared(lll_reduced_basis, vy):
    """
    Calculates l(L, y)^2, where y is a vector.
    """
    distances = []
    for v in lll_reduced_basis:
        distances.append([x - y for x, y in zip(v, vy)])
    lattice_min_distance = min(map(calculate_norm_squared, distances))
    
    print("lattice min distance: %d" % lattice_min_distance)
    # TODO: Check if vy is in the lattice.
    return lattice_min_distance

def calculate_euclidean_distance_squared(vx, vy):
    """
    Calculates the Euclidean distance
    """
    return sum([(x - y) ** 2 for x, y in zip(vx, vy)])

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

    print("Calculating constants...")
    constants = calculate_constants(alpha, beta, a, b)
    print(constants)
    print("Done.")

    print("Handling n1 = ... = nk case...")
    print("Done.")

    print("Performing first reduction case...")
    LLL_real_reduce(100, [1, 2, 3])