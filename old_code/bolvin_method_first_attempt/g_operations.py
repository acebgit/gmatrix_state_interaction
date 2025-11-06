import numpy as np
from numpy import linalg
from numpy import sqrt


def get_Hamiltonian_construction( states_ras, eigenenergies, spin_orbit_coupling ):
    """
    In Hamiltonian, SOC values are written with < B | in rows and | A > in columns.
    For example, element Hamiltonian[4,0] is < state 2 (1/2, +1/2) | Hsomf | state 1 (1/2, -1/2) >
    """

    hamiltonian = np.zeros( (len(states_ras) * 2, len(states_ras) * 2), dtype=complex)
    for i in range(0, len(states_ras) * 2):
        for j in range(0, len(states_ras) * 2):
            if i == j:
                hamiltonian[i, i] = eigenenergies[i]
            else:
                hamiltonian[i, j] = spin_orbit_coupling[j, i]

    for i in range(0, len(hamiltonian)):
        for j in range(i, len(hamiltonian)):

            if hamiltonian[i, j] != np.conj(hamiltonian[j, i]):

                print("Hamiltonian is not Hermitian. Errors between states:", i // 2 + 1, j // 2 + 1)
                print(hamiltonian[i, j], hamiltonian[j, i])
                exit()

    return hamiltonian


def hamiltonian_diagonalization(hamiltonian):

    eigenvalues, eigenvectors = linalg.eigh(hamiltonian)

    return eigenvalues, eigenvectors


def angular_matrices_obtention(eigenvalues, eigenvectors, input_angular_matrix):
    """
    First, doublets with the lowest energy are selected to be the basis set (Kramer doublets)

    Secondly, angular matrix is calculated with:
    - ci and cj: coefficients of the lineal combination of initial states.
    - spin_matrix_value: contains the spin value < B(S,Sz) | Sx | A(S',Sz') >, meaning < B(S,Sz)| corresponds to
                        state "i" (in rows) and | A(S',Sz') > to state "j" (in columns)
    """

    # Kramer doublets selection:
    minimum_energy = min(eigenvalues)
    eigenvalues_list = list(eigenvalues)
    kramer_st = eigenvalues_list.index( minimum_energy )

    # The index of the selected state must be even since Kramer doublets are [kramer_st, kramer_st+1]
    if kramer_st % 2 == 0:
        kramer_st = kramer_st
    else:
        kramer_st = kramer_st - 1

    # Matrices calculation:
    angular_matrix = np.zeros((3, 3), dtype=complex)

    for row in range(0, 3):  # For each dimension x,y,z
        for column in range(0, 3):  # For each dimension x,y,z

            for state_i in range(0, len(eigenvalues) ):
                for state_j in range(0, len(eigenvalues) ):

                    coeff_ci = np.conj( eigenvectors[state_i, kramer_st + 1])
                    coeff_cj = ( eigenvectors[state_j, kramer_st])

                    coeff_ci_2 = np.conj( eigenvectors[state_i, kramer_st + 1])
                    coeff_cj_2 = ( eigenvectors[state_j, kramer_st + 1])

                    angular_value = ( input_angular_matrix[state_i, state_j, row])

                    if column == 0:
                        element = coeff_ci * coeff_cj * angular_value
                        angular_matrix[row, column] = angular_matrix[row, column] + 2 * element.real

                    elif column == 1:
                        element = coeff_ci * coeff_cj * angular_value
                        angular_matrix[row, column] = angular_matrix[row, column] + 2 * element.imag

                    elif column == 2:
                        element = coeff_ci_2 * coeff_cj_2 * angular_value
                        angular_matrix[row, column] = angular_matrix[row, column] + 2 * element

    # Reorder the columns by the weight factor (maximums must be in the diagonal)
    change_columns = np.zeros(3, dtype=complex)
    for column in range(0, 2):  # For each dimension x,y,z

        if abs(angular_matrix[column, column]) < abs(angular_matrix[column + 1, column]):
            change_columns[:] = angular_matrix[:, column]
            angular_matrix[:, column] = angular_matrix[:, column + 1]
            angular_matrix[:, column + 1] = change_columns[:]
        
    return angular_matrix


def g_factor_calculation(lambda_matrix, sigma_matrix):

    lande_factor = 2.002319

    # G-tensor matrix:
    g_tensor = np.zeros((3, 3), dtype=complex)

    for row in range(0, 3):
        for column in range(0, 3):
            for third in range(0, 3):

                row_result = lande_factor * (sigma_matrix[row, third]) + (lambda_matrix[row, third])
                column_result = lande_factor * (sigma_matrix[column, third]) + (lambda_matrix[column, third])

                g_tensor[row, column] = g_tensor[row, column] + row_result * column_result

    # G-tensor eigenvalues computation and reorder by weight coefficients to obtain well the dimensions
    g_tensor_eigenvalue, g_tensor_eigenvector = linalg.eigh(g_tensor)

    change_order = np.zeros( len(g_tensor_eigenvalue) , dtype=complex)

    for i in range(0, 3):
        for j in range(i,3):

            if abs(g_tensor_eigenvector[i,j]) > abs(g_tensor_eigenvector[i,i]):

                change_order[:] = g_tensor_eigenvector[:, j]
                g_tensor_eigenvector[:, j] = g_tensor_eigenvector[:, i]
                g_tensor_eigenvector[:, i] = change_order[:]

                change_order.real[0] = g_tensor_eigenvalue[j]
                g_tensor_eigenvalue[j] = g_tensor_eigenvalue[i]
                g_tensor_eigenvalue[i] = change_order.real[0]

    # Obtain the G-tensor square root eigenvalues
    g_tensor_results = np.zeros( 3, dtype=complex)

    for i in range(0,3):
        g_tensor_results[i] = (sqrt( g_tensor_eigenvalue[i]) - lande_factor) * 1000

    return g_tensor, g_tensor_results
