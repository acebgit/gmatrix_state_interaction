import numpy as np
from numpy import linalg, sqrt


def hermitian_test(matrix):
    """
    Check if a matrix is Hermitian. If not, exit.
    :param: matrix
    """
    for i in range(0, len(matrix)):
        for j in range(i, len(matrix)):
            element_1 = np.round(matrix[i, j], 4)
            element_2 = np.round(np.conjugate(matrix[j, i]), 4)
            if element_1 != element_2:
                print("Hamiltonian is not Hermitian")
                print('positions: ', i // 2, 'value:', matrix[i, j])
                print('positions: ', j // 2, 'value:', matrix[j, i])
                exit()
                
                
def reordering_eigenvectors(eigenval, eigenvect):
    """
    Reorder eigenvectors (and eigenenergies) by weight coefficients
    :param: eigenvalues, eigenvectors
    :return: eigenvalues, eigenvectors
    """
    change_order = np.zeros(len(eigenvect), dtype=complex)

    for v_1 in range(0, len(eigenvect)):
        for v_2 in range(v_1, len(eigenvect)):

            if abs(eigenvect[v_1, v_2]) > abs(eigenvect[v_1, v_1]):
                change_order[:] = eigenvect[:, v_1]
                eigenvect[:, v_1] = eigenvect[:, v_2]
                eigenvect[:, v_2] = change_order[:]

                change_order.real[0] = eigenval[v_1]
                eigenval[v_1] = eigenval[v_2]
                eigenval[v_2] = change_order.real[0]
    return eigenval, eigenvect


def get_hamiltonian_construction(selected_states, eigenenergies, spin_orbit_coupling, sz_values):
    """
    Hamiltonian is written 'bra' in rows and 'ket' in columns, with spin order -1/2 , +1/2.
    :param: selected_states, eigenenergies, spin_orbit_coupling
    :return: hamiltonian
    """
    hamiltonian = np.zeros((len(selected_states) * len(sz_values),
                            len(selected_states) * len(sz_values)), dtype=complex)

    for i in range(0, len(selected_states) * len(sz_values)):
        for j in range(0, len(selected_states) * len(sz_values)):
            if i == j:
                hamiltonian[i, i] = eigenenergies[i // len(sz_values)]
            else:
                hamiltonian[i, j] = spin_orbit_coupling[i, j]

    hermitian_test(hamiltonian)
    return hamiltonian


def hamiltonian_diagonalization(hamiltonian):
    """
    1) Hamiltonian is diagonalized
    2) eigenvectors-eigenvalues are ordered by weight coefficients
    3) Doublet with the lowest energy is set as the new basis (Kramer doublet)
    :param: Hamiltonian
    :return: eigenvalues, eigenvectors, kramer_st: kramer_st is the index
    of the state set as Kramer doublet * 2
    """
    eigenvalues, eigenvectors = linalg.eigh(hamiltonian)
    eigenvalues, eigenvectors = reordering_eigenvectors(eigenvalues, eigenvectors)

    # Kramer doublets selection:
    minimum_energy = min(eigenvalues)
    eigenvalues_list = list(eigenvalues)
    kramer_st = eigenvalues_list.index(minimum_energy)

    # The index of the selected state must be even since Kramer doublets are [kramer_st, kramer_st+1]
    if (kramer_st % 2) != 0:
        kramer_st = kramer_st - 1
    return eigenvalues, eigenvectors, kramer_st


def angular_matrixes_obtention(eigenvalues, eigenvectors, kramer_st, input_angular_matrix):
    """
    Spin or orbital angular matrix calculation using:
    1) coeff_bra, coeff_ket: coefficients of the lineal combination of non-relativistic states,
    that come from Kramer doublet states eigenvectors
    2) angular_value: angular momentum between states. Depending on the column of the final matrix,
    it takes real (col 0), imaginary (col 1) or both parts (col 2).
    :param: eigenvalues, eigenvectors, kramer_st, input_angular_matrix
    :return: angular_matrix: contains the spin value < B(S,Sz) | Sx | A(S',Sz') >
    """
    angular_matrix = np.zeros((3, 3), dtype=complex)

    for row in range(0, 3):  # dimension x,y,z
        for column in range(0, 3):  # dimension x,y,z

            for bra in range(0, len(eigenvalues)):  # state <B|
                for ket in range(0, len(eigenvalues)):  # state |A>

                    coeff_bra = np.conj(eigenvectors[bra, kramer_st + 1])
                    coeff_ket = (eigenvectors[ket, kramer_st])
                    coeff_ket_2 = (eigenvectors[ket, kramer_st + 1])
                    angular_value = (input_angular_matrix[bra, ket, row])

                    if column == 0:
                        element = coeff_bra * coeff_ket * angular_value
                        angular_matrix[row, column] += 2 * element.real

                    elif column == 1:
                        element = coeff_bra * coeff_ket * angular_value
                        angular_matrix[row, column] += 2 * element.imag

                    elif column == 2:
                        element = coeff_bra * coeff_ket_2 * angular_value
                        angular_matrix[row, column] += 2 * element

    # print('SIGMA matrix with all spin angular momentums:')
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                  for row in np.round((angular_matrix[:,:]),8)]))
    # print(" ")
    return angular_matrix


def g_factor_calculation(lambda_matrix, sigma_matrix):
    """
    Calculation of the G-tensor with lambda and sigma matrices. Then, g-factors
    are calculated as square roots of the eigenvalues of the G-tensor.
    :param: lambda_matrix, sigma_matrix
    :return: upper_g_matrix, g_tensor_values
    """
    # G-tensor matrix obtention:
    lande_factor = 2.002319304363
    sigma_plus_lambda = lande_factor * sigma_matrix + lambda_matrix

    # Diagonalize and reorder by weight coefficients:
    upper_g_matrix = np.matmul(sigma_plus_lambda, np.transpose(sigma_plus_lambda))
    upper_g_matrix_diagonal, rotation_matrix = linalg.eigh(upper_g_matrix)

    change_order = np.zeros(len(upper_g_matrix_diagonal), dtype=complex)
    for i in range(0, 3):
        for j in range(i, 3):
            if abs(rotation_matrix[i, j]) > abs(rotation_matrix[i, i]):
                change_order[:] = rotation_matrix[:, j]
                rotation_matrix[:, j] = rotation_matrix[:, i]
                rotation_matrix[:, i] = change_order[:]

                change_order.real[0] = upper_g_matrix_diagonal[j]
                upper_g_matrix_diagonal[j] = upper_g_matrix_diagonal[i]
                upper_g_matrix_diagonal[i] = change_order.real[0]

    g_tensor_values = np.zeros(3, dtype=complex)
    for i in range(0, 3):
        g_tensor_values[i] = (sqrt(upper_g_matrix_diagonal[i]) - lande_factor) * 1000
    return upper_g_matrix, g_tensor_values


def from_energies_soc_to_g_values(file, states_ras, totalstates,
                                  excitation_energies_ras, soc_ras, sz_list):
    """"
    Obtention of the g-values from the eigenenergies and the SOCs.
    :param: file, states_ras, totalstates, excitation_energies_ras, soc_ras, sz_list
    :return: upper_g_matrix, g_values
    """
    from parser_init import get_spin_matrices, get_orbital_matrices

    hamiltonian_ras = get_hamiltonian_construction(states_ras, excitation_energies_ras, soc_ras, sz_list)

    eigenvalues, eigenvector, kramers_states = hamiltonian_diagonalization(hamiltonian_ras)

    spin_matrix, standard_spin_matrix = get_spin_matrices(file, states_ras, bolvin=1)

    l_matrix = get_orbital_matrices(file, totalstates, states_ras, sz_list, bolvin=1)

    sigma_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, spin_matrix)

    lambda_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, l_matrix)

    upper_g_matrix, g_values = g_factor_calculation(lambda_matrix, sigma_matrix)

    return upper_g_matrix, g_values


def print_g_calculation(file, totalstates, selected_states, symmetry_selection,
                        states_ras, upper_g_tensor_results_ras):
    print("--------------------------------------")
    print("     INPUT SECTION")
    print("--------------------------------------")
    print("File selected: ", file)
    print("Number of states: ", totalstates)
    if selected_states == 2:
        print("Symmetry: ", symmetry_selection)
        print("Selected states: ", states_ras)
    else:
        print("Selected states: ", states_ras)

    print(" ")
    print("------------------------")
    print(" ras-CI OUTPUT SECTION")
    print("------------------------")
    print('g-factor (x y z dimensions):')
    print(np.round(upper_g_tensor_results_ras.real[0], 3), np.round(upper_g_tensor_results_ras.real[1], 3),
          np.round(upper_g_tensor_results_ras.real[2], 3))
    print('')
