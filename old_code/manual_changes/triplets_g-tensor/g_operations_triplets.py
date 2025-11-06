"""
    G-TENSOR CALCULATION WITH ORBITAL AND SPIN
    ANGULAR MOMENTUM AND SOC BETWEEN STATES
 Analysis of the excited states and print of those
 orbitals involved in the configurations with the
 highest amplitudes in the excited states
"""
import numpy as np
from numpy import linalg, sqrt
import matplotlib.pyplot as plt

def get_Hamiltonian_construction(selected_states, eigenenergies, spin_orbit_coupling):
    """
    Hamiltonian is written 'bra' in rows and 'ket' in columns,
    with spin order -1/2 , +1/2.
    :param: selected_states, eigenenergies, spin_orbit_coupling
    :return: Hamiltonian
    """
    Hamiltonian = np.zeros((len(selected_states) * 2, len(selected_states) * 2), dtype=complex)

    for i in range(0, len(selected_states) * 2):
        for j in range(0, len(selected_states) * 2):
            if (i == j):
                Hamiltonian[i, i] = eigenenergies[i//2]
            else:
                Hamiltonian[i, j] = spin_orbit_coupling[i,j]

    # print('Hamiltonian (SOC in cm-1, energies in a.u):')
    # print('\n'.join([''.join(['{:^20}'.format(item) for item in row])\
    #                  for row in np.round((Hamiltonian[:,:]),5)]))
    # print(" ")
    # exit()

    return Hamiltonian

def Hamiltonian_diagonalization( Hamiltonian ):
    """
    1) Hamiltonian is diagonalized
    2) eigenvectors-eigenvalues are ordered by weight coefficients
    3) Doublet with the lowest energy is set as the new basis (Kramer doublet)
    :param: Hamiltonian
    :return eigenvalues, eigenvectors, kramer_st: kramer_st is the index
    of the state set as Kramer doublet * 2
    """
    eigenvalues, eigenvectors = linalg.eigh(Hamiltonian)

    # Reorder eigenvectors (and eigenenergies) by weight coefficients
    change_order = np.zeros( len(eigenvectors) , dtype=complex)
    for v_1 in range(0, len(eigenvectors)):
        for v_2 in range(v_1, len(eigenvectors)):

            if ( abs(eigenvectors[v_1, v_2]) > abs(eigenvectors[v_1, v_1]) ):

                change_order[:] = eigenvectors[:,v_1]
                eigenvectors[:,v_1] = eigenvectors[:,v_2]
                eigenvectors[:,v_2] = change_order[:]

                change_order.real[0] = eigenvalues[v_1]
                eigenvalues[v_1] = eigenvalues[v_2]
                eigenvalues[v_2] = change_order.real[0]

    # Kramer doublets selection:
    minimum_energy = min(eigenvalues)
    eigenvalues_list = list(eigenvalues)
    kramer_st = eigenvalues_list.index( minimum_energy )

    # The index of the selected state must be even since
    # Kramer doublets are [kramer_st, kramer_st+1]
    if ( kramer_st % 2 ) != 0: kramer_st = kramer_st - 1
    else: pass

    return eigenvalues, eigenvectors, kramer_st

def angular_matrixes_obtention( eigenvalues, eigenvectors, kramer_st, input_angular_matrix ):
    """
    Spin or orbital angular matrix calculation using:
    1) coeff_bra, coeff_ket: coefficients of the lineal combination of non-relativistic states,
    that come from Kramer doublet states eigenvectors
    2) angular_value: angular momentum between states. Depending on the
    column of the final matrix, it takes real (col 0), imaginary (col 1) or
    both parts (col 2).

    :param eigenvalues, eigenvectors, kramer_st, input_angular_matrix
    :return angular_matrix: contains the spin value < B(S,Sz) | Sx | A(S',Sz') >,
    meaning < B(S,Sz)| corresponds to state "i" (in rows) and | A(S',Sz') > to
    state "j" (in columns)
    """
    # Matrices calculation:
    angular_matrix = np.zeros((3, 3), dtype=complex)

    for row in range(0, 3):  # dimension x,y,z
        for column in range(0, 3):  # dimension x,y,z

            for bra in range(0, len(eigenvalues) ): # state <B|
                for ket in range(0, len(eigenvalues) ): # state |A>

                    # coeff_ket for 1st and 2nd column ("x", "y" dir)
                    # coeff_ket_2 for 3rd column ("z" direction)
                    coeff_bra = np.conj( eigenvectors[bra, kramer_st + 1 ] )
                    coeff_ket = ( eigenvectors[ ket, kramer_st ] )
                    coeff_ket_2 = ( eigenvectors[ket, kramer_st + 1] )

                    angular_value = ( input_angular_matrix[bra, ket, row] )

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

def g_factor_calculation( lambda_matrix, sigma_matrix ):
    """
    Calculation of the G-tensor with lambda and sigma matrices. Then, g-factors
    are calculated as square roots of the eigenvalues of the G-tensor.
    """

    # G-tensor matrix obtention:
    lande_factor = 2.002319304363
    sigma_plus_lambda = np.zeros((3, 3), dtype=complex)
    sigma_plus_lambda = lande_factor * sigma_matrix + lambda_matrix

    G_matrix = np.zeros((3, 3), dtype=complex)
    G_matrix = np.matmul((sigma_plus_lambda), np.transpose(sigma_plus_lambda))

    # Diagonalize and reorder by weight coefficients:
    G_matrix_diagonal, rotation_matrix = linalg.eigh(G_matrix)
    change_order = np.zeros( len(G_matrix_diagonal) , dtype=complex)

    for i in range(0, 3):
        for j in range(i,3):
            if ( abs(rotation_matrix[i,j]) > abs(rotation_matrix[i,i]) ):

                change_order[:] = rotation_matrix[:,j]
                rotation_matrix[:,j] = rotation_matrix[:,i]
                rotation_matrix[:,i] = change_order[:]

                change_order.real[0] = G_matrix_diagonal[j]
                G_matrix_diagonal[j] = G_matrix_diagonal[i]
                G_matrix_diagonal[i] = change_order.real[0]

    # Obtain the g-factor: g = O^(1r) * X * sqrt(Gdiag) * O
    X_mat = np.identity(3)
    g_value = np.zeros(3, dtype=complex)
    g_value = (np.transpose(rotation_matrix)).dot(X_mat).dot(sqrt(G_matrix_diagonal)).dot(rotation_matrix)
    g_value = ( (g_value) - lande_factor ) * 1000

    G_tensor_results = np.zeros(3, dtype=complex)
    for i in range(0,3):
        G_tensor_results[i] = ( sqrt( G_matrix_diagonal[i] ) - lande_factor ) * 1000

    return G_matrix, G_tensor_results

# def state_by_state_g_factor(input_file, totalstates, selected_SOC):
#     """
#     G-tensor procedure is here applied for all the states independently. G-tensor values for each state
#     in each dimension are saved to list and then matrix:
#
#     - sum_states_g_factor_x --> sum_g_x
#     - sum_states_g_factor_y --> sum_g_y
#     - sum_states_g_factor_z --> sum_g_z
#
#     - state_by_state_gvalue: contains all previous matrices results. Its dimensions are (len(sum_g_x), 3).
#     - sum_over_states_gvalue: contains total sum of all g-values in each dimension
#
#     Lastly, a presentation matrix is done to present the results.
#     """
#
#     from g_read import get_number_of_states, get_ras_eigenenergies, get_casscf_eigenenergies, get_sfdft_eigenenergies,\
#         get_spin_orbit_couplings, get_spin_of_each_state, get_spin_matrices, get_orbital_matrices
#
#     from g_operations import get_Hamiltonian_construction, Hamiltonian_diagonalization, \
#         angular_matrixes_obtention, g_factor_calculation
#
#     each_state_g_factor_x = []
#     each_state_g_factor_y = []
#     each_state_g_factor_z = []
#     state_energy_list = []
#
#     for state_number in [ i for i in range(1,totalstates+1) if i != 1 ]:
#         states_ras = [1, state_number]
#
#         eigenenergies = get_ras_eigenenergies(states_ras, input_file)
#
#         spin_orbit_coupling = get_spin_orbit_couplings(selected_SOC, totalstates, states_ras, input_file)
#
#         Ham_ras = get_Hamiltonian_construction(states_ras, eigenenergies, spin_orbit_coupling)
#
#         la_Ham_ras, v_Ham_ras, kramer_state = Hamiltonian_diagonalization(Ham_ras)
#
#         state_spin = get_spin_of_each_state(input_file)
#
#         spin_matrix = get_spin_matrices(input_file, states_ras, state_spin)
#         l_matrix = get_orbital_matrices(input_file, totalstates, states_ras)
#
#         lambda_matrix = angular_matrixes_obtention(la_Ham_ras, v_Ham_ras, kramer_state, l_matrix)
#         sigma_matrix = angular_matrixes_obtention(la_Ham_ras, v_Ham_ras, kramer_state, spin_matrix)
#
#         G_tensor, G_tensor_results = g_factor_calculation(lambda_matrix, sigma_matrix)
#
#         each_state_g_factor_x.append(G_tensor_results.real[0])
#         each_state_g_factor_y.append(G_tensor_results.real[1])
#         each_state_g_factor_z.append(G_tensor_results.real[2])
#
#         if state_number == 2: state_energy_list.append( float(eigenenergies[0]) )
#         state_energy_list.append( float(eigenenergies[2]))
#
#     all_g_x = np.array(each_state_g_factor_x, dtype=float)
#     all_g_y = np.array(each_state_g_factor_y, dtype=float)
#     all_g_z = np.array(each_state_g_factor_z, dtype=float)
#     state_energy = np.array(state_energy_list, dtype=float)
#
#     # With all g_values calculated independently for each state in "all_g_x", "all_g_y" and "all_g_z", these
#     # values are saved to "state_by_state_gvalue" and in each dimension values are summed to "sum_over_states_gvalue"
#
#     state_by_state_gvalues = np.zeros( ( len(all_g_x), 3) , dtype=complex)
#     total_sum_gvalues = np.zeros( 3 , dtype=complex)
#
#     n_dim = 0
#
#     for sum_g in [all_g_x, all_g_y, all_g_z]:
#
#         for sum_g_index in range(0, len(sum_g)):
#
#             state_by_state_gvalues[sum_g_index, n_dim] = sum_g[sum_g_index]
#             total_sum_gvalues[n_dim] += sum_g[sum_g_index]
#
#         n_dim += 1
#
#     # Prepare the matrix to show resutls by order
#
#     presentation_list = []
#     presentation_list.append(['state', 'x-dim', 'y-dim', 'z-dim', 'energies (cm-1)'])
#
#     for i in range(0, totalstates - 1):  # "totalstates-1" g-values (not account ground state)
#
#         x_result = np.round(state_by_state_gvalues.real[i, 0], 4)
#         y_result = np.round(state_by_state_gvalues.real[i, 1], 4)
#         z_result = np.round(state_by_state_gvalues.real[i, 2], 4)
#
#         energy_result = ( state_energy[ i+1 ] - state_energy[0] ) * 219474.63
#         energy = np.round(energy_result, 4)
#
#         presentation_list.append([ i+2, x_result, y_result, z_result, energy])
#
#     x_total_result = np.round(total_sum_gvalues[0].real, 4)
#     y_total_result = np.round(total_sum_gvalues[1].real, 4)
#     z_total_result = np.round(total_sum_gvalues[2].real, 4)
#     presentation_list.append(['Total:', x_total_result, y_total_result, z_total_result, '--'])
#
#     presentation_matrix = np.array(presentation_list, dtype=object)
#
#     return state_by_state_gvalues, total_sum_gvalues, presentation_matrix

def perturbative_method(totalstates, eigenenergies, spin_orbit_coupling):
    """
    The coefficients that previously were on the eigenvectors
    are now calculated using equation 27: they are called "eta". Similar
    process to the one used in Hamiltonian construction, BUT considering
    also the eigenvalues substraction in the denominator: en[a] - en[b],
    where the "a" corresponds to the row number and the "b" to the column number
    """
    coef = np.zeros((len(states_ras) * 2, len(states_ras) * 2), dtype=complex)
    k1 = 0
    k2 = 0
    for y in range(0, len(totalstates) * 2):
        k1 = 0
        for x in range(0, len(totalstates) * 2):
            if (y == x) or (x % 2 != 0 and x - 1 == y):
                coef[x][y] = 0
                k1 = k1 + 1
            elif (y % 2 == 0):
                coef[x][y] = spin_orbit_coupling[x - k1 + k2][0] / (eigenenergies[x] - eigenenergies[y])
                coef[x][y + 1] = spin_orbit_coupling[x - k1 + k2][1] / (eigenenergies[x] - eigenenergies[y + 1])
        k2 = k2 + len(states_ras) - 1

    # The procedure is the same than for the Hamiltonian case. Coefficients matrix
    # obtained by perturbation expressions is symmetric, meaning multiplication does not
    # change if it is done by columns or by rows
    # gpt_mat_ras = np.zeros((len(states_ras) * 2, len(states_ras) * 2, 3), dtype=complex)
    # for z in range(0, 3):
    #     for y in range(0, len(states_ras) * 2, 2):
    #         for x in range(0, len(states_ras) * 2):
    #             if (y == x) or (x % 2 != 0 and x - 1 == y):
    #                 gpt_mat_ras[y, x, z] = 0
    #                 gpt_mat_ras[x, y, z] = 0
    #             elif (x % 2 == 0):
    #                 gpt_mat_ras[x, y, z] = -4000 * etacoeff_ras[y, x] * l_matrix[x // 2, y // 2, z]
    #                 gpt_mat_ras[x + 1, y, z] = -4000 * etacoeff_ras[y, x + 1] * l_matrix[x // 2, y // 2, z]
    #                 gpt_mat_ras[x, y + 1, z] = -4000 * etacoeff_ras[y + 1, x] * l_matrix[x // 2, y // 2, z]
    #                 gpt_mat_ras[x + 1, y + 1, z] = -4000 * etacoeff_ras[y + 1, x + 1] * l_matrix[x // 2, y // 2, z]

    return

def plot_obtention(input, x_data, y_data):
    # "matplotlib" help: https://aprendeconalf.es/docencia/python/manual/matplotlib/
    # https://chartio.com/resources/tutorials/how-to-save-a-plot-to-a-file-using-matplotlib/
    # https://pythonspot.com/matplotlib-bar-chart/

    fuente = "Sans"
    size = 12

    y_pos = np.arange(len(x_data))
    plt.bar(y_pos, y_data, align='center', alpha=0.5, color='red')
    plt.xticks(y_pos, x_data)

    plt.title(input, fontsize=size, fontname=fuente)
    plt.ylabel('Energies (eV)', fontsize=size, fontname=fuente)
    plt.xlabel('number of state', fontsize=size, fontname=fuente)
    plt.axis([min(x_data) - 2, max(x_data), min(y_data), max(y_data) + 0.1 * max(y_data)])
    # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    # plt.grid(True)

    plt.plot()
    figure_name = input + '.png'
    plt.savefig(figure_name)

    # plt.show()
    plt.close()
