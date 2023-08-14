__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np
from numpy import linalg, sqrt
from eom_ras_program.parser_rasci_pyqchem import parser_rasci


def get_number_of_states(file):
    """
    Obtain the total number of states in ras
    :param: file
    :return: totalstates
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    element = 0
    word_search = ['Requested states: ']

    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[3]
            break
    totalstates = int(element)
    states = list(range(1, totalstates + 1))
    return totalstates, states


def get_eigenenergies(file, totalstates, selected_states):
    """
    Get energies of the RAS-CI states.
    :param: file, totalstates, selected_states
    :return: eigenenergies, excitation_energies
    """
    word_search = ' RAS-CI total energy for state  '
    element_energy = []
    elements_excitenergy = []

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                line = line.split()
                element = line[6]
                element_energy.append(element)

                next_line = next(file)
                next_line = next_line.split()
                element = next_line[4]
                elements_excitenergy.append(element)
            if len(element_energy) == totalstates:
                break

    energies_selected = []
    excited_energies_selected = []
    for i in selected_states:
        energies_selected.append(element_energy[i - 1])
        excited_energies_selected.append(elements_excitenergy[i - 1])

    eigenenergies = np.array(energies_selected, dtype=float)
    excitation_energies_ev = np.array(excited_energies_selected, dtype=float)
    excitation_energies = excitation_energies_ev / 27.211399  # From eV to a.u.
    return eigenenergies, excitation_energies


def get_spin_orbit_couplings(file, totalstates, selected_states):
    """
    Spin-orbit coupling values are written in matrix with 'bra' in rows
    and 'ket' in columns, with spin order -Ms , +Ms.
    :param: file, totalstates, selected_states
    :return: doublet_soc, sz_list
    """
    with open(file, encoding="utf8") as f:
        output = f.read()
    output = parser_rasci(output)
    data = output['interstate_properties']

    def get_states_sz(qchem_file):
        """
        Get S² and Sz of all states
        :param: output
        :return: all_multiplicities, all_sz, ground_sz
        """
        read_multip = []
        for i, state in enumerate(qchem_file['excited_states']):
            read_multip.append(np.round(state['multiplicity'], 2))

        all_multiplicities = []
        for i in range(0, len(read_multip)):
            # Making multiplicity numbers multiples of doublets (s2=0.75).
            element = read_multip[i]
            n_times = np.round(element / 0.75)
            new_multip = 0.75 * n_times
            all_multiplicities.append(new_multip)

        # Obtain maximum and ground state multiplicity and s values
        s2_max = max(all_multiplicities)
        s2_ground = all_multiplicities[0]
        s_max = 0.5 * (-1 + np.sqrt(1 + 4 * s2_max))
        s_ground = 0.5 * (-1 + np.sqrt(1 + 4 * s2_ground))

        # Making multiplicity list from -n to +n in 1/2 intervals
        s_max = 0.5 * np.round(s_max / 0.5)
        s_ground = 0.5 * np.round(s_ground / 0.5)
        all_sz = list(np.arange(-s_max, s_max + 1, 1))
        ground_sz = list(np.arange(-s_ground, s_ground + 1, 1))

        return all_multiplicities, all_sz, ground_sz

    def get_all_socs(line, n_states, multiplicities, all_sz, soc_selection):
        """
        Get SOC matrix. For all the states it put values from maximum Sz to -Sz.
        If Sz does not exist (i.e., we consider Sz=-1.5 and Sz of the state is 0.5),
        then the SOC value is 0.
        :param: data, state_multiplicities, sz_list
        :return: soc_matrix
        """
        all_soc = np.zeros((n_states * len(all_sz), n_states * len(all_sz)), dtype=complex)
        # print('Multip:', state_multiplicities)
        # print('Sz:', sz_list)
        # print(len(soc_matrix[0,:]), len(soc_matrix[:,0]))
        # exit()

        for i in range(0, n_states):
            for j in range(0, n_states):

                if i != j:

                    i_multip = -0.5 * (-1 + np.sqrt(1 + 4 * multiplicities[i]))
                    i_multip = 0.5 * np.round(i_multip / 0.5)
                    # Making multiplicity numbers multiples of doublets (s=0.5)
                    j_multip = -0.5 * (-1 + np.sqrt(1 + 4 * multiplicities[j]))
                    j_multip = 0.5 * np.round(j_multip / 0.5)
                    # Making multiplicity numbers multiples of doublets (s=0.5)

                    i_index = all_sz.index(i_multip)
                    j_index = all_sz.index(j_multip)

                    i_position = i_index + len(all_sz) * i
                    j_position = j_index + len(all_sz) * j

                    i_j_soc_matrix = line[(i + 1, j + 1)][soc_selection]

                    for sz_1 in range(0, len(i_j_soc_matrix)):
                        for sz_2 in range(0, len(i_j_soc_matrix[0])):
                            all_soc[j_position + sz_1, i_position + sz_2] = i_j_soc_matrix[sz_1][sz_2]
                            # print('j i positions:', j_position+sz_1, i_position+sz_2, 'sz1 2:', sz_1, sz_2)

        # print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) \
        #                  for row in np.round((soc_matrix[:, :]), 5)]))
        # print(" ")
        # exit()

        # print('---------------')
        # state_A = 10
        # state_B = 8
        # for i in range(0, len(sz_list)):
        #     element_A = (state_A-1)*3+i
        #     element_B = (state_B-1)*3
        #     print(soc_matrix[element_A, element_B], soc_matrix[element_A, element_B+1], \
        #           soc_matrix[element_A, element_B+2])
        # exit()

        return all_soc

    def get_selected_states_socs(n_states, all_sz, socs):
        """
        Get SOC matrix between selected states. For all the states it put values
        from maximum Sz to -Sz. If Sz does not exist (i.e., we consider Sz=-1.5 and
        Sz of the state is 0.5), then the SOC value is 0.
        :param: selected_states, sz_list, all_soc
        :return: soc_matrix
        """
        selected_soc = np.zeros((len(n_states) * len(all_sz), len(n_states) * len(all_sz)), dtype=complex)

        for i, all_i in enumerate(n_states):
            for j, all_j in enumerate(n_states):

                for sz_1 in range(0, len(all_sz)):
                    for sz_2 in range(0, len(all_sz)):
                        i_index = i * len(all_sz) + sz_1
                        j_index = j * len(all_sz) + sz_2
                        all_i_index = (all_i - 1) * len(all_sz) + sz_1
                        all_j_index = (all_j - 1) * len(all_sz) + sz_2

                        selected_soc[i_index][j_index] = socs[all_i_index][all_j_index]

        # print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) \
        #                  for row in np.round((soc[:, :]), 5)]))
        # print(" ")
        # print(len(soc[0,:]), len(soc[:,0]))
        # exit()
        return selected_soc

    def get_doublets_soc(n_states, all_soc):
        """
        Get SOC matrix between selected states in doublets,
        meaning Sz = -0.5, 0.5.
        :param: selected_states, sz_list, all_soc
        :return: soc_matrix
        """
        doublets_soc = np.zeros((len(n_states) * 2, len(n_states) * 2), dtype=complex)

        for i, selected_i in enumerate(n_states):
            for j, selected_j in enumerate(n_states):

                for sz_1 in range(0, 2):
                    for sz_2 in range(0, 2):
                        i_index = i * 2 + sz_1
                        j_index = j * 2 + sz_2
                        all_i_index = i * len(sz_list) + (len(sz_list) // 2 - 1) + sz_1
                        all_j_index = j * len(sz_list) + (len(sz_list) // 2 - 1) + sz_2

                        doublets_soc[i_index][j_index] = all_soc[all_i_index][all_j_index]
        # print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) \
        #                  for row in np.round((doublet_soc[:, :]), 5)]))
        # print(" ---")
        # exit()
        return doublets_soc

    soc_search = 'total_soc_mat'
    state_multiplicities, sz_list, sz_ground = get_states_sz(output)
    all_socs = get_all_socs(data, totalstates, state_multiplicities, sz_list, soc_search)
    selected_socs = get_selected_states_socs(selected_states, sz_list, all_socs)
    doublet_soc = get_doublets_soc(selected_states, selected_socs)

    sz_list = [-0.5, 0.5]
    doublet_soc = doublet_soc / 219474.63068  # From cm-1 to a.u.
    return doublet_soc, sz_list, sz_ground


def get_spin_matrices(file, n_states):
    """
    Obtain the 3 dimensions of spin (s_x, s_y, s_z) from s2 of each state and put them in a 3-dim array.
    Spin is written with 'bra' in rows and 'ket' in columns, with spin order -Ms , +Ms.
    :param: file, n_states
    :return: spin_matrix: matrix with spins < m' | S | m >  in 3-D
    """

    def s2_from_file(file_qchem):
        """
        get s2 of each state from Q-Chem otuput
        :param: file
        :return: s2
        """
        search = ['  <S^2>      : ']
        elements = []

        with open(file_qchem, encoding="utf8") as file_qchem:
            for line in file_qchem:
                if any(ii in line for ii in search):
                    element = line[16:]
                    elements.append(element.split())

        s2_each_states = np.array(elements, dtype=float)
        return s2_each_states

    def s2_single_values(states_s2):
        """
        get s2 values from the s2 of all the states, to obtain all values of
        s2 only one time.
        :param: s2_states
        :return: s2_values
        """
        s2_values_list = []
        for ii in range(0, len(states_s2)):
            if states_s2[ii] not in s2_values_list:
                s2_values_list.append(states_s2[ii])

        s2_values = np.array(s2_values_list, dtype=float)
        return s2_values

    def s2_to_s(s2):
        """
        get total spin (s) from s^2
        :param: s2
        :return: total spin (s)
        """
        return 0.5 * (-1 + np.sqrt(1 + 4 * s2))

    def s_to_s2(spin):
        """
        get s2 from total spin (s)
        :param: spin: total spin (s)
        :return: s2
        """
        return spin * (spin + 1)

    def long_spin_matrices(s_x, s_y, s_z, multip_max, multip_state):
        """
        from sx, sy, sz with dimension of the state multiplicity to
        sx, sy, sz with dimension of the maximum multiplicity
        :param: s_x, s_y, s_z, multip_max, multip_state
        :return: s_x, s_y, s_z
        """
        long_sx = np.zeros((max_multiplicity, max_multiplicity), dtype=complex)
        long_sy = np.zeros((max_multiplicity, max_multiplicity), dtype=complex)
        long_sz = np.zeros((max_multiplicity, max_multiplicity), dtype=complex)

        multipl_difference = int((multip_max - multip_state) // 2)

        if multipl_difference != 0:
            for n_row in range(0, len(sx)):
                for n_column in range(0, len(sx)):
                    ii = n_row + multipl_difference
                    jj = n_column + multipl_difference
                    long_sx[ii, jj] = s_x[n_row, n_column]
                    long_sy[ii, jj] = s_y[n_row, n_column]
                    long_sz[ii, jj] = s_z[n_row, n_column]
        else:
            long_sx = s_x
            long_sy = s_y
            long_sz = s_z

        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
        #                  for row in np.round((long_sx[:, :]), 5)]))
        # exit()
        return long_sx, long_sy, long_sz

    # Abel function to determine the spin matrix:
    # https://github.com/abelcarreras/PyQchem/blob/1b1a0291f2737474955a5045ffbc56a2efd50911/pyqchem/tools/spin.py#L14

    def spin_matrices(spin):
        """
        Get spin matrices s_x, s_y, s_z between two spin states (s,m) and (s,m') such that
        sx = < m' | s_x | m >, sy = < m' | s_y | m > and sz = < m' | s_z | m >
        :param: spin: total spin (s)
        :return: s_x, s_y, s_z
        """

        def are_equal(a, b, thresh=1e-4):
            return abs(a - b) < thresh

        def sz_values(spinn):
            return np.arange(-spinn, spinn + 1)

        # spin-multiplicities
        multiplicity = len(sz_values(spin))

        # initialize s_x, s_y, s_z
        s_x = np.zeros((multiplicity, multiplicity), dtype=complex)
        s_y = np.zeros((multiplicity, multiplicity), dtype=complex)
        s_z = np.zeros((multiplicity, multiplicity), dtype=complex)

        # build spin matrices
        for ii, sz_bra in enumerate(sz_values(spin)):
            for g, sz_ket in enumerate(sz_values(spin)):

                if are_equal(sz_bra, sz_ket):
                    s_z[ii, g] = sz_ket

                if are_equal(sz_bra, sz_ket + 1):
                    s_x[ii, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[ii, g] = -0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)

                if are_equal(sz_bra, sz_ket - 1):
                    s_x[ii, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[ii, g] = 0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)

        return s_x, s_y, s_z

    s2_states = s2_from_file(file)  # s2 of each of the states
    single_s2_values = s2_single_values(s2_states)  # s2 of each of the states, without repetition
    number_of_spins = len(single_s2_values)

    # Spin matrix of non-relativistic states:
    max_multiplicity = int(2 * s2_to_s(max(single_s2_values)) + 1)
    ground_multiplicity = int(2 * s2_to_s(single_s2_values[0]) + 1)

    standard_spin_matrix = np.zeros((ground_multiplicity, ground_multiplicity, 3), dtype=complex)
    spin_matrix = np.zeros((len(n_states) * max_multiplicity, len(n_states) * max_multiplicity, 3), dtype=complex)

    for i in range(0, number_of_spins):
        for j in n_states:
            if s2_states[j - 1] == single_s2_values[i]:
                s = s2_to_s(single_s2_values[i])
                sx, sy, sz = spin_matrices(s)

                s_dim = len(sx) // 2 - 1
                total_dim = n_states.index(j) * 2

                spin_matrix[total_dim, total_dim, 0] = sx[s_dim, s_dim]
                spin_matrix[total_dim + 1, total_dim, 0] = sx[s_dim + 1, s_dim]
                spin_matrix[total_dim, total_dim + 1, 0] = sx[s_dim, s_dim + 1]
                spin_matrix[total_dim + 1, total_dim + 1, 0] = sx[s_dim + 1, s_dim + 1]

                spin_matrix[total_dim, total_dim, 1] = sy[s_dim, s_dim]
                spin_matrix[total_dim + 1, total_dim, 1] = sy[s_dim + 1, s_dim]
                spin_matrix[total_dim, total_dim + 1, 1] = sy[s_dim, s_dim + 1]
                spin_matrix[total_dim + 1, total_dim + 1, 1] = sy[s_dim + 1, s_dim + 1]

                spin_matrix[total_dim, total_dim, 2] = sz[s_dim, s_dim]
                spin_matrix[total_dim + 1, total_dim, 2] = sz[s_dim + 1, s_dim]
                spin_matrix[total_dim, total_dim + 1, 2] = sz[s_dim, s_dim + 1]
                spin_matrix[total_dim + 1, total_dim + 1, 2] = sz[s_dim + 1, s_dim + 1]

        # Standard spin matrix
        multip_difference = (max_multiplicity - ground_multiplicity) // 2
        for k in range(0, 3):
            for i in range(0, ground_multiplicity):
                for j in range(0, ground_multiplicity):
                    standard_spin_matrix[i, j, k] = spin_matrix[i + multip_difference, j + multip_difference, k]

    # print('Spin Matrices:')
    # for k in range(0,3):
    #    print('Dimension: ', k)
    #    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                     for row in np.round((spin_matrix[:,:,k]),5)]))
    #    print(" ")
    # exit()
    return spin_matrix, standard_spin_matrix


def get_orbital_matrices(file, totalstates, selected_states, sz_list):
    """
    Orbital angular momentum values are written in matrix with 'bra' in rows and 'ket' in columns,
    with spin order -Ms , +Ms. Third dimension is the direction.
    :param: file, totalstates, selected_states, sz_list
    :return: orbital_matrix
    """
    def get_all_momentum(line, n_states):
        """
        Get Lk between all the states in x,y,z dimensions. | A > in columns, < B | in rows.
        :param: line, n_states
        :return: all_orbital_momentum
        """
        all_orbital_momentum = np.zeros((n_states, n_states, 3), dtype=complex)

        for i in range(0, n_states):
            for j in range(0, n_states):
                if i != j:
                    element = line[(i + 1, j + 1)]['angular_momentum']

                    for k in range(0, 3):
                        all_orbital_momentum[j, i, k] = element[k]
        return all_orbital_momentum

    def get_selected_states_momentum(n_states, all_momentum):
        """
        Get Lk between the selected states in x,y,z dimensions.
        :param: n_states, all_momentum
        :return: selected_momentum
        """
        selected_momentum = np.zeros((len(n_states), len(n_states), 3), dtype=complex)

        for k in range(0, 3):
            for i, all_i in enumerate(n_states):
                for j, all_j in enumerate(n_states):
                    # print('i, j; ', i, j, 'all_i, all_j:', all_i-1, all_j-1)
                    selected_momentum[i][j][k] = all_momentum[all_i-1][all_j-1][k]
        return selected_momentum

    def get_doublets_momentum(n_states, selected_momentum):
        """
        Get Lk between the selected states in x,y,z dimensions for doublets,
        i.e. in each row (column) there are state A,-1/2 and A,+1/2.
        :param: n_states, selected_momentum
        :return: doublets_momentum
        """
        doublets_momentum = np.zeros((len(selected_momentum) * 2, len(n_states) * 2, 3), dtype=complex)

        for k in range(0, 3):
            for i in range(0, len(selected_momentum) * 2, 2):
                for j in range(0, len(selected_momentum) * 2, 2):

                    doublets_momentum[i, j][k] = selected_momentum[i // 2][j // 2][k]
                    doublets_momentum[i + 1, j + 1][k] = selected_momentum[i // 2][j // 2][k]
        return doublets_momentum

    def get_all_multip_momentum(all_momentums, all_sz):
        """
        Get Lk between the selected states in x,y,z dimensions for doublets,
        i.e. in each row (column) there are state A,-1/2 and A,+1/2.
        :param: all_momentums, all_sz
        :return: lk_values
        """
        lk_values = np.zeros((len(all_momentums) * len(all_sz), len(selected_states) * len(all_sz), 3), dtype=complex)

        for k in range(0, 3):
            for i in range(0, len(all_momentums) * len(all_sz), len(all_sz)):
                for j in range(0, len(all_momentums) * len(all_sz), len(all_sz)):

                    for multip in range(0, len(all_sz)):
                        lk_values[i + multip, j + multip][k] = all_momentums[i // len(all_sz)][j // len(all_sz)][k]
        return lk_values

    with open(file, encoding="utf8") as f:
        output = f.read()
    output = parser_rasci(output)
    data = output['interstate_properties']

    all_lk = get_all_momentum(data, totalstates)
    selected_lk = get_selected_states_momentum(selected_states, all_lk)
    all_multip_lk = get_all_multip_momentum(selected_lk, sz_list)

    doublets_lk = get_doublets_momentum(selected_states, selected_lk)
    all_multip_lk = doublets_lk
    return all_multip_lk


def get_hamiltonian_construction(selected_states, eigenenergies, spin_orbit_coupling, sz_values):
    """
    Hamiltonian is written 'bra' in rows and 'ket' in columns, with spin order -1/2 , +1/2.
    :param: selected_states, eigenenergies, spin_orbit_coupling
    :return: hamiltonian
    """
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


def diagonalization(hamiltonian):
    """
    1) Hamiltonian is diagonalized
    2) eigenvectors-eigenvalues are ordered by weight coefficients
    3) Doublet with the lowest energy is set as the new basis (Kramer doublet)
    :param: Hamiltonian
    :return: eigenvalues, eigenvectors, kramer_st: kramer_st is the index
    of the state set as Kramer doublet * 2
    """
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


def gfactor_calculation(lambda_matrix, sigma_matrix):
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
    return g_tensor_values


def from_energies_soc_to_g_values(file, states_ras, totalstates,
                                         excitation_energies_ras, soc_ras, sz_list):
    """"
    Obtention of the g-values from the eigenenergies and the SOCs.
    :param: file, states_ras, totalstates, excitation_energies_ras, soc_ras, sz_list
    :return: upper_g_matrix, g_values
    """
    hamiltonian_ras = get_hamiltonian_construction(states_ras, excitation_energies_ras, soc_ras, sz_list)

    eigenvalues, eigenvector, kramers_states = diagonalization(hamiltonian_ras)

    spin_matrix, standard_spin_matrix = get_spin_matrices(file, states_ras)

    l_matrix = get_orbital_matrices(file, totalstates, states_ras, sz_list)

    sigma_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, spin_matrix)

    lambda_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, l_matrix)

    g_values = gfactor_calculation(lambda_matrix, sigma_matrix)

    return g_values



