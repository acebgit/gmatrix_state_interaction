"""
Calculation of the g-tensor using Q-Chem output with RAS-CI
"""
import sys
import numpy as np
from numpy import linalg, sqrt
from pyqchem.parsers.parser_rasci import parser_rasci


def get_number_of_states(file):
    """
    Obtain the total number of states selected in ras
    :param: file_ms_notnull
    :return: selected_states
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
    return totalstates


def get_symmetry_states(file, totalstates):
    """
    Create two lists: first with the symmetry of each state (A1,A2,A3...) and second with
    the order of these symmetries (1A1,2A2,3A3...)
    :param: file_ras, selected_states
    :return: all_state_symmetries, ordered_state_symmetries
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    word_search = ['State symmetry: ']
    elements = []
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[2]
            element = element.replace('*', '')  # '*' means mix of symmetries
            elements.append(element)

        if len(elements) == totalstates:
            break  # All symmetries appended
    all_state_symmetries = np.array(elements)

    ordered_state_symmetries = []  # State symmetries with the order (1,2,3..)
    for i in range(0, len(all_state_symmetries)):
        number = 1
        for j in range(0, i):
            if all_state_symmetries[i] == all_state_symmetries[j]:
                number += 1
        element = str(number) + all_state_symmetries[i]
        ordered_state_symmetries.append(element)
    return all_state_symmetries, ordered_state_symmetries


def get_selected_states(file, totalstates, selected_states, states_option, symmetry_selection):
    """
    Select the states selected used depending on "states_option" value:
    0: use "state_ras" ; 1: use all states selected ; 2: use states selected by selected symmetry
    :param: file_ms_notnull, selected_states, selected_states, states_option, symmetry_selection
    :return: selected_states
    """
    all_symmetries, ordered_symmetries = get_symmetry_states(file, totalstates)

    if states_option == 0:  # Se
        for i in selected_states:
            if i <= 0 or i > totalstates:
                raise ValueError("The number of states selected selected must be among the total number of states selected calculated in QChem. "
                                 "Select a different number of states selected")

    elif states_option == 1:
        selected_states = list(range(1, totalstates + 1))

    elif states_option == 2:
        states_selected_by_symmetry = [1]  # GS is always added

        for nstate in range(1, totalstates + 1):
            if (all_symmetries[nstate - 1] == symmetry_selection) and (nstate != 1):
                states_selected_by_symmetry.append(nstate)
        selected_states = states_selected_by_symmetry

        if states_selected_by_symmetry == [1]:
            raise ValueError("There is not this symmetry. Change the symmetry selection.")
    return selected_states


def order_by_states(input_list, selected_states):
    """
    Order the list by the selected.
    :return:
    """
    output_list = []
    for i in selected_states:
        output_list.append(input_list[i - 1])
    return output_list


def get_eigenenergies(file, totalstates, selected_states):
    """
    Get energies of the RAS-CI states selected.
    :param: file_ms_notnull, selected_states, selected_states
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

    energies_selected = order_by_states(element_energy, selected_states)
    excited_energies_selected = order_by_states(elements_excitenergy, selected_states)

    eigenenergies = np.array(energies_selected, dtype=float)
    excitation_energies_ev = np.array(excited_energies_selected, dtype=float)
    excitation_energies = excitation_energies_ev / 27.211399  # From eV to a.u.
    return eigenenergies, excitation_energies


def get_spin_orbit_couplings(file, totalstates, selected_states, soc_option):
    """
    Spin-orbit coupling values are written in matrix with 'bra' in rows
    and 'ket' in columns, with spin order -Ms , +Ms from first to last selected states selected.
    :param: file_ms_notnull, selected_states, selected_states, soc_option
    :return: doublet_soc, sz_list
    """
    def get_states_sz(qchem_file, states_selected):
        """
        Get S² and Sz of all states selected
        :param: output
        :return: all_multip, all_sz, ground_sz
        """
        def from_s2_to_sz_list(s2):
            # Obtain s values
            s = 0.5 * (-1 + np.sqrt(1 + 4 * s2))
            # Making multiplicity list from -n to +n in 1/2 intervals
            s = 0.5 * np.round(s / 0.5)
            sz = list(np.arange(-s, s + 1, 1))
            return sz

        read_multip = []
        for i, state in enumerate(qchem_file['excited_states']):
            read_multip.append(np.round(state['multiplicity'], 2))

        all_multip = []
        for i in range(0, len(read_multip)):
            # Making multiplicity numbers multiples of doublets (s2=0.75).
            element = read_multip[i]
            n_times = np.round(element / 0.75)
            new_multip = 0.75 * n_times
            all_multip.append(new_multip)

        ordered_multip = []
        for i in range(0, len(states_selected)):
            selected_mult = states_selected[i] - 1
            ordered_multip.append(all_multip[selected_mult])

        # print(all_multip)
        # print(ordered_multip)
        # exit()

        all_sz = from_s2_to_sz_list(max(all_multip))
        ground_sz = from_s2_to_sz_list(ordered_multip[0])
        if len(ground_sz) == 1:
            raise ValueError("Warning! It is not allowed the calculation of the g-tensor in a singlet ground state. "
                             "Ground state corresponds to the first of the included states selected.")
        return all_multip, all_sz, ground_sz

    def get_all_socs(line, n_states, multiplicities, all_sz, soc_selection):
        """
        Get SOC matrix. For all the states selected it put values from maximum Sz to -Sz.
        If Sz does not exist (i.e., we consider Sz=-1.5 and Sz of the state is 0.5),
        then the SOC value is 0.
        :param: data, all_multip, sz_list
        :return: soc_matrix
        """
        all_soc = np.zeros((n_states * len(all_sz), n_states * len(all_sz)), dtype=complex)
        # print('Multip:', all_multip)
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
        Get SOC matrix between selected states selected. For all the states selected it put values
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

    with open(file, encoding="utf8") as f:
        output = f.read()
    output = parser_rasci(output)
    data = output['interstate_properties']

    soc_search = 'total_soc_mat'
    if soc_option == 0:
        pass
    elif soc_option == 1:
        soc_search = '1e_soc_mat'
    elif soc_option == 2:
        soc_search = '2e_soc_mat'

    all_multiplicities, maximum_sz_all_states, sz_ground_state = get_states_sz(output, selected_states)
    all_socs = get_all_socs(data, totalstates, all_multiplicities, maximum_sz_all_states, soc_search)
    selected_socs = get_selected_states_socs(selected_states, maximum_sz_all_states, all_socs)

    # print('SOC:')
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                 for row in np.round((selected_socs[:,:]),5)]))
    # exit()

    selected_socs = selected_socs / 219474.63068  # From cm-1 to a.u.
    return selected_socs, maximum_sz_all_states, sz_ground_state


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
                print('positions: ', i // 2, 'value:', matrix[i, j])
                print('positions: ', j // 2, 'value:', matrix[j, i])
                raise ValueError("Matrix is not Hermitian: see the elements shown above")


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
    # print('Hamiltonian:')
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                 for row in np.round((hamiltonian[:,:]),5)]))
    # exit()
    return hamiltonian


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


def diagonalization(initial_matrix):
    """
    1) Hamiltonian is diagonalized
    2) eigenvectors-eigenvalues are ordered by weight coefficients
    3) Doublet with the lowest energy is set as the new basis (Kramer doublet)
    :param: initial_matrix
    :return: eigenvalues, eigenvectors, diagonal_matrix
    """
    eigenvalues, eigenvectors = linalg.eigh(initial_matrix)
    eigenvalues, eigenvectors = reordering_eigenvectors(eigenvalues, eigenvectors)

    rotation_inverse = np.linalg.inv(eigenvectors)
    diagonal_matrix = np.matmul(np.matmul(rotation_inverse, initial_matrix), eigenvectors)

    # for i in range(0, len(eigenvalues)):
    #     print('eigenvalue:', eigenvalues[i])
    #     print('eigenvector:', eigenvectors[:, i])
    #     print()
    # exit()
    return eigenvalues, eigenvectors, diagonal_matrix


def get_spin_matrices(file, n_states):
    """
    Obtain the 3 dimensions of spin (s_x, s_y, s_z) from s2_all of each state and put them in a 3-dim array.
    Spin is written with 'bra' in rows and 'ket' in columns, with spin order -Ms , +Ms.
    Abel functions to determine the spin matrix: https://github.com/abelcarreras/PyQchem/blob/1b1a0291f2737474955a5045ffbc56a2efd50911/pyqchem/tools/spin.py#L14
    :param: file_ms_notnull, n_states
    :return: spin_matrix: matrix with spins < m' | S | m >  in 3-D
    """

    def get_all_s2(file_qchem, n_states):
        """
        get s2_all of each state from Q-Chem otuput
        :param: file_ms_notnull
        :return: s2_all
        """
        search = ['  <S^2>      : ']
        s2_all_list = []

        with open(file_qchem, encoding="utf8") as file_qchem:
            for line in file_qchem:
                if any(a in line for a in search):
                    line = line.split()
                    element = line[2]
                    s2_all_list.append(element)

        s2_nstates = np.zeros(len(n_states), dtype=float)
        for i in range(0, len(n_states)):
            selected_mult = n_states[i] - 1
            s2_nstates[i] = s2_all_list[selected_mult]

        if s2_nstates[0] == 0:
            raise ValueError("Warning! It is not allowed the calculation of the g-tensor in a singlet ground state. "
                             "Ground state corresponds to the first of the included states selected.")

        s2_all = np.array(s2_all_list, dtype=float)
        return s2_all, s2_nstates

    def get_s2_single_values(list_s2):
        """
        get s2_all values from the s2_all of all the states selected, to obtain all values of
        s2_all only one time.
        :param: s2_all_states
        :return: s2_values
        """
        set_s2 = set(list_s2)
        list_s2 = list(set_s2)
        # matrix_s2 = np.array(s2_all, dtype=float)
        return list_s2

    def s2_to_s(s2):
        """
        get total spin (s) from s^2
        :param: s2_all
        :return: total spin (s)
        """
        return 0.5 * (-1 + np.sqrt(1 + 4 * s2))

    def s_to_s2(spin):
        """
        get s2_all from total spin (s)
        :param: spin: total spin (s)
        :return: s2_all
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
                    iii = n_row + multipl_difference
                    jj = n_column + multipl_difference
                    long_sx[iii, jj] = s_x[n_row, n_column]
                    long_sy[iii, jj] = s_y[n_row, n_column]
                    long_sz[iii, jj] = s_z[n_row, n_column]
        else:
            long_sx = s_x
            long_sy = s_y
            long_sz = s_z

        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
        #                  for row in np.round((long_sx[:, :]), 5)]))
        # exit()
        return long_sx, long_sy, long_sz

    def spin_matrices(spin):
        """
        Get spin matrices s_x, s_y, s_z between two spin states selected (s,m) and (s,m') such that
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
        for iii, sz_bra in enumerate(sz_values(spin)):
            for g, sz_ket in enumerate(sz_values(spin)):

                if are_equal(sz_bra, sz_ket):
                    s_z[iii, g] = sz_ket

                if are_equal(sz_bra, sz_ket + 1):
                    s_x[iii, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[iii, g] = -0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)

                if are_equal(sz_bra, sz_ket - 1):
                    s_x[iii, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[iii, g] = 0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)

        return s_x, s_y, s_z

    s2_all_states, s2_selected_states = get_all_s2(file, n_states)  # s2_all of each of the states selected
    single_s2_values = get_s2_single_values(s2_all_states)  # s2_all of each of the states selected, without repetition

    max_multiplicity = int(2 * s2_to_s(max(s2_all_states)) + 1)
    ground_multiplicity = int(2 * s2_to_s(s2_selected_states[0]) + 1)

    standard_spin_matrix = np.zeros((ground_multiplicity, ground_multiplicity, 3), dtype=complex)
    spin_matrix = np.zeros((len(n_states) * max_multiplicity, len(n_states) * max_multiplicity, 3), dtype=complex)

    for i in range(0, len(single_s2_values)):
        for j in n_states:
            if s2_all_states[j - 1] == single_s2_values[i]:
                s = s2_to_s(single_s2_values[i])
                sx, sy, sz = spin_matrices(s)
                sx, sy, sz = long_spin_matrices(sx, sy, sz, max_multiplicity, s2_all_states[j - 1])

                s_dim = 0
                total_dim = n_states.index(j) * max_multiplicity

                for row in range(0, max_multiplicity):
                    for column in range(0, max_multiplicity):
                        spin_matrix[total_dim, total_dim, 0] = sx[s_dim, s_dim]
                        spin_matrix[total_dim, total_dim + column, 0] = sx[s_dim, s_dim + column]
                        spin_matrix[total_dim + row, total_dim, 0] = sx[s_dim + row, s_dim]
                        spin_matrix[total_dim + row, total_dim + column, 0] = sx[s_dim + row, s_dim + column]

                        spin_matrix[total_dim, total_dim, 1] = sy[s_dim, s_dim]
                        spin_matrix[total_dim, total_dim + column, 1] = sy[s_dim, s_dim + column]
                        spin_matrix[total_dim + row, total_dim, 1] = sy[s_dim + row, s_dim]
                        spin_matrix[total_dim + row, total_dim + column, 1] = sy[s_dim + row, s_dim + column]

                        spin_matrix[total_dim, total_dim, 2] = sz[s_dim, s_dim]
                        spin_matrix[total_dim, total_dim + column, 2] = sz[s_dim, s_dim + column]
                        spin_matrix[total_dim + row, total_dim, 2] = sz[s_dim + row, s_dim]
                        spin_matrix[total_dim + row, total_dim + column, 2] = sz[s_dim + row, s_dim + column]

        # Standard spin matrix
        multip_difference = (max_multiplicity - ground_multiplicity) // 2
        for k in range(0, 3):
            for ii in range(0, ground_multiplicity):
                for j in range(0, ground_multiplicity):
                    standard_spin_matrix[ii, j, k] = spin_matrix[ii + multip_difference, j + multip_difference, k]

    # print('Spin Matrices:')
    # for k in range(0,3):
    #    print('Dimension: ', k)
    #    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                     for row in np.round((standard_spin_matrix[:,:,k]),5)]))
    #    print(" ")
    # exit()
    return spin_matrix, standard_spin_matrix


def get_orbital_matrices(file, totalstates, selected_states, sz_list):
    """
    Orbital angular momentum values are written in matrix with 'bra' in rows and 'ket' in columns,
    with spin order -Ms , +Ms. Third dimension is the direction.
    :param: file_ms_notnull, selected_states, selected_states, sz_list
    :return: orbital_matrix
    """
    def get_all_momentum(line, n_states):
        """
        Get Lk between all the states selected in x,y,z dimensions. | A > in columns, < B | in rows.
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
        Get Lk between the selected states selected in x,y,z dimensions.
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

    def get_all_multip_momentum(all_momentums, all_sz):
        """
        Get Lk between the selected states selected in x,y,z dimensions for doublets,
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

    # print('Orbital Matrices:')
    # for k in range(0,3):
    #    print('Dimension: ', k)
    #    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                     for row in np.round((all_multip_lk[:,:,k]),5)]))
    #    print(" ")
    # exit()
    return all_multip_lk


def angular_matrixes_obtention(eigenvectors, input_angular_matrix, sz_list):
    """
    Spin or orbital angular matrix calculation. The matrix has the dimension
    len(sz_values)
    :param: sz_values, eigenvectors, input_angular_matrix
    :return: angular_matrix: contains the spin value < B(S,Sz) | Sx | A(S',Sz') >,
    meaning < B(S,Sz)| corresponds to state "i" (in rows) and | A(S',Sz') > to
    state "j" (in columns)
    """
    angular_matrix = np.zeros((len(sz_list), len(sz_list), 3), dtype=complex)

    for k in range(0, 3):
        for row in range(0, len(sz_list)):  # state i
            for column in range(0, len(sz_list)):  # state j

                for bra in range(0, len(eigenvectors)):  # state <B|
                    for ket in range(0, len(eigenvectors)):  # state |A>

                        coeff_bra = np.conj(eigenvectors[bra, row])
                        coeff_ket = (eigenvectors[ket, column])
                        angular_value = input_angular_matrix[bra, ket, k]

                        element = coeff_bra * coeff_ket * angular_value
                        angular_matrix[row, column, k] += element

    # print('Angular matrix with all spin angular momentums:')
    # for k in range(0, 3):
    #     print('Dimension: ', k)
    #     print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
    #                      for row in np.round((angular_matrix[:, :, k]), 5)]))
    #     print(" ")
    # exit()
    return angular_matrix


def g_factor_calculation(standard_spin_matrix, s_matrix, l_matrix, sz_list, ground_sz):
    """
    Calculation of the g-shift with orbital and spin angular momentum matrices. J-values
    are reorganized using sz_list and ground_sz.
    :param: standard_spin_matrix, s_matrix, l_matrix, sz_list, ground_sz
    :return: g_shifts
    """
    lande_factor = 2.002319304363

    def j_diagonalization(initial_matrix, sz_list, ground_sz):
        """
        J-matrix diagonalization including the reorganization of this diagonal values.
        WARNING!: most of the problems occur when the final diagonal matrix do not have
        the values in the center of the diagonal:
        - if dimension 4, the two values in 1,1 and 2,2
        - if dimension 5, values in 2,2 and 3,3
        :param: initial_matrix
        :return: eigenvalues, eigenvectors, diagonal_matrix
        """
        eigenvalues, eigenvectors = linalg.eigh(initial_matrix)
        # reordering_eigenvectors(eigenvalues, eigenvectors)

        # for i in range(0, len(eigenvectors)):
        #     print(np.round(eigenvectors[i,:]))
        # print()
        #
        # lista = [0, 2, 1, 2, 2, 3, 3, 4]
        # nuevo = np.zeros(len(eigenvectors), dtype=complex)
        # for i in range(0, len(lista), 2):
        #     nuevo[:] = eigenvectors[lista[i+1],:]
        #     eigenvectors[lista[i+1],:] = eigenvectors[lista[i],:]
        #     eigenvectors[lista[i],:] = nuevo
        #
        # for i in range(0, len(eigenvectors)):
        #     print(np.round(eigenvectors[i,:]))
        # print()
        # exit()

        rotation_inverse = np.linalg.inv(eigenvectors)
        diagonal_matrix = np.matmul(np.matmul(rotation_inverse, initial_matrix), eigenvectors)
        # print('Rotation matrix:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
        #                  for row in np.round((diagonal_matrix[:, :]), 5)]))
        # print(" ")

        # Obtain the difference between the maximum sz, that determines J dimensions, and the ground
        # state sz to set the positions of j-values in the diagonal
        # sz_max = max(sz_list)
        # if sz_max - ground_sz != 0:
        #     s_difference = int(sz_max - ground_sz)
        #
        # if abs(diagonal_matrix[s_difference, s_difference]) <= 0.0001:
        #     len_diag = len(diagonal_matrix) - 1
        #
        #     position_1 = s_difference
        #     position_2 = len_diag - s_difference
        #     # Move the values to correct positions
        #     diagonal_matrix[position_1, position_1] = diagonal_matrix[0, 0]
        #     diagonal_matrix[position_2, position_2] = diagonal_matrix[len_diag, len_diag]
        #     # Set the previous j-values positions to zero
        #     diagonal_matrix[0, 0] = 0
        #     diagonal_matrix[len_diag, len_diag] = 0
        # print('Rotation matrix:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
        #                  for row in np.round((diagonal_matrix[:, :]), 5)]))
        # print("------- ")
        # exit()
        return eigenvalues, eigenvectors, diagonal_matrix

    def trace_g_values(total_j, spin_matr):
        a = np.matmul(total_j, spin_matr)
        b = np.matmul(spin_matr, spin_matr)
        g_value = np.trace(a) / np.trace(b)
        return g_value

    def j_matrix_formation(landefactor, spin, orbital, list_sz, sz_ground):
        j_big_matrix = landefactor * spin + orbital
        sz_difference = (len(list_sz) - len(sz_ground)) // 2
        j_matrix = np.zeros((len(sz_ground), len(sz_ground), 3), dtype=complex)
        for k in range(0, 3):
            for i in range(0, len(j_matrix)):
                for j in range(0, len(j_matrix)):
                    j_matrix[i, j, k] = j_big_matrix[i + sz_difference, j + sz_difference, k]
            hermitian_test(j_matrix[:, :, k])

        # print('J-matrix:')
        # for k in range(0,3):
        #    print('Dimension: ', k)
        #    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                     for row in np.round((j_matrix[:,:,k]),5)]))
        #    print(" ")
        # exit()
        return j_matrix

    j_matrix = j_matrix_formation(lande_factor, s_matrix, l_matrix, sz_list, ground_sz)

    # 1) g-value zz
    j_eigenvalues_z, j_matrix_rotation, j_matrix_diagonal_z = j_diagonalization(j_matrix[:, :, 2], sz_list, ground_sz)
    g_matrix_triangular = np.zeros((3, 3), dtype=complex)
    g_matrix_triangular[2, 2] = trace_g_values(j_matrix_diagonal_z, standard_spin_matrix[:, :, 2])

    # 2) pseudospin z
    pseudospin_matrix = np.zeros((len(j_matrix[0, :]), len(j_matrix[:, 0]), 3), dtype=complex)
    pseudospin_matrix[:, :, 2] = j_matrix[:, :, 2] / g_matrix_triangular[2, 2]

    # 3) g-value yz
    g_matrix_triangular[1, 2] = trace_g_values(j_matrix[:, :, 1], pseudospin_matrix[:, :, 2])

    # 4) Residue
    residue = j_matrix[:, :, 1] - g_matrix_triangular[1, 2] * pseudospin_matrix[:, :, 2]

    # 5) g-value yy
    j_eigenvalues_z, j_matrix_rotation, j_matrix_transformed_y = \
        j_diagonalization(j_matrix[:, :, 1], sz_list, ground_sz)
    g_matrix_triangular[1, 1] = trace_g_values(j_matrix_transformed_y, standard_spin_matrix[:, :, 2])

    # 6) pseudospin y
    pseudospin_matrix[:, :, 1] = residue[:, :] / g_matrix_triangular[1, 1]

    # 7) g-value xy and xz
    g_matrix_triangular[0, 1] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 1])
    g_matrix_triangular[0, 2] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 2])

    # 8) g-value xx
    j_eigenvalues_z, j_matrix_rotation, j_matrix_transformed_x = \
        j_diagonalization(j_matrix[:, :, 0], sz_list, ground_sz)
    g_matrix_triangular[0, 0] = trace_g_values(j_matrix_transformed_x, standard_spin_matrix[:, :, 2])
    # print('g-matrix:')
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                 for row in np.round((g_matrix_triangular[:,:]),5)]))
    # exit()

    # 9) g-shifts
    upper_g_matrix = np.matmul(g_matrix_triangular, np.transpose(g_matrix_triangular))
    upper_g_matrix_eigenvalues, rotation_g_matrix, upper_g_matrix_diagonal = diagonalization(upper_g_matrix)

    g_shifts = np.zeros(3, dtype=complex)
    for i in range(0, 3):
        g_shifts[i] = (sqrt(upper_g_matrix_diagonal[i, i]) - lande_factor) * 1000
    return g_shifts


def from_energies_soc_to_g_values(file, states_ras, totalstates,
                                  excitation_energies_ras, soc_ras, sz_list, ground_sz):
    """"
    Obtention of the g-values from the eigenenergies and the SOCs.
    :param:file_ms_notnull, states_ras, selected_states, excitation_energies_ras, soc_ras, sz_list, ground_sz
    :return: g_shift
    """
    hamiltonian_ras = get_hamiltonian_construction(states_ras, excitation_energies_ras, soc_ras, sz_list)

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian_ras)

    spin_matrix, standard_spin_matrix = get_spin_matrices(file, states_ras)

    orbital_matrix = get_orbital_matrices(file, totalstates, states_ras, sz_list)

    combination_spin_matrix = angular_matrixes_obtention(eigenvector, spin_matrix, sz_list)

    combination_orbital_matrix = angular_matrixes_obtention(eigenvector, orbital_matrix, sz_list)

    g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                                   sz_list, ground_sz)

    return g_shift


def print_g_calculation(file, totalstates, selected_states,
                        states_ras, upper_g_tensor_results_ras, symmetry_selection):
    print("--------------------------------------")
    print("     INPUT SECTION")
    print("--------------------------------------")
    print("File selected: ", file)
    print("Number of states selected: ", totalstates)
    if selected_states == 2:
        print("Symmetry: ", symmetry_selection)
        print("Selected states selected: ", states_ras)
    else:
        print("Selected states selected: ", states_ras)

    print(" ")
    print("------------------------")
    print(" RAS-CI RESULTS")
    print("------------------------")
    print('g-factor (x y z dimensions):')
    print(np.round(upper_g_tensor_results_ras.real[0], 3), np.round(upper_g_tensor_results_ras.real[1], 3),
          np.round(upper_g_tensor_results_ras.real[2], 3))
    print('')


def gfactor_presentation(ras_input, states_ras, selected_states, symmetry_selection, soc_options, ppm):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, states_ras, selected_states, symmetry_selection, soc_options
    :return: g-shifts
    """
    totalstates = get_number_of_states(ras_input)

    states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

    selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options)

    hamiltonian_ras = get_hamiltonian_construction(states_ras, excitation_energies_ras, selected_socs, sz_list)

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian_ras)

    spin_matrix, standard_spin_matrix = get_spin_matrices(ras_input, states_ras)

    orbital_matrix = get_orbital_matrices(ras_input, totalstates, states_ras, sz_list)

    combination_spin_matrix = angular_matrixes_obtention(eigenvector, spin_matrix, sz_list)

    combination_orbital_matrix = angular_matrixes_obtention(eigenvector, orbital_matrix, sz_list)

    g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                                   sz_list, sz_ground)

    if ppm == 1:
        g_shift = g_shift * 1000

    print_g_calculation(ras_input, totalstates, selected_states, states_ras, g_shift, symmetry_selection)


def from_gvalue_to_shift(lista):
    """
    Obtain the g-shifts from the g-values
    :param: list
    :return: g_shift
    """
    lande_factor = 2.002319304363
    g_shift = []
    for i in range(0, len(lista)):
        value = (lista[i] - lande_factor) * 10**6
        g_shift.append(value)
    print(np.round(g_shift, 3))


# def gfactor_two_files(file_ms_notnull, file_ms_null, states_ras, states_option):
#     """
#     Returns the g-shifts for doublet ground state molecules.
#     :param: file_ms_notnull, selected_states, selected_states, symmetry_selection, soc_options
#     :return: g-shifts
#     """
#     def get_energies_socs(file, selected_states, states_option):
#         """
#         Having the selected states selected, get the energy and SOCs between them
#         :param: file, selected_states, states_option
#         :return: totalstates, eigenenergies_ras, selected_socs, sz_list, sz_ground
#         """
#         totalstates = get_number_of_states(file)
#
#         symmetry_selections = 'None'
#         selected_states = get_selected_states(file, totalstates, selected_states, states_option, symmetry_selections)
#
#         eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, selected_states)
#         selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(file, totalstates, selected_states, soc_option=0)
#         return totalstates, eigenenergies_ras, selected_socs, sz_list, sz_ground
#
#     def mapping_between_states(ms_notnull_energies, ms_null_energies):
#         """
#         Comparing the energies, mapping between states selected with Ms = 0, that do not have coupling between states selected triplets,
#         and states selected with Ms not 0, that do have this coupling.
#         :return:
#         """
#         mapping_list = []
#
#         for i in range(0, len(ms_notnull_energies)):
#             for j in range(0, len(ms_null_energies)):
#                 ener_ms_notnull = np.round(ms_notnull_energies[i], 5)
#                 ener_ms_null = np.round(ms_null_energies[j], 5)
#
#                 if ener_ms_notnull == ener_ms_null:
#                     mapping_dict = {'state ms not null': i, 'state ms null': j}
#                     mapping_list.append(mapping_dict)
#
#         # for mapping_dict in mapping_list:
#         #     a = mapping_dict['state ms not null']
#         #     b = mapping_dict['state ms null']
#         #     print('file_ms_notnull', states_ras[a], 'file_ms_null', states_ras[b])
#         # exit()
#         return mapping_dict, mapping_list
#
#     def exchange_coupling(mapping_list, socs_msnull, socs_msnotnull, sz_list):
#         """
#         Put SOCs between states selected with Ms different than 0 (that are obtained in the output) in the
#         SOC matrix of states selected with Ms 0 (that are not obtained since Clebsh-Gordan coefficient is too small)
#         :return:
#         """
#         for i in mapping_list:
#             for j in mapping_list:
#                 if i['state ms not null'] != j['state ms not null']:  # If states selected are not the same (in Ms not null list)
#                     i_state_ms_notnull = i['state ms not null']
#                     j_state_ms_notnull = j['state ms not null']
#
#                     i_state_ms_null = i['state ms null']
#                     j_state_ms_null = j['state ms null']
#
#                     # print('State Ms not null', states_ras[i_state_ms_notnull], '(', i_state_ms_notnull, ')',
#                     #       states_ras[j_state_ms_notnull], '(', j_state_ms_notnull, ')', '; ',
#                     #       'State Ms null', states_ras[i_state_ms_null], '(', i_state_ms_null, ')',
#                     #       states_ras[j_state_ms_null], '(', j_state_ms_null, ')',)
#
#                     for sz_1 in range(0, len(sz_list)):
#                         for sz_2 in range(0, len(sz_list)):
#                             i_msnotnull = i_state_ms_notnull * len(sz_list) + sz_1
#                             j_msnotnull = j_state_ms_notnull * len(sz_list) + sz_2
#
#                             i_msnull = i_state_ms_null * len(sz_list) + sz_1
#                             j_msnull = j_state_ms_null * len(sz_list) + sz_2
#
#                             # print(socs_msnotnull[i_msnull, j_msnull], '<--->', socs_msnull[i_msnotnull, j_msnotnull])
#                             socs_msnotnull[i_msnull, j_msnull] = socs_msnull[i_msnotnull, j_msnotnull]
#                     # print('--END LOOP--')
#         # print('SOC:')
#         # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
#         #                 for row in np.round((socs_msnotnull[:,:]),5)])) # * 219474.63068
#         # exit()
#         return socs_msnotnull
#
#     # print('File ms not null: ', file_ms_notnull)
#     # print('File ms null: ', file_ms_null)
#     # print('States: ', states_ras)
#
#     totalstates_ms_notnull, energies_ms_notnull, socs_msnull, sz_list_ms_notnull, sz_ground_ms_notnull = get_energies_socs(file_ms_notnull, states_ras, states_option)
#
#     totalstates_ms_null, energies_ms_null, socs_msnotnull, sz_list_ms_null, sz_ground_ms_null = get_energies_socs(file_ms_null, states_ras, states_option)
#
#     mapping_dict, mapping_list = mapping_between_states(energies_ms_notnull, energies_ms_null)
#
#     socs_msnotnull = exchange_coupling(mapping_list, socs_msnull, socs_msnotnull, sz_list_ms_notnull)
#
#     hamiltonian_ras = get_hamiltonian_construction(states_ras, energies_ms_null, socs_msnotnull, sz_list_ms_null)
#
#     eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian_ras)
#
#     spin_matrix, standard_spin_matrix = get_spin_matrices(file_ms_notnull, states_ras)
#
#     orbital_matrix = get_orbital_matrices(file_ms_notnull, totalstates_ms_null, states_ras, sz_list_ms_null)
#
#     combination_spin_matrix = angular_matrixes_obtention(eigenvector, spin_matrix, sz_list_ms_null)
#
#     combination_orbital_matrix = angular_matrixes_obtention(eigenvector, orbital_matrix, sz_list_ms_null)
#
#     g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
#                                    sz_list_ms_null, sz_ground_ms_null)
#
#     print_g_calculation(file_ms_notnull, totalstates_ms_null, states_ras, states_ras, g_shift * 1000, symmetry_selection=0)
