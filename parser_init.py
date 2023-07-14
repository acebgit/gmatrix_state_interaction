import sys
import numpy as np
from pyqchem.parsers.parser_rasci import parser_rasci

def get_scf_energy(file):
    """
    Obtain SCF energy
    :param: ras_input
    :return: scf_energy
    """
    word_search = ['SCF   energy in the final basis set']

    with open(file, encoding="utf8") as data:
        for line in data:
            if any(i in line for i in word_search):
                element = line[39:]
                scf_energy = float(element)
    return scf_energy


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
    return totalstates


def get_selected_states(file, totalstates, selected_states, states_option, symmetry_selection):
    """
    Select the states used depending on "states_option" value:
    0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
    :param: file, totalstates, selected_states, states_option, symmetry_selection
    :return: selected_states
    """
    all_symmetries, ordered_symmetries = get_symmetry_states(file, totalstates)

    if states_option == 0:  # Se
        for i in selected_states:
            if i <= 0 or i > totalstates:
                print("The number of states selected must be among the total number of states calculated in QChem.")
                print("Select a different number of states")
                sys.exit()

    elif states_option == 1:
        selected_states = list(range(1, totalstates + 1))

    elif states_option == 2:
        states_selected_by_symmetry = [1]  # GS is always added

        for nstate in range(1, totalstates + 1):
            if (all_symmetries[nstate - 1] == symmetry_selection) and (nstate != 1):
                states_selected_by_symmetry.append(nstate)
        selected_states = states_selected_by_symmetry

        if states_selected_by_symmetry == [1]:
            print('There is not this symmetry.')
            print('Change the symmetry selection:')
            for nstate in range(0, totalstates):
                print('- State', nstate + 1, ', symmetry', all_symmetries[nstate])
            sys.exit()
    return selected_states


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


def get_symmetry_states(file, totalstates):
    """
    Create two lists: first with the symmetry of each state (A1,A2,A3...) and second with
    the order of these symmetries (1A1,2A2,3A3...)
    :param: file_ras, totalstates
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


def get_hole_part_contributions(file, totalstates):
    """
    Take the hole and particle contributions of each state.
    :param: file, totalstates
    :return: hole_contributions, part_contributions
    """
    with open(file, encoding="utf-8") as file:
        data = file.readlines()

    word_search = ['   Hole: ']
    elements = []

    for line in data:
        if any(i in line for i in word_search):
            element = line[39:]
            elements.append(element.split())
        if len(elements) == totalstates:
            break

    hole_contributions = np.array(elements, dtype=float)

    word_search = ['   Part: ']
    elements = []

    for line in data:
        if any(i in line for i in word_search):
            element = line[39:]
            elements.append(element.split())
        if len(elements) == totalstates:
            break
    part_contributions = np.array(elements, dtype=float)
    return hole_contributions, part_contributions


def get_mulliken_spin(file, totalstates, states):
    """
    Get Mulliken charge and spin of the first atom, i.e. a metal in transition atom complexes.
    :param: file, totalstates, states
    :return: charge_mulliken, spin_mulliken
    """
    word_search = '    Mulliken population analysis '
    element_charge = []
    elements_spin = []

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                for i in range(0, 4):  # Take the 5th line
                    next_line = next(file)
                    i += 1

                next_line = next_line.split()

                element = next_line[2]
                element_charge.append(element)

                element = next_line[3]
                elements_spin.append(element)
            if len(element_charge) == totalstates:
                break

    charge_mulliken_selected = []
    spin_mulliken_selected = []

    for i in states:
        charge_mulliken_selected.append(element_charge[i - 1])
        spin_mulliken_selected.append(elements_spin[i - 1])

    charge_mulliken = np.array(charge_mulliken_selected, dtype=float)
    spin_mulliken = np.array(spin_mulliken_selected, dtype=float)
    return charge_mulliken, spin_mulliken


def get_spin_orbit_couplings(file, totalstates, selected_states, soc_option, bolvin):
    """
    Spin-orbit coupling values are written in matrix with 'bra' in rows
    and 'ket' in columns, with spin order -Ms , +Ms.
    :param: file, totalstates, selected_states, soc_option
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
    if soc_option == 0:
        pass
    elif soc_option == 1:
        soc_search = '1e_soc_mat'
    elif soc_option == 2:
        soc_search = '2e_soc_mat'

    state_multiplicities, sz_list, sz_ground = get_states_sz(output)
    all_socs = get_all_socs(data, totalstates, state_multiplicities, sz_list, soc_search)
    selected_socs = get_selected_states_socs(selected_states, sz_list, all_socs)

    if bolvin == 1:
        doublet_soc = get_doublets_soc(selected_states, selected_socs)
        selected_socs = doublet_soc
        sz_list = [-0.5, 0.5]

    selected_socs = selected_socs / 219474.63068  # From cm-1 to a.u.
    return selected_socs, sz_list, sz_ground


def get_socc_values(file, totalstates):
    """
    Get spin-orbit coupling constant between states
    :param: file, totalstates
    :return: socc
    """
    with open(file, encoding="utf-8") as file:
        data = file.readlines()

    word_search = ['Mean-Field SOCC']
    elements = []
    n_states = 0

    elements.append('0.000000')  # SOCC between state 1 and 1 is zero
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[3]
            elements.append(element)
            n_states += 1
        if n_states == totalstates - 1:
            break

    socc = np.array(elements, dtype=float)
    return socc


def get_spin_matrices(file, n_states, bolvin):
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
                if any(k in line for k in search):
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
        for k in range(0, len(states_s2)):
            if states_s2[k] not in s2_values_list:
                s2_values_list.append(states_s2[k])

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

        multip_difference = int((multip_max - multip_state) // 2)

        if multip_difference != 0:
            for n_row in range(0, len(sx)):
                for n_column in range(0, len(sx)):
                    ii = n_row + multip_difference
                    jj = n_column + multip_difference
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
        for k, sz_bra in enumerate(sz_values(spin)):
            for g, sz_ket in enumerate(sz_values(spin)):

                if are_equal(sz_bra, sz_ket):
                    s_z[k, g] = sz_ket

                if are_equal(sz_bra, sz_ket + 1):
                    s_x[k, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[k, g] = -0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)

                if are_equal(sz_bra, sz_ket - 1):
                    s_x[k, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[k, g] = 0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)

        return s_x, s_y, s_z

    s2_states = s2_from_file(file)  # s2 of each of the states
    single_s2_values = s2_single_values(s2_states)  # s2 of each of the states, without repetition
    number_of_spins = len(single_s2_values)

    # Spin matrix of non-relativistic states:
    max_multiplicity = int(2 * s2_to_s(max(single_s2_values)) + 1)
    ground_multiplicity = int(2 * s2_to_s(single_s2_values[0]) + 1)

    standard_spin_matrix = np.zeros((ground_multiplicity, ground_multiplicity, 3), dtype=complex)
    spin_matrix = np.zeros((len(n_states) * max_multiplicity, len(n_states) * max_multiplicity, 3), dtype=complex)

    if bolvin == 1:

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

    else:

        for i in range(0, number_of_spins):
            for j in n_states:
                if s2_states[j - 1] == single_s2_values[i]:
                    s = s2_to_s(single_s2_values[i])
                    sx, sy, sz = spin_matrices(s)
                    sx, sy, sz = long_spin_matrices(sx, sy, sz, max_multiplicity, s2_states[j - 1])

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


def get_ground_state_orbital_momentum(file, totalstates):
    """
    Obtaining the orbital angular momentum between the ground state and all excited states.
    :param: file, totalstates
    :return: orbital_momentum
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    word_search = ['< B | Lx | A >', '< B | Ly | A >', '< B | Lz | A >']
    elements = ['0.000', '0.000', '0.000']  # L between ground state and it-self does not exist (x,y,z)

    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[8]
            element = element.replace('i', '')
            elements.append(element)

        if len(elements) == totalstates * 3:
            break  # There are all the Lk values with between ground and excited states

    orbital_momentum = []
    nstate = 0

    for i in range(0, len(elements), 3):
        one_state_momentum = []

        x_orb = abs(float(elements[i]))
        one_state_momentum.append(x_orb)

        y_orb = abs(float(elements[i+1]))
        one_state_momentum.append(y_orb)

        z_orb = abs(float(elements[i+2]))
        one_state_momentum.append(z_orb)

        # Select the orbital angular momentum that is higher, to be shown
        index_max_momentum = one_state_momentum.index(max(one_state_momentum))
        orbital_momentum.append(float(elements[nstate*3 + index_max_momentum]))

        nstate += 1
    return orbital_momentum


def get_orbital_matrices(file, totalstates, selected_states, sz_list, bolvin):
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

    if bolvin == 1:
        doublets_lk = get_doublets_momentum(selected_states, selected_lk)
        all_multip_lk = doublets_lk
    return all_multip_lk


def get_orbital_matrices_manual(file, totalstates, n_states):
    """
    Orbital angular momentum values are written in matrix with 'bra' in rows and 'ket' in columns,
    with spin order -Ms , +Ms. Third dimension is the direction.
    :param: file, totalstates, n_states
    :return: orbital_matrix
    """
    word_search = ['< B | Lx | A >', '< B | Ly | A >', '< B | Lz | A >']
    elements = []

    # take L values from output
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    for line in data:
        if any(i in line for i in word_search):
            element = line[19:32]
            elements.append(element.split())
    # ----------------------------------------------------------------------------------
    # put orbital angular momentum in the array
    angular_selection_list = []

    for ket in n_states:  # | A >
        for bra in n_states:  # < B |

            n_dim = 0
            while n_dim < 3:

                if ket == bra:  # momentum between same state values zero
                    angular_selection_list.append(['0.000000'])
                    # print('eq',ket,bra)

                elif ket < bra:  # In Q-Chem, L written when ket < bra
                    # print('diff', ket, bra)
                    ket_position = 0
                    for i in range(1, ket):
                        ket_position = ket_position + 3 * (totalstates - i)
                    bra_position = 3 * (bra - ket - 1)

                    orbital_position = ket_position + bra_position + n_dim
                    angular_selection_list.append(elements[orbital_position])

                n_dim = n_dim + 1

    angular_selection = np.array(angular_selection_list, dtype=float)
    # ----------------------------------------------------------------------------------------------
    orbital_matrix = np.zeros((len(n_states) * 2, len(n_states) * 2, 3), dtype=complex)
    nel = 0

    for ket in n_states:  # |A, -1/2 >, |A, +1/2 >
        for bra in n_states:  # <B, -1/2|, <B, +1/2|
            ket_index = n_states.index(ket) * 2
            bra_index = n_states.index(bra) * 2

            if ket == bra:
                for k in range(0, 3):
                    # print('eq', bra, ket, '-->', nel + k)
                    orbital_matrix[bra_index, ket_index, k] = complex(0, angular_selection[nel + k])
                    orbital_matrix[bra_index + 1, ket_index + 1, k] = complex(0, angular_selection[nel + k])

                    orbital_matrix[ket_index, bra_index, k] = (-1) * orbital_matrix[bra_index, ket_index, k]
                    orbital_matrix[ket_index + 1, bra_index + 1, k] = (-1) * orbital_matrix[
                        bra_index + 1, ket_index + 1, k]
                nel = nel + 3

            if ket < bra:  # In Q-Chem, SOCs written when ket < bra
                for k in range(0, 3):
                    # print('diff', bra, ket, '-->', nel + k)
                    orbital_matrix[bra_index, ket_index, k] = complex(0, angular_selection[nel + k])
                    orbital_matrix[bra_index + 1, ket_index + 1, k] = complex(0, angular_selection[nel + k])

                    orbital_matrix[ket_index, bra_index, k] = (-1) * orbital_matrix[bra_index, ket_index, k]
                    orbital_matrix[ket_index + 1, bra_index + 1, k] = (-1) * orbital_matrix[
                        bra_index + 1, ket_index + 1, k]
                nel = nel + 3
    return orbital_matrix
