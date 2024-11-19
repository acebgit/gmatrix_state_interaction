"""
Calculation of the g-tensor using Q-Chem output with RAS-CI
"""
__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np
from numpy import linalg, sqrt
from scipy import constants
import pandas as pd
from projection_method.parsers.gtensor_calculation import get_hamiltonian_construction, hermitian_test, \
    units_gshift, diagonalization, angular_matrices_obtention

# Get the absolute path to local folder 'PyQChem'
import sys 
import os 
parent_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), 'C:\\Users\\HP.LAPTOP-F127N3L3\\Desktop\\my_programs\\PyQChem'))
if parent_directory not in sys.path: # Add 'folder2' to the Python path
    sys.path.append(parent_directory)
from pyqchem.parsers.parser_rasci import parser_rasci

lande_factor = 2.002319304363

def get_number_of_states(file):
    """
    Obtain the total number of states selected in RAS-CI in Q-Chem output.
    :param: file
    :return: totalstates
    """   
    with open(file, encoding="utf8") as file:
        for line in file:
            if 'RAS-CI total energy for state' in line:
                totalstates = int(line.split(':')[0].split()[5])
    return totalstates


def get_symmetry_states(file, totalstates):
    """
    Return two lists: i) Symmetry of each state (A1,A2,A3...), ii) Symmetries ordered (1A1,2A2,3A3...)
    :param: file, totalstates
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
            element = element.replace('*', '')
            elements.append(element)
        if len(elements) == totalstates:
            break
    all_state_symmetries = np.array(elements)

    ordered_state_symmetries = []
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
    Depending on "states_option" it returns states: 0) in "state_ras" 1) all states selected 2) States selected by
    selected symmetry "symmetry_selection"
    :param: file, totalstates, selected_states, states_option, symmetry_selection
    :return: selected_states
    """
    all_symmetries, ordered_symmetries = get_symmetry_states(file, totalstates)

    if states_option == 0:
        for i in selected_states:
            if i <= 0 or i > totalstates:
                raise ValueError("The number of states selected selected must be among the total number of states "
                                 "selected calculated in QChem. "
                                 "Select a different number of states selected")

    elif states_option == 1:
        selected_states = list(range(1, totalstates + 1))

    elif states_option == 2:
        states_selected_by_symmetry = [1]  # ground state is always added

        for nstate in range(1, totalstates + 1):
            if (all_symmetries[nstate - 1] == symmetry_selection) and (nstate != 1):
                states_selected_by_symmetry.append(nstate)
        selected_states = states_selected_by_symmetry

        if states_selected_by_symmetry == [1]:
            raise ValueError("There is not this symmetry. Change the symmetry selection.")
    return selected_states


def take_selected_states_values(input_list, selected_states):
    """
    From "input_list", take only the elements in "selected_states" positions and put them in "output_list"
    :param: input_list, selected_states
    :return: output_list
    """
    output_list = []
    for i in selected_states:
        output_list.append(input_list[i - 1])
    return output_list


def get_eigenenergies(file, totalstates, selected_states):
    """
    Get energies of RAS-CI states selected in Q-Chem output. Energies in a.u., excitation energies in eV.
    :param: file, totalstates, selected_states
    :return: eigenenergies, excitation_energies
    """
    word_search = 'RAS-CI total energy for state'
    element_energy = []
    elements_excitenergy = []

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                element_energy.append(float(line.split(":")[1]))
                next_line = next(file)

                if "*******" in next_line.split("=")[1]:
                    ener = (float(line.split(":")[1]) - element_energy[0]) * constants.physical_constants['Hartree energy in eV'][0]
                    elements_excitenergy.append(ener)
                else:
                    elements_excitenergy.append(float(next_line.split("=")[1]))


            if len(element_energy) == totalstates:
                break

    energies_selected = take_selected_states_values(element_energy, selected_states)
    excited_energies_selected = take_selected_states_values(elements_excitenergy, selected_states)

    eigenenergies = np.array(energies_selected)
    excitation_energies_ev = np.array(excited_energies_selected)
    excitation_energies = excitation_energies_ev / constants.physical_constants['Hartree energy in eV'][0]
    return eigenenergies, excitation_energies


def get_spin_orbit_couplings(output, totalstates, selected_states, soc_option, soc_order):
    """
    Get: i) Spin-orbit matrix with dimensions 'bra' x 'ket', with spin order (-Ms , +Ms) in the order of
    "selected_states", ii) [-sz,..,+sz] of the highest s2 state, iii) [-sz,..,+sz] of ground state.
    :param: file, totalstates, selected_states, soc_option.
    :return: selected_socs, sz_max_allstates, sz_ground_state.
    """

    def select_soc_search(soc_options):
        soc_searchs = 'total_soc_mat'
        if soc_options == 1:
            soc_searchs = '1e_soc_mat'
        elif soc_options == 2:
            soc_searchs = '2e_soc_mat'
        else:
            pass
        return soc_searchs

    def get_states_sz(qchem_file, states_selected):
        """
        Get a list of: i) s2 of all the states ii) [-sz,..,+sz] of the highest s2 state iii) [-sz,..,+sz] of ground
        state.
        :param: qchem_file, states_selected
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

        all_sz = from_s2_to_sz_list(max(all_multip))
        ground_sz = from_s2_to_sz_list(ordered_multip[0])

        if len(ground_sz) == 1:
            raise ValueError("Warning! It is not allowed the calculation of the g-tensor in a singlet ground state. "
                             "Ground state corresponds to the first of the included states selected.")
        return all_multip, all_sz, ground_sz

    def get_all_socs(line, n_states, multiplicities, all_sz, soc_selection):
        """
        Get SOC matrix with spin order (-Ms , +Ms) in the order of "selected_states". If there is no value in that
        |I,S,Ms> state, SOC = 0.
        :param: line, n_states, multiplicities, all_sz, soc_selection
        :return: all_soc
        """
        all_soc = np.zeros((n_states * len(all_sz), n_states * len(all_sz)), dtype=complex)

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
        return all_soc

    def get_selected_states_socs(n_states, all_sz, socs, soc_orders):
        """
        Get SOC matrix of the selected states "n_states" from all the socs "all_sz".
        :param: n_states, all_sz, socs
        :return: selected_soc
        """
        selected_soc = np.zeros((len(n_states) * len(all_sz), len(n_states) * len(all_sz)), dtype=complex)

        if soc_orders == 0: # include all orders
            for i, all_i in enumerate(n_states):
                for j, all_j in enumerate(n_states):
                    for sz_1 in range(0, len(all_sz)):
                        for sz_2 in range(0, len(all_sz)):
                            i_index = i * len(all_sz) + sz_1
                            j_index = j * len(all_sz) + sz_2
                            all_i_index = (all_i - 1) * len(all_sz) + sz_1
                            all_j_index = (all_j - 1) * len(all_sz) + sz_2
                            selected_soc[i_index][j_index] = socs[all_i_index][all_j_index]

        elif soc_orders == 1: # include only first orders 
            for i, all_i in enumerate(n_states): # Fist row
                for j, all_j in enumerate(n_states):
                    if i==0 or j==0:
                        # print('i', i, 'j', j)
                        for sz_1 in range(0, len(all_sz)):
                            for sz_2 in range(0, len(all_sz)):
                                i_index = i * len(all_sz) + sz_1
                                j_index = j * len(all_sz) + sz_2
                                all_i_index = (all_i - 1) * len(all_sz) + sz_1
                                all_j_index = (all_j - 1) * len(all_sz) + sz_2
                                selected_soc[i_index][j_index] = socs[all_i_index][all_j_index]

        elif soc_orders == 2: # include only second orders 
            for i, all_i in enumerate(n_states): # Fist row
                for j, all_j in enumerate(n_states):
                    if i!=0 and j!=0:
                        # print('i', i, 'j', j)
                        for sz_1 in range(0, len(all_sz)):
                            for sz_2 in range(0, len(all_sz)):
                                i_index = i * len(all_sz) + sz_1
                                j_index = j * len(all_sz) + sz_2
                                all_i_index = (all_i - 1) * len(all_sz) + sz_1
                                all_j_index = (all_j - 1) * len(all_sz) + sz_2
                                selected_soc[i_index][j_index] = socs[all_i_index][all_j_index]
        return selected_soc

    data = output['interstate_properties']
    soc_search = select_soc_search(soc_option)

    all_multiplicities, sz_max_allstates, sz_ground_state = get_states_sz(output, selected_states)
    all_socs = get_all_socs(data, totalstates, all_multiplicities, sz_max_allstates, soc_search)
    selected_socs = get_selected_states_socs(selected_states, sz_max_allstates, all_socs, soc_order)
    selected_socs = selected_socs / (constants.physical_constants['hartree-inverse meter relationship'][0] / 100)

    # print('SOC:')
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                 for row in np.round((selected_socs[:,:])*(constants.physical_constants['hartree-inverse meter relationship'][0] / 100))]))
    # exit()
    return selected_socs, sz_max_allstates, sz_ground_state


def get_spin_matrices(file, selected_states):
    """
    Get spin matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the order of
    "selected_state". Functions "s2_to_s", "s_to_s2" and "spin_matrices" are taken from inside PyQChem.
    :param: file, selected_state.
    :return: spin_matr, standard_spin_mat.
    """

    def get_all_s2(file_qchem, n_states):
        """
        Get i) list with s2 of all states ii) list with s2 all of selected states
        :param: file_qchem, n_states
        :return: s2_all, s2_selected_list
        """
        search = ['  <S^2>      : ']

        # Take s2 of all the states
        s2_all_list = []
        with open(file_qchem, encoding="utf8") as file_qchem:
            for line in file_qchem:
                if any(a in line for a in search):
                    line = line.split()
                    element = line[2]
                    s2_all_list.append(element)
        s2_all = np.array(s2_all_list, dtype=float)

        # Make a list of dictionaries with s2 of selected states
        s2_selected_list = []
        for i in n_states:
            s2_selected_dict = {'st': i, 's2': s2_all_list[i - 1]}
            s2_selected_list.append(s2_selected_dict)

        if s2_selected_list[0]['s2'] == 0:
            raise ValueError("Warning! It is not allowed the calculation of the g-tensor in a singlet ground st. "
                             "Ground st corresponds to the first of the included states selected.")
        return s2_all, s2_selected_list

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

    def spin_matrices(spin):
        """
        Get spin matrices s_x, s_y, s_z between two spin states selected (s,m) and (s,m') such that
        sxx = < m' | s_x | m >, syy = < m' | s_y | m > and szz = < m' | s_z | m >
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

    def expand_spin_matrices(s_x, s_y, s_z, max_mult, state_mult):
        """
        Expand (sxx, syy, szz) matrices with dimension of the st multiplicity to (sxx, syy, szz) with dimension of the
        maximum multiplicity of states.
        :param: s_x, s_y, s_z, max_mult, state_mult
        :return: long_sx, long_sy, long_sz
        """
        long_sx = np.zeros((max_mult, max_mult), dtype=complex)
        long_sy = np.zeros((max_mult, max_mult), dtype=complex)
        long_sz = np.zeros((max_mult, max_mult), dtype=complex)

        multipl_difference = int((max_mult - state_mult) // 2)

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
        return long_sx, long_sy, long_sz

    def form_big_spin_matrix(st, selected_state, max_mult, sxx, syy, szz, spin_matr):
        """
        Get big spin matrix with dimensions "(len(selected_state) * max_mult, len(selected_state) * max_mult, 3)"
        with s_x, s_y, s_z.
        :param: st, selected_state, max_mult, sxx, syy, szz, spin_matr
        :return: spin_matr
        """
        s_dim = 0
        total_dim = selected_state.index(selected_state[st]) * max_mult
        for row in range(0, max_mult):
            for column in range(0, max_mult):
                spin_matr[total_dim, total_dim, 0] = sxx[s_dim, s_dim]
                spin_matr[total_dim, total_dim + column, 0] = sxx[s_dim, s_dim + column]
                spin_matr[total_dim + row, total_dim, 0] = sxx[s_dim + row, s_dim]
                spin_matr[total_dim + row, total_dim + column, 0] = sxx[s_dim + row, s_dim + column]

                spin_matr[total_dim, total_dim, 1] = syy[s_dim, s_dim]
                spin_matr[total_dim, total_dim + column, 1] = syy[s_dim, s_dim + column]
                spin_matr[total_dim + row, total_dim, 1] = syy[s_dim + row, s_dim]
                spin_matr[total_dim + row, total_dim + column, 1] = syy[s_dim + row, s_dim + column]

                spin_matr[total_dim, total_dim, 2] = szz[s_dim, s_dim]
                spin_matr[total_dim, total_dim + column, 2] = szz[s_dim, s_dim + column]
                spin_matr[total_dim + row, total_dim, 2] = szz[s_dim + row, s_dim]
                spin_matr[total_dim + row, total_dim + column, 2] = szz[s_dim + row, s_dim + column]
        return spin_matr

    def get_standard_spin_matrix(s2_selected_dict, max_mult, spin_mat):
        """
        Construct Standard Spin matrix from the previous Spin Matrix.
        :param: s2_selected_dict, max_mult, spin_mat
        :return: standard_spin_mat
        """
        s2_ground = (float(s2_selected_dict[0]['s2']))
        ground_multip = int(2 * s2_to_s(s2_ground) + 1)
        standard_spin_mat = np.zeros((ground_multip, ground_multip, 3), dtype=complex)

        multip_difference = (max_mult - ground_multip) // 2
        for k in range(0, 3):
            for ii in range(0, ground_multip):
                for jj in range(0, ground_multip):
                    standard_spin_mat[ii, jj, k] = spin_mat[ii + multip_difference, jj + multip_difference, k]
        return standard_spin_mat

    # Take s2 of all states and each of the selected states.
    s2_all_states, s2_selected_dicts = get_all_s2(file, selected_states)

    # Take maximum multiplicity, which determines the dimension of the spin matrix (max_mult x max_mult x 3).
    max_multip = int(2 * s2_to_s(max(s2_all_states)) + 1)
    spin_matrix = np.zeros((len(selected_states) * max_multip, len(selected_states) * max_multip, 3), dtype=complex)

    for state in range(0, len(selected_states)):
        # Form the spin matrix of the state
        s2_state = (float(s2_selected_dicts[state]['s2']))
        s = s2_to_s(s2_state)
        sx, sy, sz = spin_matrices(s)

        # Expand the spin matrix of the st to the dimension of the maximum multiplicity (max_mult)
        state_multip = int(2 * s2_to_s(s2_state) + 1)
        sx, sy, sz = expand_spin_matrices(sx, sy, sz, max_multip, state_multip)

        # Mix (sxx,syy,szz) in one spin matrix:
        spin_matrix = form_big_spin_matrix(state, selected_states, max_multip, sx, sy, sz, spin_matrix)

    # Take Standard Spin Matrix from Spin Matrix
    standard_spin_matrix = get_standard_spin_matrix(s2_selected_dicts, max_multip, spin_matrix)

    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #                 for row in np.round((spin_matrix[:,:,0]))]))

    return spin_matrix, standard_spin_matrix


def get_orbital_matrices(output, totalstates, selected_states, sz_list):
    """
    Get orbitals angular momentum matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the
    order of "selected_state".
    :param: file, totalstates, selected_states, sz_list
    :return: all_multip_lk
    """

    def get_all_momentum(dataa, n_states):
        """
        Get Lk between all the states selected in x,y,z dimensions. | A > in columns, < B | in rows.
        :param: line, selected_states
        :return: all_orbital_momentum
        """
        all_orbital_momentum = np.zeros((n_states, n_states, 3), dtype=complex)

        for i in range(0, n_states):
            for j in range(0, n_states):
                if i != j:
                    element = dataa[(i + 1, j + 1)]['angular_momentum']

                    for k in range(0, 3):
                        all_orbital_momentum[j, i, k] = element[k]
        return all_orbital_momentum

    def get_selected_states_momentum(n_states, all_momentum):
        """
        Get Lk between the selected states selected in x,y,z dimensions.
        :param: selected_states, all_momentum
        :return: selected_momentum
        """
        selected_momentum = np.zeros((len(n_states), len(n_states), 3), dtype=complex)

        for k in range(0, 3):
            for i, all_i in enumerate(n_states):
                for j, all_j in enumerate(n_states):
                    selected_momentum[i][j][k] = all_momentum[all_i - 1][all_j - 1][k]
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

    data = output['interstate_properties']

    all_lk = get_all_momentum(data, totalstates)
    selected_lk = get_selected_states_momentum(selected_states, all_lk)
    all_multip_lk = get_all_multip_momentum(selected_lk, sz_list)
    return all_multip_lk


def from_qchem_to_sto(g_shift):
    """
        Change orientation  from Q-Chem to Standard Nuclear Orientation
        :param g_shift:
        :return:
        """
    a = g_shift[0]
    g_shift[0] = g_shift[1]
    g_shift[1] = a
    return g_shift


def from_angmoments_to_gshifts(standard_spin_matrix, s_matrix, l_matrix, sz_list, ground_sz, ppms=None):
    """
    g-shift with orbitals and spin angular momentum matrices.
    :param: standard_spin_matrix, s_matrix, l_matrix, sz_list, ground_sz
    :return: g_shift
    """

    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                 for row in np.round((standard_spin_matrix[:,:,0]))]))

    def j_matrix_formation(spin, orbital, list_sz, sz_ground):
        """
        Get the total angular momentum matrix from the orbitals and spin angular momentums. Then, expand it to the
        multiplicity of the ground state multiplicity "sz_ground".
        :param: spin, orbitals, list_sz, sz_ground
        :return: j_matr
        """
        j_big_matrix = lande_factor * spin + orbital

        sz_difference = (len(list_sz) - len(sz_ground)) // 2
        j_matr = np.zeros((len(sz_ground), len(sz_ground), 3), dtype=complex)
        for k in range(0, 3):
            for ii in range(0, len(j_matr)):
                for j in range(0, len(j_matr)):
                    j_matr[ii, j, k] = j_big_matrix[ii + sz_difference, j + sz_difference, k]

            hermitian_test(j_matr[:, :, k], list_sz)
        return j_matr

    def j_diagonalization(matrix):
        """
        J-matrix diagonalization giving i) eigenvectors ii) diagonal matrix.
        :param: matrix
        :return: eigenvalues, eigenvectors, diagonal_matrix
        """
        eigenvalues, eigenvectors = linalg.eigh(matrix)
        rotation_inverse = np.linalg.inv(eigenvectors)
        diagonal_matrix = np.matmul(np.matmul(rotation_inverse, matrix), eigenvectors)
        return eigenvectors, diagonal_matrix

    def trace_g_values(total_j, spin_matr):
        """
        Obtaining g_value with projection of the traces.
        :param: total_j, spin_matr
        :return: g_value
        """
        a = np.matmul(total_j, spin_matr)
        b = np.matmul(spin_matr, spin_matr)
        g_value = np.trace(a) / np.trace(b)
        return g_value

    j_matrix = j_matrix_formation(s_matrix, l_matrix, sz_list, ground_sz)

    # PROJECTION TECHNIQUE TO OBTAIN THE TRIANGULAR G-MATRIX:
    # 1) g-value zz
    j_matrix_rotation, j_matrix_diagonal_z = j_diagonalization(j_matrix[:, :, 2])
    g_matrix_triangular = np.zeros((3, 3), dtype=complex)
    g_matrix_triangular[2, 2] = trace_g_values(j_matrix_diagonal_z, standard_spin_matrix[:, :, 2])

    # 2) pseudospin z
    pseudospin_matrix = np.zeros((len(j_matrix[0, :]), len(j_matrix[:, 0]), 3), dtype=complex)
    pseudospin_matrix[:, :, 2] = j_matrix[:, :, 2] / g_matrix_triangular[2, 2]

    # 3) g-value yz
    g_matrix_triangular[1, 2] = trace_g_values(j_matrix[:, :, 1], pseudospin_matrix[:, :, 2])
    residue = j_matrix[:, :, 1] - g_matrix_triangular[1, 2] * pseudospin_matrix[:, :, 2]

    # 4) g-value yy
    j_matrix_rotation, j_matrix_transformed_y = j_diagonalization(j_matrix[:, :, 1])
    g_matrix_triangular[1, 1] = trace_g_values(j_matrix_transformed_y, standard_spin_matrix[:, :, 2])

    # 5) pseudospin y
    pseudospin_matrix[:, :, 1] = residue[:, :] / g_matrix_triangular[1, 1]

    # 6) g-value xy and xz
    g_matrix_triangular[0, 1] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 1])
    g_matrix_triangular[0, 2] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 2])

    # 7) g-value xx
    j_matrix_rotation, j_matrix_transformed_x = j_diagonalization(j_matrix[:, :, 0])
    g_matrix_triangular[0, 0] = trace_g_values(j_matrix_transformed_x, standard_spin_matrix[:, :, 2])

    # 8) from g_matrix_triangular to g-shifts
    g_matrix = np.matmul(g_matrix_triangular, np.transpose(g_matrix_triangular))
    g_matrix_eigenvalues, rotation_g_matrix, g_matrix_diagonal = diagonalization(g_matrix)
    # print('Lande factor: ', lande_factor)
    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #             for row in np.round((g_matrix[:,:]), 10)]))
    # print("")
    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #             for row in np.round((g_matrix[:,:]), 10)]))
    # exit()

    g_shifts = np.zeros(3, dtype=complex)
    for i in range(0, 3):
        g_shifts[i] = (sqrt(g_matrix_diagonal[i, i]) - lande_factor) * 1000

    g_shifts = units_gshift(g_shifts, ppms)
    # g_shifts = from_qchem_to_sto(g_shifts)

    return g_matrix, g_shifts


def from_energies_soc_to_g_values(file, states_ras, totalstates,
                                  excitation_energies_ras, soc_ras, sz_list, ground_sz, ppm):
    """"
    Obtention of the g-values from the eigenenergies and the SOCs.
    :param:file_ms_notnull, states_msnull, states_option, excitation_energies_ras, soc_ras, list_sz, ground_sz
    :return: g_shift
    """
    with open(file, encoding="utf8") as f:
        output = f.read()
    output_parsered = parser_rasci(output)

    hamiltonian_ras = get_hamiltonian_construction(excitation_energies_ras, soc_ras, sz_list)

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian_ras)

    spin_matrix, standard_spin_matrix = get_spin_matrices(file, states_ras)

    orbital_matrix = get_orbital_matrices(output_parsered, totalstates, states_ras, sz_list)

    combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, sz_list)

    combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, sz_list)

    gmatrix, gshift = from_angmoments_to_gshifts(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                                   sz_list, ground_sz, ppm)

    return gshift


def print_g_calculation(file, totalstates, selected_states,
                        states_ras, g_shifts, g_matrix, symmetry_selection, soc_orders, soc_option):
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

    soc_options_dict = {0: "All BP Hamiltonian", 1: "One electron term", 2: "Bielectornic term"}
    soc_order_dict = {0: "All orders", 1: "First-order", 2: "Second-order"}
    print("SOC: ", soc_order_dict[soc_orders], ',', soc_options_dict[soc_option])

    print(" ")
    print("------------------------")
    print(" RAS-CI RESULTS")
    print("------------------------")
    print('g-factor (x y z dimensions):')
    print(np.round(g_shifts[0].real, 3), np.round(g_shifts[1].real, 3),
          np.round(g_shifts[2].real, 3))
    # print('')
    # print('g-matrix:')
    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row]) 
    #                  for row in np.round((g_matrix[:,:]), 4)]))


def get_states_notnull_soc(output, totalstatess, states_initial):
    """
    Take the initial states and create a list with all the states that have not null SOC 
    with the first state, that is the considered ground state. 
    :param: states_initial
    :return: states_notnull_soc
    """
    data = output['interstate_properties']
    for i in range(0, totalstatess):
        for j in range(0, totalstatess):
                element = data[(i + 1, j + 1)]['mf_socc']
                print(element)

    # # Take SOCCs between the considered ground state and all other states.
    # ground_state_soccs = []
    # initial_socc = (states_selected[0]-1) * totalstates
    # final_socc = (states_selected[0]-1) * totalstates + totalstates
    # for i in range(initial_socc, final_socc):
    #     ground_state_soccs.append(all_soccs[i])

    # # Take SOCCs between the considered ground state and the selected states.
    # socc_list = take_selected_states_values(ground_state_soccs, states_selected)
    # socc = np.array(socc_list, dtype=float)
    # return socc


def gfactor_presentation(ras_input, states_ras, states_option, symmetry_selection, soc_options, soc_order, ppm):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, states_msnull, states_option, symmetry_selection, soc_options
    :return: g-shifts
    """
    with open(ras_input, encoding="utf8") as f:
        output = f.read()
    output_parsered = parser_rasci(output)

    totalstates = get_number_of_states(ras_input)

    states_selected = get_selected_states(ras_input, totalstates, states_ras, states_option, symmetry_selection)

    # states_selected = get_states_notnull_soc(output_parsered, totalstates, states_selected)

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_selected)

    selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(output_parsered, totalstates, states_selected, soc_options, soc_order)

    hamiltonian_ras = get_hamiltonian_construction(excitation_energies_ras, selected_socs, sz_list)

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian_ras)

    # def change_eigenvectors(eignvector):
    #     row1 = 0
    #     row2 = 1
    #     for column in range(0, len(eigenvector[:,:]), len(sz_ground)):
    #         a = eignvector[:, column+row1].copy()
    #         eignvector[:, column+row1] = eignvector[:, column+row2]
    #         eignvector[:, column+row2] = a
    #     exit()
        
    #     # pd.set_option('display.max_rows', None)
    #     # pd.set_option('display.max_columns', None)
    #     # df = pd.DataFrame(eignvector)
    #     # # df = pd.DataFrame(eigenvector, index=states_list,
    #     #         # columns=['state', 'config.', 'symmetry', 'unpaired orb', 'gxx', 'gyy', 'gzz'])    
    #     # df_rounded = df.round(2)
    #     # print(df)
    #     # exit()

    #     for column in range(len(eignvector)):
    #         rounded_col = [complex(round(num.real, 2), round(num.imag, 2)) for num in eignvector[:, column]]
    #         print('State ', (column//len(sz_ground))+1, 'substate', column%len(sz_ground)+1)
    #         print(rounded_col)
    #         print()
    #     exit()
    #     return eignvector

    # eigenvector = change_eigenvectors(eigenvector)

    spin_matrix, standard_spin_matrix = get_spin_matrices(ras_input, states_selected)

    orbital_matrix = get_orbital_matrices(output_parsered, totalstates, states_selected, sz_list)

    combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, sz_list)

    combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, sz_list)

    gmatrix, gshift = from_angmoments_to_gshifts(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                                   sz_list, sz_ground, ppm)

    print_g_calculation(ras_input, totalstates, states_option, states_selected, gshift, gmatrix, symmetry_selection, soc_order, soc_options)


def from_gvalue_to_shift(lista):
    """
    Obtain the g-shifts from the g-values
    :param: list
    :return: g_shift
    """
    g_shift = []
    for i in range(0, len(lista)):
        value = (lista[i] - lande_factor) * 10**3
        g_shift.append(value)
    print(np.round(g_shift, 2))
