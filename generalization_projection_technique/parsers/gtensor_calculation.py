import json
import numpy as np
from scipy import constants
from procedure_projection_technique.parsers.parser_gtensor import get_hamiltonian_construction, diagonalization, \
    angular_matrices_obtention, g_factor_calculation

file = '../../generalization_projection_technique/test/ncl_12_9_triplets_diisgdm.out.json'
ppm = 0


def extract_data_from_json(filee):
    """
    Get lists with the information from "json" output filee.
    :param filee:
    :return: total_energy, excitation_energy_list, spin_list, soc_list, orbital_momentum_list
    """
    with open(filee, 'r') as f:
        object_text = f.read()
    input_dict = json.loads(object_text)

    total_energy_list = []
    excitation_energy_list = []
    for i in input_dict['selected_energy_dict']:
        total_energy_list.append(float(input_dict['selected_energy_dict'][i][0]))
        excit_energy = float(input_dict['selected_energy_dict'][i][1]) / constants.physical_constants['Hartree '
                                                                                                      'energy in eV'][0]
        excitation_energy_list.append(excit_energy)

    spin_list = []
    for i in input_dict['spin_dict']:
        spin_list.append(int(float(input_dict['spin_dict'][i])))

    soc_list = []
    for i in input_dict['soclist_dict']:
        soc_list.append((input_dict['soclist_dict'][i]))

    orbital_momentum_list = []
    for i in input_dict['orbitalmomentlist_dict']:
        orbital_momentum_list.append([complex(input_dict['orbitalmomentlist_dict'][i][0]),
                                      complex(input_dict['orbitalmomentlist_dict'][i][1]),
                                      complex(input_dict['orbitalmomentlist_dict'][i][2])])
    return total_energy_list, excitation_energy_list, spin_list, soc_list, orbital_momentum_list


def s2_to_s(s2):
    """
    get total spin (s) from s^2
    :param: s2_all
    :return: total spin (s)
    """
    return 0.5 * (-1 + np.sqrt(1 + 4 * s2))


def get_input_data(spin_state):
    """
    Get: i) number of states, ii) a list from -sz to sz, where sz = maximum multiplicity,
    iii) a list from -sz to sz, where sz = ground state spin multiplicity
    :param spin_state:
    :return:
    """
    nstate = len(spin_state)

    max_multip = int(s2_to_s(max(spin_state)))
    max_multip_szlist = list(range(-max_multip, max_multip + 1))

    s_ground = int(s2_to_s(spin_state[0]))
    grst_szlist = list(range(-s_ground, s_ground + 1))

    return nstate, max_multip_szlist, grst_szlist


def from_soclist_socmatrix(soc_list, maxsz_list):
    """
    Construct the SOC matrix from the SOC list of json filee.
    :param soc_list:
    :param maxsz_list:
    :return: socmatrix
    """
    len_sz = len(maxsz_list)  # dimension determined by the maximum multiplicity
    nstate = len(soc_list) + 1  # +1 is the ground state
    socmatrix = np.zeros((nstate * len_sz, nstate * len_sz), dtype=complex)

    for i in range(1, nstate):  # ground state does not have soc, and there is no state j
        soc_list_sz_state1 = len(soc_list[i-1][0])
        soc_list_sz_state2 = len(soc_list[i-1])

        for sz_1 in range(0, soc_list_sz_state1):
            for sz_2 in range(0, soc_list_sz_state2):
                matrix_row = sz_1 + ((len_sz-soc_list_sz_state1)//2)
                matrix_col = i * len_sz + sz_2 + ((len_sz-soc_list_sz_state2)//2)

                # print('State', i, 'Matrix:', matrix_row, matrix_col, '--> List:', i-1, sz_2, sz_1)
                value = complex(soc_list[i-1][sz_2][sz_1])
                socmatrix[matrix_row][matrix_col] = np.conj(value)
                socmatrix[matrix_col][matrix_row] = value

    # print('SOC:')
    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #                 for row in np.round((socmatrix[:,:]))]))
    # exit()
    socmatrix = socmatrix / (constants.physical_constants['hartree-inverse meter relationship'][0]/100)
    return socmatrix


def get_spin_matrices(states_spin, maxsz_list):
    """
    Get spin matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the order of
    "selected_state". Functions "s2_to_s", "s_to_s2" and "spin_matrices" are taken from inside PyQChem.
    :param: filee, selected_state.
    :return: spin_matr, standard_spin_mat.
    """
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
        :param: s_x, s_y, s_z, max_sz_listt, state_mult
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

    def form_big_spin_matrix(statee, max_mult, sxx, syy, szz, spin_matr):
        """
        Get big spin matrix with dimensions "(len(selected_state) * max_sz_listt, len(selected_state) *
        max_sz_listt, 3)" with s_x, s_y, s_z.
        :param: st, selected_state, max_sz_listt, sxx, syy, szz, spin_matr
        :return: spin_matr
        """
        s_dim = 0
        initial_pos = statee * max_mult
        for row in range(0, max_mult):
            for column in range(0, max_mult):
                spin_matr[initial_pos, initial_pos, 0] = sxx[s_dim, s_dim]
                spin_matr[initial_pos, initial_pos + column, 0] = sxx[s_dim, s_dim + column]
                spin_matr[initial_pos + row, initial_pos, 0] = sxx[s_dim + row, s_dim]
                spin_matr[initial_pos + row, initial_pos + column, 0] = sxx[s_dim + row, s_dim + column]

                spin_matr[initial_pos, initial_pos, 1] = syy[s_dim, s_dim]
                spin_matr[initial_pos, initial_pos + column, 1] = syy[s_dim, s_dim + column]
                spin_matr[initial_pos + row, initial_pos, 1] = syy[s_dim + row, s_dim]
                spin_matr[initial_pos + row, initial_pos + column, 1] = syy[s_dim + row, s_dim + column]

                spin_matr[initial_pos, initial_pos, 2] = szz[s_dim, s_dim]
                spin_matr[initial_pos, initial_pos + column, 2] = szz[s_dim, s_dim + column]
                spin_matr[initial_pos + row, initial_pos, 2] = szz[s_dim + row, s_dim]
                spin_matr[initial_pos + row, initial_pos + column, 2] = szz[s_dim + row, s_dim + column]
        return spin_matr

    def get_standard_spin_matrix(spin_states, max_sz_listt, spin_mat):
        """
        Construct Standard Spin matrix from the previous Spin Matrix.
        :param: spin_states, max_sz_listt, spin_mat
        :return: standard_spin_mat
        """
        ground_multiplicity = int(2 * s2_to_s(spin_states[0]) + 1)
        standard_spin_mat = np.zeros((ground_multiplicity, ground_multiplicity, 3), dtype=complex)

        multip_difference = (max_sz_listt - ground_multiplicity) // 2
        for k in range(0, 3):
            for ii in range(0, ground_multiplicity):
                for jj in range(0, ground_multiplicity):
                    standard_spin_mat[ii, jj, k] = spin_mat[ii + multip_difference, jj + multip_difference, k]
        return standard_spin_mat

    len_sz = len(maxsz_list)  # dimension determined by the maximum multiplicity
    nstate = len(states_spin)
    spinmatrix = np.zeros((nstate * len_sz, nstate * len_sz, 3), dtype=complex)

    for state in range(0, nstate):
        # Form the spin matrix of the state
        s2_state = states_spin[state]
        s = s2_to_s(s2_state)
        sx, sy, sz = spin_matrices(s)

        # Expand the spin matrix of the st to the dimension of the maximum multiplicity (max_sz_listt)
        state_multip = int(2 * s + 1)
        sx, sy, sz = expand_spin_matrices(sx, sy, sz, len_sz, state_multip)

        # Mix (sxx,syy,szz) in one spin matrix:
        spinmatrix = form_big_spin_matrix(state, len_sz, sx, sy, sz, spinmatrix)

    # Take Standard Spin Matrix from Spin Matrix
    standardspin_matrix = get_standard_spin_matrix(states_spin, len_sz, spinmatrix)

    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #                 for row in np.round((standardspin_matrix[:,:,0]))]))
    # exit()
    return spinmatrix, standardspin_matrix


def get_orbital_matrices(states_momentum, maxsz_list):
    """
    Get orbitals angular momentum matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the
    order of "selected_state".
    :param: filee, totalstates, selected_states, sz_list
    :return: all_multip_lk
    """
    def get_selected_states_momentum(momentum_states):
        """
        Get Lk between the selected states selected in x,y,z dimensions.
        :param: selected_states, all_momentum
        :return: selected_momentum
        """
        n_states = len(momentum_states) + 1  # +1 because of the ground state
        selected_momentum = np.zeros((n_states, n_states, 3), dtype=complex)

        for i in range(1, n_states):
            for ndim in range(0, 3):
                selected_momentum[i][0][ndim] = momentum_states[i-1][ndim]
                selected_momentum[0][i][ndim] = np.conj(momentum_states[i - 1][ndim])
        return selected_momentum

    def get_all_multip_momentum(selected_momentums, max_szlist):
        """
        Get Lk between the selected states selected in x,y,z dimensions for doublets,
        i.e. in each row (column) there are state A,-1/2 and A,+1/2.
        :param: selected_momentums, max_szlist
        :return: big_orbit_matrix
        """
        n_states = len(selected_momentums)
        multip = len(max_szlist)
        big_orbit_matrix = np.zeros((n_states * multip, n_states * multip, 3), dtype=complex)

        for k in range(0, 3):
            for i in range(0, len(selected_momentums) * len(max_szlist), len(max_szlist)):
                for j in range(0, len(selected_momentums) * len(max_szlist), len(max_szlist)):

                    for multip in range(0, len(max_szlist)):
                        big_orbit_matrix[i + multip, j + multip][k] = \
                            selected_momentums[i // len(max_szlist)][j // len(max_szlist)][k]
        return big_orbit_matrix

    selected_lk = get_selected_states_momentum(states_momentum)
    all_multip_lk = get_all_multip_momentum(selected_lk, maxsz_list)

    # for ndim in range(0, 3):
    #     print()
    #     print('\n'.join([''.join(['{:^8}'.format(item) for item in row]) \
    #                      for row in np.round(all_multip_lk[:, :, ndim], 3)]))
    # exit()
    return all_multip_lk


def print_g_calculation(filee, totalstates, upper_g_tensor_results_ras):
    print("--------------------------------------")
    print("     INPUT SECTION")
    print("--------------------------------------")
    print("File selected: ", filee)
    print("Number of states: ", totalstates)
    print('g-factor (x y z dimensions):')
    print(np.round(upper_g_tensor_results_ras[0].real, 3), np.round(upper_g_tensor_results_ras[1].real, 3),
          np.round(upper_g_tensor_results_ras[2].real, 3))
    print('')


energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json = extract_data_from_json(file)

nstates, max_sz_list, sz_ground = get_input_data(spin_json)

soc_matrix = from_soclist_socmatrix(soc_json, max_sz_list)

hamiltonian = get_hamiltonian_construction(excitenergies_json, soc_matrix, max_sz_list)

eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

spin_matrix, standard_spin_matrix = get_spin_matrices(spin_json, max_sz_list)

orbital_matrix = get_orbital_matrices(orbitmoment_json, max_sz_list)

combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, max_sz_list)

combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, max_sz_list)

g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                               max_sz_list, sz_ground, ppm)

print_g_calculation(file, nstates, g_shift)
