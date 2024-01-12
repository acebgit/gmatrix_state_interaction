import json
import numpy as np
from scipy import constants
from procedure_projection_technique.parsers.parser_gtensor import get_hamiltonian_construction, diagonalization, \
    angular_matrices_obtention, g_factor_calculation, print_g_calculation

file = '../../generalization_projection_technique/test/ncl_12_9_triplets_diisgdm.out.json'
ppm = 0


def extract_data_from_json(filee):
    """
    Get lists from the information from "json" output file.
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
        excitation_energy_list.append(float(input_dict['selected_energy_dict'][i][1]))

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


def get_max_sz(spin_state):
    def s2_to_s(s2):
        """
        get total spin (s) from s^2
        :param: s2_all
        :return: total spin (s)
        """
        return int(0.5 * (-1 + np.sqrt(1 + 4 * s2)))

    max_spin = max(spin_state)
    max_multip = s2_to_s(max_spin)
    maxsz_list = list(range(-max_multip, max_multip + 1))
    return maxsz_list


def from_soclist_socmatrix(soc_list, maxsz_list):
    """
    Construct the SOC matrix from the SOC list of json file.
    :param soc_list:
    :param maxsz_list:
    :return: soc_matrix
    """
    len_sz = len(maxsz_list)  # dimension determined by the maximum multiplicity
    nstates = len(soc_list) + 1  # +1 is the ground state
    soc_matrix = np.zeros((nstates * len_sz, nstates * len_sz), dtype=complex)

    for i in range(1, nstates):  # ground state does not have soc, and there is no state j
        len_list_sz1 = len(soc_list[i-1])
        len_list_sz2 = len(soc_list[i-1][0])

        for sz_1 in range(0, len_list_sz1):
            for sz_2 in range(0, len_list_sz2):
                matrix_row = sz_1 + ((len_sz-len_list_sz1)//2)
                matrix_col = i * len_sz + sz_2 + ((len_sz-len_list_sz2)//2)
                # print('State', i, 'Matrix:', matrix_row, matrix_col, '--> List:', list_state, list_dim2, list_dim3)
                soc_matrix[matrix_row][matrix_col] = soc_list[i-1][sz_1][sz_2]
                soc_matrix[matrix_col][matrix_row] = np.conj(soc_matrix[matrix_row][matrix_col])

    # print('SOC:')
    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #                 for row in np.round((soc_matrix[:,:]))]))
    # exit()
    soc_matrix = soc_matrix / (constants.physical_constants['hartree-inverse meter relationship'][0]/100)
    return soc_matrix


energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json = extract_data_from_json(file)

max_sz_list = get_max_sz(spin_json)

soc_matrix = from_soclist_socmatrix(soc_json, max_sz_list)

hamiltonian = get_hamiltonian_construction(excitenergies_json, soc_matrix, max_sz_list)

eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

exit()

spin_matrix, standard_spin_matrix = get_spin_matrices(file, states_ras)

orbital_matrix = get_orbital_matrices(file, totalstates, states_ras, max_sz_list)

combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, max_sz_list)

combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, max_sz_list)

g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                               max_sz_list, sz_ground, ppm)

print_g_calculation(file, totalstates, states_option, states_ras, g_shift, symmetry_selection)


