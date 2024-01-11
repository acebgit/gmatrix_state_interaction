import json
import numpy as np
from procedure_projection_technique.parsers.parser_gtensor import get_number_of_states, get_hamiltonian_construction, diagonalization

file = '../../generalization_projection_technique/test/benzophenone_10_7_triplets.out.json'


def extract_list_from_json(filee):
    """
    Get lists from the information from "json" output file.
    :param filee:
    :return: total_energy, excitation_energy, spin, soc, orbital_momentum
    """
    with open(filee, 'r') as f:
        object_text = f.read()
    input_dict = json.loads(object_text)

    total_energy = []
    excitation_energy = []
    for i in input_dict['selected_energy_dict']:
        total_energy.append(float(input_dict['selected_energy_dict'][i][0]))
        excitation_energy.append(float(input_dict['selected_energy_dict'][i][1]))

    spin = []
    for i in input_dict['spin_dict']:
        spin.append(int(float(input_dict['spin_dict'][i])))

    soc = []
    for i in input_dict['soclist_dict']:
        soc.append((input_dict['soclist_dict'][i]))

    orbital_momentum = []
    for i in input_dict['orbitalmomentlist_dict']:
        orbital_momentum.append([complex(input_dict['orbitalmomentlist_dict'][i][0]),
                                 complex(input_dict['orbitalmomentlist_dict'][i][1]),
                                 complex(input_dict['orbitalmomentlist_dict'][i][2])])
    return total_energy, excitation_energy, spin, soc, orbital_momentum


def s2_to_s(s2):
    """
    get total spin (s) from s^2
    :param: s2_all
    :return: total spin (s)
    """
    return int(0.5 * (-1 + np.sqrt(1 + 4 * s2)))

total_energies, excitation_energies, spin_states, soc_states, orbitmoment_states = extract_list_from_json(file)

max_spin = max(spin_states)
max_multip = s2_to_s(max_spin)
max_sz_list = list(range(-max_multip, max_multip+1))
print(max_sz_list)
exit()

hamiltonian = get_hamiltonian_construction(excitation_energies, selected_socs, sz_list)


exit()

eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

spin_matrix, standard_spin_matrix = get_spin_matrices(file, states_ras)

orbital_matrix = get_orbital_matrices(file, totalstates, states_ras, sz_list)

combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, sz_list)

combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, sz_list)

g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                               sz_list, sz_ground, ppm)

print_g_calculation(file, totalstates, states_option, states_ras, g_shift, symmetry_selection)


