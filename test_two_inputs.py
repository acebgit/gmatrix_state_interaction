#####################################
#          MODULES SELECTION
#####################################
import numpy as np

from parser_gtensor import get_hamiltonian_construction, diagonalization, angular_matrixes_obtention, \
    g_factor_calculation, print_g_calculation
from parser_mixing_inputs import mapping_between_states, get_input_values, totalstates_mix, eigenenergy_mix, socs_mix, \
    angular_momentums_mix, sos_analysis_and_plot
from parser_excitstates import get_excited_states_analysis

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 1
excited_states_analysis = 0
sos_analysis = 0
ppm = 1

file_msnull = '\
triplets_molecules/o2_11_9_allmultip.out'
file_ms_notnull = '\
triplets_molecules/o2_11_9_triplets.out'

states_option = 1  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [2, 1, 3, 4, 5,6,7,8,9,10]  # States to be included when "selected_states = 0"
states_sym = 'A1'

#####################################
#      G-TENSOR CALCULATION
#####################################
dict_mapping, list_mapping = mapping_between_states(file_msnull, file_ms_notnull, states_ras, states_option, states_sym)

totalstates_1, states_ras_1, eigenenergies_ras_1, selected_socs_1, sz_list_1, sz_ground_1, \
spin_matrix_1, standard_spin_matrix_1, orbital_matrix_1 = get_input_values(file_msnull, states_ras, states_option,
                                                                                                             states_sym, soc_options=0)

totalstates_2, states_ras_2, eigenenergies_ras_2, selected_socs_2, sz_list_2, sz_ground_2, \
spin_matrix_2, standard_spin_matrix_2, orbital_matrix_2= get_input_values(file_ms_notnull, states_ras, states_option,
                                                                                                             states_sym, soc_options=0)

totalstates = totalstates_mix(states_ras_1, states_ras_1, list_mapping)

eigenenergy = eigenenergy_mix(eigenenergies_ras_1, eigenenergies_ras_2, list_mapping)

socs = socs_mix(selected_socs_1, selected_socs_2, list_mapping, sz_list_1, totalstates)

hamiltonian = get_hamiltonian_construction(states_ras, eigenenergy, socs, sz_list_1)
# print('Hamiltonian:')
# print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
#                  for row in np.round((hamiltonian[:, :] * 219474.63068), 5)]))  # * 219474.63068
# print('---')

eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

spin_matrix = angular_momentums_mix(spin_matrix_1, spin_matrix_2, list_mapping, sz_list_1, totalstates)

orbital_matrix = angular_momentums_mix(orbital_matrix_1, orbital_matrix_2, list_mapping, sz_list_1, totalstates)

combination_spin_matrix = angular_matrixes_obtention(eigenvector, spin_matrix, sz_list_1)

combination_orbital_matrix = angular_matrixes_obtention(eigenvector, orbital_matrix, sz_list_1)

g_shift = g_factor_calculation(standard_spin_matrix_1, combination_spin_matrix, combination_orbital_matrix,
                               sz_list_1, sz_ground_1)

print_g_calculation(file_msnull, totalstates_1, states_option, states_ras, g_shift, ppm, states_sym)

if excited_states_analysis == 1:
    get_excited_states_analysis(file_ms_notnull, cutoff=0.9)

if sos_analysis == 1:
    sos_analysis_and_plot(file_msnull, file_ms_notnull, states_ras, states_option, states_sym)
