#####################################
#          MODULES SELECTION
#####################################
import numpy as np

from parser_gtensor import get_hamiltonian_construction, diagonalization
from parser_mixing_inputs import mapping_between_states, get_input_values, totalstates_mix, eigenenergy_mix, socs_mix, \
    angular_momentums_mix
from parser_excitstates import get_excited_states_analysis

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 1
excited_states_analysis = 0
sos_analysis = 0

file_msnull = '\
roberto_molecules/C54H24B4/ms_null.out'
file_ms_notnull = '\
roberto_molecules/C54H24B4/ms_notnull.out'

states_option = 0  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [2, 1, 3, 4, 5]  # States to be included when "selected_states = 0"
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
print('Hamiltonian:')
print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
                 for row in np.round((hamiltonian[:, :] * 219474.63068), 5)]))  # * 219474.63068
print('---')

eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

spin_matrix = angular_momentums_mix(spin_matrix_1, spin_matrix_2, sz_list_1, totalstates)

exit()

energies_ms_2, selected_socs_ms_2, sz_list_ms_2, ground_sz, totalstates \
    = gfactor_exchange_energies_socs(file_msnull, file_ms_notnull, states_ras, selected_states)

g_shift = from_energies_soc_to_g_values(file_ms_notnull, states_ras, totalstates, energies_ms_2, selected_socs_ms_2, sz_list_ms_2, ground_sz)

print_g_calculation(file_ms_notnull, totalstates, states_ras, states_ras, g_shift * 1000, symmetry_selection=0)

if excited_states_analysis == 1:
    get_excited_states_analysis(file_ms_notnull, cutoff=0.9)

if sos_analysis == 1:
    sos_analysis_and_plot(file_msnull, file_ms_notnull, states_ras, selected_states, order_symmetry=1)
