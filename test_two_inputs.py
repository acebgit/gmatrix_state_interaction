#####################################
#          MODULES SELECTION
#####################################
from parser_mixing_inputs import gfactor_exchange_energies_socs, from_energies_soc_to_g_values, \
    print_g_calculation, sos_analysis_and_plot
from parser_excitstates import get_excited_states_analysis

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 1
excited_states_analysis = 0
sos_analysis = 0

file_ms_notnull = '\
triplets_molecules/o2_11_9_triplets.out'
file_ms_null = '\
triplets_molecules/o2_11_9_allmultip.out'

selected_states = 0  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1, 2, 3, 4, 5]  # States to be included when "selected_states = 0"

#####################################
#      G-TENSOR CALCULATION
#####################################
energies_ms_null, selected_socs_ms_null, sz_list_ms_null, ground_sz, totalstates \
    = gfactor_exchange_energies_socs(file_ms_notnull, file_ms_null, states_ras, selected_states)

g_shift = from_energies_soc_to_g_values(file_ms_null, states_ras, totalstates, energies_ms_null, selected_socs_ms_null, sz_list_ms_null, ground_sz)

print_g_calculation(file_ms_null, totalstates, states_ras, states_ras, g_shift * 1000, symmetry_selection=0)

if excited_states_analysis == 1:
    get_excited_states_analysis(file_ms_null, cutoff=0.9)

if sos_analysis == 1:
    sos_analysis_and_plot(file_ms_notnull, file_ms_null, states_ras, selected_states, order_symmetry=1)
