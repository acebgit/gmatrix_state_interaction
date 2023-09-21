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
roberto_molecules/C54H24B4/ms_notnull.out'
file_ms_null = '\
roberto_molecules/C54H24B4/ms_null.out'

selected_states = 0  # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [2, 1, 3, 4, 5]  # States to be included when "selected_states = 0"
# [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
symmetry_selection = 'A2'  # Symmetry selected states
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

#####################################
#      G-TENSOR CALCULATION
#####################################
energies_ms_null, selected_socs_ms_null, sz_list_ms_null, ground_sz, totalstates \
    = gfactor_exchange_energies_socs(file_ms_notnull, file_ms_null, states_ras, soc_options)

g_shift = from_energies_soc_to_g_values(file_ms_null, states_ras, totalstates, energies_ms_null, selected_socs_ms_null, sz_list_ms_null, ground_sz)

print_g_calculation(file_ms_null, totalstates, states_ras, states_ras, g_shift * 1000, symmetry_selection=0)

if excited_states_analysis == 1:
    get_excited_states_analysis(file_ms_null, cutoff=0.9)

if sos_analysis == 1:
    sos_analysis_and_plot(file_ms_null, states_ras, selected_states, order_symmetry=1)
