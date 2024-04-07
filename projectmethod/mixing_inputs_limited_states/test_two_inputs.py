#####################################
#          MODULES SELECTION
#####################################
from projectmethod.parsers.parser_gtensor import print_g_calculation
from projectmethod.mixing_inputs_limited_states.parser_mixing_inputs import gfactor_exchange_energies_socs, from_energies_soc_to_g_values, sos_analysis_and_plot
from projectmethod.parsers.parser_excitstates import get_excited_states_analysis

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 0
excited_states_analysis = 1
sos_analysis = 0
ppm = 1

file_ms_notnull = '\
../triplets_molecules/o2_8_6_triplets.out'
file_ms_null = '\
../triplets_molecules/o2_8_6_allmultip.out'

selected_states = 0  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1,9,10]  # States to be included when "states_option = 0"

#####################################
#      G-TENSOR CALCULATION
#####################################
if gfactor_two_inputs == 1:
    energies_ms_null, selected_socs_ms_null, sz_list_ms_null, ground_sz, totalstates \
        = gfactor_exchange_energies_socs(file_ms_notnull, file_ms_null, states_ras, selected_states)

    g_shift = from_energies_soc_to_g_values(file_ms_null, states_ras, totalstates, energies_ms_null, selected_socs_ms_null, sz_list_ms_null, ground_sz)

    print_g_calculation(file_ms_null, totalstates, states_ras, states_ras, g_shift, ppm, symmetry_selection=0)

if excited_states_analysis == 1:
    cut_off = 0.9
    get_excited_states_analysis(file_ms_null, selected_states, states_ras, cut_off, plots=1, save_pict=0)

if sos_analysis == 1:
    sos_analysis_and_plot(file_ms_notnull, file_ms_null, states_ras, selected_states, ppm, order_symmetry=0)
