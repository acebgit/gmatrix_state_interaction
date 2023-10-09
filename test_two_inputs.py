#####################################
#          MODULES SELECTION
#####################################
from parser_mixing_inputs import gfactor_presentation_mixinputs, sos_analysis_and_plot_mixinputs, excited_states_analysis_mixinputs

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 0
excited_states_analysis = 0
sos_analysis = 1
ppm = 1

file_msnull = '\
triplets_molecules/o2_11_9_allmultip.out'
file_ms_notnull = '\
triplets_molecules/o2_11_9_triplets.out'

states_option = 1  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # States to be included when "selected_states = 0"
states_sym = 'A1'

#####################################
#      G-TENSOR CALCULATION
#####################################
if gfactor_two_inputs == 1:
    gfactor_presentation_mixinputs(file_msnull, file_ms_notnull, states_ras, states_option, states_sym, ppm)

if excited_states_analysis == 1:
    excited_states_analysis_mixinputs(file_msnull, file_ms_notnull, states_ras, states_option, states_sym, plots = 1, save_pict=0)

if sos_analysis == 1:
    sos_analysis_and_plot_mixinputs(file_msnull, file_ms_notnull, states_ras, states_option, states_sym, ppm)
