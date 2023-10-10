#####################################
#          MODULES SELECTION
#####################################
from parser_mixing_inputs import gfactor_presentation_mixinputs, sos_analysis_and_plot_mixinputs, excited_states_analysis_mixinputs

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 1
excited_states_analysis = 0
sos_analysis = 0
ppm = 0

file_msnull = '\
triplets_molecules/nf_12_9_allmultip.out'
file_ms_notnull = '\
triplets_molecules/nf_12_9_enerproc_triplets.out'

states_option = 0  # 0: use "state_ras" ; 1: use all states_selected
states_msnull = [1, 4]  # States to be included when "selected_states = 0"
states_msnotnull = [1, 2, 3]

#####################################
#      G-TENSOR CALCULATION
#####################################
if gfactor_two_inputs == 1:
    gfactor_presentation_mixinputs(file_msnull, file_ms_notnull, states_option, states_msnull, states_msnotnull, ppm)

if excited_states_analysis == 1:
    excited_states_analysis_mixinputs(file_msnull, file_ms_notnull, states_msnull, states_option, plots = 0, save_pict=0)

if sos_analysis == 1:
    sos_analysis_and_plot_mixinputs(file_msnull, file_ms_notnull, states_msnull, states_option, ppm)
