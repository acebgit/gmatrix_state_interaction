#####################################
#          MODULES SELECTION
#####################################
from procedure_projection_technique.parsers.parser_mixing_inputs import gfactor_presentation_mixinputs, sos_analysis_and_plot_mixinputs, \
    excited_states_analysis_mixinputs

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 1
excited_states_analysis = 0
sos_analysis = 0
ppm = 0

file_msnull = '\
triplets_molecules/nf_4_4_sto3g_allmultip.out'
file_ms_notnull = '\
triplets_molecules/nf_4_4_sto3g_triplets.out'

states_option = 1  # 0: use "state_ras" ; 1: use all states_selected
states_msnull = [1, 6]  # States to be included when "states_option = 0"
states_msnotnull = [1, 3]
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ,15, 16, 17, 18, 19, 20
#####################################
#      G-TENSOR CALCULATION
#####################################
if gfactor_two_inputs == 1:
    gfactor_presentation_mixinputs(file_msnull, file_ms_notnull, states_option, states_msnull, states_msnotnull, ppm)

if excited_states_analysis == 1:
    excited_states_analysis_mixinputs(file_msnull, file_ms_notnull, states_msnull, states_option,
                                      plots=0, save_pict=0)

if sos_analysis == 1:
    sos_analysis_and_plot_mixinputs(file_msnull, file_ms_notnull, states_msnull, states_option, ppm)
