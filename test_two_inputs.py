#####################################
#          MODULES SELECTION
#####################################
from parser_mixing_inputs import gfactor_presentation, sos_analysis_and_plot
from parser_excitstates import get_excited_states_analysis

#####################################
#            INPUT VALUES
#####################################
gfactor_two_inputs = 1
excited_states_analysis = 0
sos_analysis = 0
ppm = 0

file_msnull = '\
triplets_molecules/o2_11_9_allmultip.out'
file_ms_notnull = '\
triplets_molecules/o2_11_9_triplets.out'

states_option = 0  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # States to be included when "selected_states = 0"
states_sym = 'A1'

#####################################
#      G-TENSOR CALCULATION
#####################################
if gfactor_two_inputs == 1:
    gfactor_presentation(file_msnull, file_ms_notnull, states_ras, states_option, states_sym, ppm)

if excited_states_analysis == 1:
    get_excited_states_analysis(file_ms_notnull, cutoff=0.9)

if sos_analysis == 1:
    sos_analysis_and_plot(file_msnull, file_ms_notnull, states_ras, states_option, states_sym)
