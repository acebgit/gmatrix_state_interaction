__author__ = 'Antonio Cebreiro-Gallardo'

######################
# Program to obtain extended information from QChem RASCI outputs
######################

#####################################
#          MODULES SELECTION
#####################################
import sys 
from projection_method.parsers.parser_gtensor import gfactor_presentation
from projection_method.parsers.parser_excitstates import get_excited_states_analysis, get_gtensor_analysis, improved_active_space
from projection_method.parsers.parser_plots import sos_analysis_and_plot, gfactor_all_states, compare_gcalculation_gestimation

#####################################
#            INPUT VALUES
#####################################
# ras_input= str(sys.argv[1])
ras_input = '../molecules/doublets/cucl4_2-_def2tzvp_11_6_d10_2.out'

######## G-TENSOR CALCULATION ########
g_calculation = 0 
ppm = 0 # 0: ppt; 1: ppm 
state_selection = 1 # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1,2,3]
symmetry_selection = 'B3u'  # Symmetry selected states_selected
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix


######## G-TENSOR EXCITED STATES ANALYSIS ########
excitanalysis_gvalue_cut = 0 # 0: not calculate; ≠0: cut-off between ground-excited states (% of maximum g-value in each dim)
gestimation = 0 # 0: g-tensor calculation (projection procedure); 1: g-tensor estimation (g = -4 L SOC / E)


######## EXCITED STATES ANALYSIS ########
excitanalysis = 0
excitanalysis_config_cut = 0.9 # cut-off for configurations amplitude (% of maximum amplitude)
excitanalysis_soc_cut = 0 # cut-off for soccs (cm-1)
excitanalysis_angmoment_cut = 0 # cut-off for orbital angular momentum (cm-1)
excit_plot = 1 # 0: not show plot, 1: show plot 


######## SOS PLOTS ########
sos_analysis = 0 # SOS g-tensor plot: g-tensor calculation with n states
gestimation_comparison = 0 # 1: SOS comparison between g-shift calculated and estimated

#  --------------------------------------------------------
gfactor_excited_states = 1


######## AUTOMATIC ACTIVE SPACE ANALYSIS ####
new_active_space = 0

########################################
#      FUNCTIONS CALLED 
########################################
# Function to calculate the g-shift 
if g_calculation == 1:
    gfactor_presentation(ras_input, states_ras, state_selection, symmetry_selection, soc_options, ppm)

if excitanalysis_gvalue_cut != 0:   
    get_gtensor_analysis(ras_input, state_selection, states_ras, symmetry_selection, excitanalysis_gvalue_cut, ppm, gestimation, cut_off=0.5)
    
if excitanalysis == 1:
    get_excited_states_analysis(ras_input, state_selection, states_ras, symmetry_selection, excitanalysis_config_cut,
                                excitanalysis_soc_cut, excitanalysis_angmoment_cut, plots=excit_plot, save_pict=0)

if sos_analysis == 1:
    sos_analysis_and_plot(ras_input, states_ras, state_selection, ppm, gestimation, order_symmetry=1, save_option=0)

if gestimation_comparison == 1:
    compare_gcalculation_gestimation(ras_input, states_ras, state_selection, ppm, plotting=1)

if gfactor_excited_states == 1:
    gfactor_all_states(ras_input, states_ras, ppm)


if new_active_space == 1:
    improved_active_space(ras_input)
