__author__ = 'Antonio Cebreiro-Gallardo'

##################################################################
# Program to obtain extended information from QChem RASCI outputs
##################################################################

#####################################
#          MODULES SELECTION
#####################################
import sys 
from projection_method.parsers.parser_gtensor import gfactor_presentation 
from projection_method.parsers.parser_excitstates import get_excited_states_analysis, gtensor_state_pairs_analysis, improved_active_space
from projection_method.parsers.parser_plots import sos_analysis_and_plot, gfactor_change_ground_state, compare_gcalculation_gestimation


#####################################
#            INPUT VALUES
#####################################
ras_input= str(sys.argv[1])
# ras_input = '../../molecules/doublets/cucl4_2-_def2tzvp_11_6_d10_2.out'


######## G-TENSOR CALCULATION ########
calculate_gshift = 1 
ppm = 1 # 0: ppt; 1: ppm 
state_selection = 1 # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry

states_ras = list(range(1, 526)) 
# states_ras.extend([148, 149, 151, 152, 153, 154, 155]) 
# states_ras.remove(2) 
# states_ras.insert(0, 2) 

symmetry_selection = 'B2'  # Symmetry selected states_selected
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix


######## G-TENSOR CALCULATION BY PAIRS ########
gshift_estimation_by_state_pairs = 0.5
cut_off_config = 0.75
# ≠0: cut-off between ground-excited states (% of maximum g-value in each dim), where g = -4 L SOC / E


######## SOS PLOTS ########
sum_over_states_analysis = 0 # SOS g-tensor plot: g-tensor calculation with n states
gestimation_gsos_comparison = 0 # 1: SOS comparison between g-shift calculated and estimated


######## EXCITED STATES ANALYSIS ########
excited_states_analysis = 0
excitanalysis_config_cut = 0.5 # cut-off for configurations amplitude (% of maximum amplitude)
excitanalysis_soc_cut = 0.05 # cut-off for soccs (% of maximum SOCC)
excitanalysis_angmoment_cut = 10E-10 # cut-off for orbital angular momentum (% of maximum L)
excit_plot = 0 # 0: not show plot, 1: show plot 


######## G-TENSOR CALCULATION USING DIFFERENT STATES AS GROUND STATE ########
gshift_changing_ground_state = 0


#################################
#      FUNCTIONS CALLED 
#################################
if calculate_gshift == 1:
    gfactor_presentation(ras_input, states_ras, state_selection, symmetry_selection, soc_options, ppm)

if gshift_estimation_by_state_pairs != 0:
    gtensor_state_pairs_analysis(ras_input, state_selection, states_ras, symmetry_selection, gshift_estimation_by_state_pairs, ppm, cut_off_config)
elif gshift_estimation_by_state_pairs == 0 and excited_states_analysis == 1:
    get_excited_states_analysis(ras_input, state_selection, states_ras, symmetry_selection, excitanalysis_config_cut,
                                excitanalysis_soc_cut, excitanalysis_angmoment_cut, plots=excit_plot, save_pict=0)

if sum_over_states_analysis == 1:
    sos_analysis_and_plot(ras_input, states_ras, state_selection, ppm, order_symmetry=0, save_option=0)

if gestimation_gsos_comparison == 1:
    compare_gcalculation_gestimation(ras_input, states_ras, state_selection, ppm, plotting=1)

if gshift_changing_ground_state == 1:
    gfactor_change_ground_state(ras_input, states_ras, state_selection, symmetry_selection, soc_options, ppm)
