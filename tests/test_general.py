__author__ = 'Antonio Cebreiro-Gallardo'

#####################################
#          MODULES SELECTION
#####################################
from projectmethod.parsers.parser_gtensor import gfactor_presentation, from_gvalue_to_shift
from projectmethod.parsers.parser_excitstates import get_excited_states_analysis, get_gtensor_analysis
from projectmethod.parsers.parser_plots import sos_analysis_and_plot, gfactor_all_states, compare_gcalculation_gestimation
# import cProfile

#####################################
#            INPUT VALUES
#####################################
# ras_input='molecules/roberto_molecules_2/C74H32B4/C74H32B4_8_8_singletref_concat_triplets/C74H32B4_8_8_singletref_concat_triplets.out'
ras_input = '../molecules/roberto_molecules_2/C94H40B4/C94H40B4_8_8_singletref_triplets.out'

######## G-TENSOR CALCULATION ########
g_calculation = 1
ppm = 1 # 0: ppt; 1: ppm
state_selection = 1 # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20]
symmetry_selection = 'B1u'  # Symmetry selected states_selected
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix
#  [2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]

######## G-TENSOR ANALYSIS ########
excitanalysis_gvalue_cut = 0.1 # =0: not calculate; ≠0: cut-off between ground-excited states (% of maximum g-value in each dim)
gestimation = 0 # 0: g-tensor calculation (projection procedure); 1: g-tensor estimation (g = -4 L SOC / E)

######## EXCITED STATES ANALYSIS ########
excitanalysis = 1
excitanalysis_config_cut = 0.5 # cut-off for configurations amplitude (% of maximum amplitude)
excitanalysis_soc_cut = 0 # cut-off for soccs (cm-1)
excitanalysis_angmoment_cut = 0 # cut-off for orbital angular momentum (cm-1)

######## SOS PLOTS ########
sos_analysis = 0 # SOS g-tensor plot: g-tensor calculation with n states
gestimation_comparison = 0 # 1: SOS comparison between g-shift calculated and estimated

#  --------------------------------------------------------
gfactor_excited_states = 0

########################################
#      FUNCTIONS CALLED 
########################################
# cProfile.run('gfactor_presentation(ras_input, states_ras, state_selection, symmetry_selection, soc_options, ppm)')
# exit()

if g_calculation == 1:
    gfactor_presentation(ras_input, states_ras, state_selection, symmetry_selection, soc_options, ppm)

if excitanalysis_gvalue_cut != 0:   
    get_gtensor_analysis(ras_input, state_selection, states_ras, symmetry_selection, excitanalysis_gvalue_cut, ppm, gestimation, cut_off=0.5)
    
if excitanalysis == 1:
    get_excited_states_analysis(ras_input, state_selection, states_ras, symmetry_selection, excitanalysis_config_cut,
                                excitanalysis_soc_cut, excitanalysis_angmoment_cut, plots=0, save_pict=0)

if sos_analysis == 1:
    sos_analysis_and_plot(ras_input, states_ras, state_selection, ppm, gestimation, order_symmetry=1, save_option=0)

if gfactor_excited_states == 1:
    gfactor_all_states(ras_input, states_ras, ppm)

if gestimation_comparison == 1:
    compare_gcalculation_gestimation(ras_input, states_ras, state_selection, ppm, plotting=1)

###########################################
#      ACTING IN SEVERAL MOLECULES
###########################################
# several_molecules = 0
# path = "roberto_molecules"
# Import Module
# https://www.geeksforgeeks.org/how-to-read-multiple-text-files-from-folder-in-python/
# if several_molecules == 1:
#     os.chdir(path)  # Change the directory
#
#     g_list = []
#     # iterate through all file_ms_notnull
#     for file_ms_notnull in os.listdir():
#         qchem_file = f"{file_ms_notnull}"
#
#         totalstates = get_number_of_states(qchem_file)
#         states_msnull = get_selected_states(qchem_file, totalstates, states_msnull, states_option, symmetry_selection)
#         eigenenergies_ras, excitation_energies_ras = get_eigenenergies(qchem_file, totalstates, states_msnull)
#         selected_socs, list_sz, sz_ground = get_spin_orbit_couplings
#         (qchem_file, totalstates, states_msnull, soc_options)
#         g_shift = from_energies_soc_to_g_values(qchem_file, states_msnull,
#                                                 totalstates, excitation_energies_ras,
#                                                 selected_socs, list_sz, sz_ground)
#
#         qchem_file = qchem_file.replace('.out', '')
#         g_list.append([qchem_file, np.round(g_shift.real[0]*1000, 3),
#                        np.round(g_shift.real[1]*1000, 3), np.round(g_shift.real[2]*1000, 3)])
#
#     sorted_list = sorted(g_list)
#     # for i in range(0, len(sorted_list)):
#     #     print(sorted_list[i, :])
#     print(tabulate(sorted_list, headers=["molecule", "gxx", "gyy", "gzz"]))
#     exit()
