#####################################
#          MODULES SELECTION
#####################################
__author__ = 'Antonio Cebreiro-Gallardo'

import sys

from parser_gtensor import gfactor_presentation
from parser_excitstates import get_excited_states_analysis, improved_active_space
from parser_plots import sos_analysis_and_plot, gfactor_all_states

# from_gvalue_to_shift([2.00292, 2.00323, 2.00246])
# exit()
#####################################
#            INPUT VALUES
#####################################
ras_input = '\
different_multiplicities/h2o_doublet.out'

g_calculation = 1
sos_analysis = 0
gfactor_excited_states = 0
ppm = 0

excited_states_analysis = 0
improve_as = 0

state_selection = 1 # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1, 2, 10]  # States to be included when "states_option = 0"
# [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
symmetry_selection = 'A2'  # Symmetry selected states_selected
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

# OUTPUT
write_ras_input = 0  # 0: write results directly; 1: write in output qchem_file
output_ras_input = ras_input + '-gvalues.txt'
if write_ras_input == 1:
    sys.stdout = open(output_ras_input, "w")

#      G-VALUE CALCULATION
if g_calculation == 1:
    gfactor_presentation(ras_input, states_ras, state_selection, symmetry_selection, soc_options, ppm)

if excited_states_analysis == 1:
    get_excited_states_analysis(ras_input, state_selection, states_ras, cut_off=0.9, plots=0, save_pict=0)

if improve_as == 1:
    improved_active_space(ras_input, cut_off=0.9, see_soc=1)

if sos_analysis == 1:
    sos_analysis_and_plot(ras_input, states_ras, state_selection, ppm, order_symmetry=1, save_option=0)

if gfactor_excited_states == 1:
    gfactor_all_states(ras_input, states_ras, ppm)

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
#         selected_socs, list_sz, sz_ground = get_spin_orbit_couplings(qchem_file, totalstates, states_msnull, soc_options)
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
