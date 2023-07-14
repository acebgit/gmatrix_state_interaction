#####################################
#          MODULES SELECTION
#####################################
# module load Python/3.7.4-Anaconda3-2019.10
import sys

from parser_init import get_eigenenergies, get_selected_states, get_number_of_states, get_spin_orbit_couplings

from parser_gtensor import from_energies_soc_to_g_values, print_g_calculation

from parser_plots import sos_analysis_and_plot

#####################################
#            INPUT VALUES
#####################################
# G-TENSOR CALCULATION
g_calculation = 1
ras_input = '\
doublets_molecules/h2o/h2o_6-31Gd_5_5.out'

selected_states = 0  # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [1, 2, 3, 4, 5]  # States to be included when "selected_states = 0"
symmetry_selection = 'A2'  # Symmetry selected states
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

sos_analysis = 0

# OUTPUT
write_ras_input = 0  # 0: write results directly; 1: write in output ras_input
output_ras_input = ras_input + '-gvalues.txt'
if write_ras_input == 1:
    sys.stdout = open(output_ras_input, "w")

###########################################
#      ACTING IN SEVERAL MOLECULES
###########################################
# Import Module
# https://www.geeksforgeeks.org/how-to-read-multiple-text-files-from-folder-in-python/
# import os
# path = "roberto_molecules"
# os.chdir(path)  # Change the directory
#
# g_list = []
# # iterate through all file
# for file in os.listdir():
#     ras_input = f"{file}"
#
#     totalstates = get_number_of_states(ras_input)
#     states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)
#     eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)
#     selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options)
#     g_shift = from_energies_soc_to_g_values(ras_input, states_ras,
#                                             totalstates, excitation_energies_ras,
#                                             selected_socs, sz_list, sz_ground)
#
#     ras_input = ras_input.replace('.out', '')
#     g_list.append([ras_input, np.round(g_shift.real[0]*1000, 3), np.round(g_shift.real[1]*1000, 3),
#           np.round(g_shift.real[2]*1000, 3)])
#
# sorted_list = sorted(g_list)
# print (tabulate(sorted_list, headers=["molecule", "gxx", "gyy", "gzz"]))
# exit()

#####################################
#      G-VALUE CALCULATION
#####################################
if g_calculation == 1:

    totalstates = get_number_of_states(ras_input)

    states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

    selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options,
                                                                 bolvin=0)

    g_shift = from_energies_soc_to_g_values(ras_input, states_ras,
                                            totalstates, excitation_energies_ras,
                                            selected_socs, sz_list, sz_ground)

    print_g_calculation(ras_input, totalstates, selected_states, symmetry_selection, states_ras, g_shift)

#####################################
#        PLOT ANALYSIS
#####################################
if sos_analysis == 1:
    sos_analysis_and_plot(ras_input)
