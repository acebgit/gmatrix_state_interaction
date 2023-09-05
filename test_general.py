#####################################
#          MODULES SELECTION
#####################################
import sys
import os
import numpy as np
from tabulate import tabulate

from parser_init import *  #get_eigenenergies, get_selected_states, get_number_of_states, get_spin_orbit_couplings
from parser_gtensor import *  #from_energies_soc_to_g_values, print_g_calculation
from parser_plots import *  #sos_analysis_and_plot
from parser_excitstates import *  #get_excited_states_analysis

# from_gvalue_to_shift([2.004])

#####################################
#            INPUT VALUES
#####################################
g_calculation = 1
ras_input = '\
3_David_molecules/triangulene_10_10_noholenopart.out'

several_molecules = 0
path = "2_roberto_molecules"

selected_states = 1  # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [1, 2, 3, 4, 5]  # States to be included when "selected_states = 0"
symmetry_selection = 'A2'  # Symmetry selected states
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

excited_states_analysis = 0
sos_analysis = 0

# OUTPUT
write_ras_input = 0  # 0: write results directly; 1: write in output qchem_file
output_ras_input = ras_input + '-gvalues.txt'
if write_ras_input == 1:
    sys.stdout = open(output_ras_input, "w")

###########################################
#      ACTING IN SEVERAL MOLECULES
###########################################
# Import Module
# https://www.geeksforgeeks.org/how-to-read-multiple-text-files-from-folder-in-python/
if several_molecules == 1:
    os.chdir(path)  # Change the directory

    g_list = []
    # iterate through all file
    for file in os.listdir():
        qchem_file = f"{file}"

        totalstates = get_number_of_states(qchem_file)
        states_ras = get_selected_states(qchem_file, totalstates, states_ras, selected_states, symmetry_selection)
        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(qchem_file, totalstates, states_ras)
        selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(qchem_file, totalstates, states_ras, soc_options, bolvin=0)
        g_shift = from_energies_soc_to_g_values(qchem_file, states_ras,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, sz_ground)

        qchem_file = qchem_file.replace('.out', '')
        g_list.append([qchem_file, np.round(g_shift.real[0]*1000, 3), np.round(g_shift.real[1]*1000, 3),
              np.round(g_shift.real[2]*1000, 3)])

    sorted_list = sorted(g_list)
    # for i in range(0, len(sorted_list)):
    #     print(sorted_list[i, :])
    print(tabulate(sorted_list, headers=["molecule", "gxx", "gyy", "gzz"]))
    exit()

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
if excited_states_analysis == 1:
    get_excited_states_analysis(ras_input, cutoff=0.9)

if sos_analysis == 1:
    sos_analysis_and_plot(ras_input)
