#########################################################
# PROGRAM TO OBTAIN THE GAP BETWEEN THE
# GROUND STATE SINGLET AND THE FIRST TRIPLET
# EXCITED  STATEs
#########################################################
import numpy as np
import os
from gaps_parsers import get_energy_and_spins, obtain_gaps, ordering_by_actspace, plots
from parser_excitstates import get_excited_states_analysis
from tabulate import tabulate

archivo = 'srdft_triplet_bia16_concat/bia16_4_4_srdft_triplet.out'
get_excited_states_analysis(archivo, cutoff=0.99)
exit()

single_file = 'srdft_triplet/bia16_2_2_srdft_triplet.out'
path = 'srdft_triplet_bia16_concat' # "srdft_triplet"
my_order = ['2_2', '4_4', '6_6', '8_8']

gaps_list = []
corrected_gaps_list = []
if path:
    os.chdir(path)  # Change the directory
    for archivos in os.listdir():
        archivo = f"{archivos}"
        with open(archivo, encoding="utf8") as file:
            content = file.readline()
            # With 'content = files.read()' does not work https://pynative.com/python-search-for-a-string-in-text-files/
            excit_energy_list, spin_list, corrected_energy_list = get_energy_and_spins(file)

            if corrected_energy_list:
                gap_t1, gap_s1 = obtain_gaps(spin_list, excit_energy_list)
                gaps_list.append([archivo, gap_s1, gap_t1, np.round(gap_s1 - gap_t1, 3)])

                gap_t1_corrected, gap_s1_corrected = obtain_gaps(spin_list, corrected_energy_list)
                corrected_gaps_list.append([archivo, gap_s1_corrected, gap_t1_corrected,
                                            np.round(gap_s1_corrected - gap_t1_corrected, 3)])
            else:
                gap_t1, gap_s1 = obtain_gaps(spin_list, excit_energy_list)
                gaps_list.append([archivo, gap_s1, gap_t1, np.round(gap_s1 - gap_t1, 3)])

    gap_list = ordering_by_actspace(gaps_list, my_order)
    corrected_gaps_list = ordering_by_actspace(corrected_gaps_list, my_order)
    gap_matrix = np.array(gap_list, dtype=object)
    corrected_gap_matrix = np.array(corrected_gaps_list, dtype=object)

    print(tabulate(gap_matrix, headers=["S1", "T1", "GAP", "S1", "T1",
                                        "GAP", "S1", "T1", "GAP", "S1", "T1", "GAP",
                                        "S1", "T1", "GAP", "S1", "T1", "GAP"]))
    print()
    plots(gap_matrix, 'Active spaces', 'Excitation energy(eV)', path)

    if corrected_energy_list:
        print(tabulate(corrected_gap_matrix, headers=["S1", "T1", "GAP", "S1", "T1", "GAP", "S1", "T1", "GAP", "S1", "T1",
                                                      "GAP", "S1", "T1", "GAP", "S1", "T1", "GAP"]))
        print()
        # plots(corrected_gap_matrix, 'Active spaces', 'Excitation energy(eV)', path)

        difference_matrix = gap_matrix
        for i in range(0, len(gap_matrix[:, 0])):
            for j in range(1, len(gap_matrix[0, :])):
                difference_matrix[i, j] = gap_matrix[i, j] - corrected_gap_matrix[i, j]

        # print(tabulate(difference_matrix, headers=["S1", "T1", "GAP", "S1", "T1", "GAP", "S1", "T1", "GAP", "S1", "T1",
        #                                            "GAP", "S1", "T1", "GAP", "S1", "T1", "GAP"]))
        plots(difference_matrix, 'Active spaces', 'Excitation energy(eV)', path)

else:
    with open(single_file, encoding="utf8") as file:
        content = file.readline()
        excit_energy_list, spin_list, corrected_energy_list = get_energy_and_spins(file)

        gap_t1, gap_s1 = obtain_gaps(spin_list, excit_energy_list)
        gaps_list.append([single_file, gap_s1, gap_t1, np.round(gap_s1 - gap_t1, 3)])
        print(gaps_list)