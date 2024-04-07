#####################################
#          MODULES SELECTION
#####################################
# module load Python/3.7.4-Anaconda3-2019.10
import numpy as np
# from tabulate import tabulate

from projectmethod_allparsers.parsers.parser_eom import get_eom_type, get_scf_energy, get_irreps_energies, \
    get_maximum_amplitude_orbitals, get_eom_socc_values, second_smallest_number

#####################################
#            EOM COMPARISON
#####################################
eom_input = '../Sven_EOM-MP2/vo_h2o5_2+_eomea_mp2_def2-TZVP.in.out'
threshold_excitation_energy = 10  # in eV, energy difference with respect to ground state

eom_type = get_eom_type(eom_input)
scf_energy = get_scf_energy(eom_input)

print("------------------------")
print("    EOM-CC ANALYSIS ")
print("------------------------")
print("SCF   reference energy: ", scf_energy)

# if (EOM_type == 'EOMIP'):
irreps, optimized_state_index, energies, excit_energies = get_irreps_energies(eom_input)
eom_socc = get_eom_socc_values(eom_input, optimized_state_index)
orbitals = get_maximum_amplitude_orbitals(eom_input, eom_type)
presentation_list = [['Symmetry', 'Transition', 'Energy (eV)',
                      'Excitation energy (eV)', 'SOCC (cm-1)', 'Occ. 1', 'Occ. 2', 'Virtual 1']]

if eom_type == 'EOMEA':
    irreps, optimized_state_index, energies, excit_energies = get_irreps_energies(eom_input)
    eom_socc = get_eom_socc_values(eom_input, optimized_state_index)
    orbitals = get_maximum_amplitude_orbitals(eom_input, eom_type)
    presentation_list = [['Symmetry', 'Transition', 'Energy (eV)', 'Excitation energy (eV)',
                          'SOCC (cm-1)', 'Occ. 1', 'Virtual 1', 'Virtual 2']]

eom_state = 0
transition = 0
eom_excitation_energies_ev_list = []

for i in range(0, len(irreps)):
    symmetry = irreps[i][1]
    # transition = irreps[i][0]
    total_energy = np.round(float(energies[i]), 3)
    excitation_energy = np.round(float(excit_energies[i]), 2)
    soccs = np.round(float(eom_socc[i]), 2)
    # SOCC = 0

    if excitation_energy <= threshold_excitation_energy:
        eom_excitation_energies_ev_list.append(excit_energies[i])
        transition += 1
        presentation_list.append([symmetry, transition, total_energy,
                                  excitation_energy, soccs, orbitals[i * 3], orbitals[i * 3 + 1], orbitals[i * 3 + 2]])
        eom_state += 1

    if (i < len(irreps)-1) and (irreps[i][1] != irreps[i+1][1]):
        presentation_list.append(['---'])

second_smallest_excit_energy = second_smallest_number(excit_energies)
if second_smallest_excit_energy > threshold_excitation_energy:
    print('Increase the threshold_excitation_energy')
    exit()

eom_excitation_energies_ev = np.array(eom_excitation_energies_ev_list, dtype=float)
eom_excitation_energies = eom_excitation_energies_ev / 27.211399

# print(tabulate(presentation_list, headers='firstrow'))

print('\n'.join([','.join(['{:^20}'.format(item) for item in row]) for row in presentation_list]))
