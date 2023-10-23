"""
Obtain an improved active space of RAS-CI Q-Chem output including orbitals
with unpaired electrons in relevant hole/particle configurations.
"""
import numpy as np
import sys

from parser_gtensor import get_number_of_states, get_symmetry_states, get_selected_states, get_eigenenergies
from parser_excitstates import get_hole_part_contributions, get_ras_spaces, get_alpha_beta, \
    get_highest_amplitudes, get_orbital



# def improved_active_space(file):
#     """
#     Obtain an improved active space of RAS-CI Q-Chem output including orbitals with unpaired electrons in relevant
#     hole/particle configurations.
#     :param file:
#     :return:
#     """
#     totalstates = get_number_of_states(file)
#
#     state_symmetries, ordered_state_symmetries = get_symmetry_states(file, totalstates)
#
#     states_ras = get_selected_states(file, totalstates, selected_states='None', states_option=1, symmetry_selection='None')
#
#     eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
#
#     hole_contributions, part_contributions = get_hole_part_contributions(file, totalstates)
#
#     orbital_momentum = get_ground_state_orbital_momentum(file, totalstates)
#
#     socc_values = get_socc_values(file, totalstates)
#
#     mulliken_charge, mulliken_spin = get_mulliken_spin(file, totalstates, states_ras)
#
#     initial_active_orbitals, ras_occ = get_ras_spaces(file)
#
#     alpha, homo_orbital = get_alpha_beta(file)
#
#     # TAKE STATES WITH HIGHEST CONTRIBUTION
#     new_space_list = []
#     final_scf_ordered_space = []
#     excited_states_presentation_list = ['State', 'Symmetry', 'Hole', 'Part',
#                                         'Excitation energy (eV)', 'Orbitals', 'SOCC (cm-1)',
#                                         'Orbital momentum']
#
#     word_search = ' | HOLE  | '
#     n_states = 0
#
#     with open(file, encoding="utf-8") as file:
#         for line in file:
#             if word_search in line:  # Go to configurations line
#
#                 index_max_amplitudes, state_orbitals = get_highest_amplitudes(file, cutoff)
#
#                 for i in index_max_amplitudes:
#                     new_orbital = get_orbital(homo_orbital, state_orbitals[i], initial_active_orbitals)
#
#                     if new_orbital not in new_space_list:
#                         new_space_list.append(new_orbital)
#                         excited_states_presentation_list, soc = print_excited_states(excited_states_presentation_list,
#                                                                                      n_states, hole_contributions,
#                                                                                      part_contributions, socc_values,
#                                                                                      excitation_energies_ras
#                                                                                      * 27.211399,
#                                                                                      state_symmetries, new_orbital,
#                                                                                      orbital_momentum, mulliken_spin)
#
#                         if soc != 0 or i == 0:
#                             final_scf_ordered_space.append(new_orbital)
#
#                 n_states += 1
#
#     print("------------------------")
#     print(" IMPROVED ACTIVE SPACE ")
#     print("------------------------")
#
#     print('Initial active space (HOMO =', homo_orbital, '):')
#     initial_active_orbitals_list = initial_active_orbitals.tolist()
#     electrons = get_new_active_space_electrons(initial_active_orbitals_list, homo_orbital)
#     print('[', electrons, ',', len(initial_active_orbitals_list), '] ;', initial_active_orbitals_list)
#
#     # print('')
#     # new_active_space = np.array(new_space_list, dtype=int)
#     # print('New active space (SOCC not zero, HOMO singly occupied):')
#     # electrons = get_new_active_space_electrons(final_SCF_ordered_space, homo_orbital)
#     # final_SCF_ordered_space.sort()
#     # print('[',electrons, ',', len(final_SCF_ordered_space), '] ;', final_SCF_ordered_space)
#
#     print('')
#     print('Initial active space + New active space:')
#
#     initial = set(initial_active_orbitals)
#     final = set(final_scf_ordered_space)
#     in_final_but_not_in_initial = final - initial
#     initial_plus_final_space = initial_active_orbitals + list(in_final_but_not_in_initial)
#
#     electrons = get_new_active_space_electrons(initial_plus_final_space, homo_orbital)
#     initial_plus_final_space.sort()
#     print('[', electrons, ',', len(initial_plus_final_space), '] ;', initial_plus_final_space)
