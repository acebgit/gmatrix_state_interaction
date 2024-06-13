"""
Calculation of the g-tensor using Q-Chem output with RAS-CI
"""
__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np
import matplotlib.pyplot as plt

from projection_method.parsers.parser_gtensor import get_number_of_states, get_symmetry_states, get_selected_states, get_eigenenergies, \
    get_spin_orbit_couplings, get_spin_matrices, get_orbital_matrices, hermitian_test, \
    get_hamiltonian_construction, diagonalization, angular_matrices_obtention, g_factor_calculation
from projection_method.parsers.parser_excitstates import s2_from_file, get_hole_part_contributions, get_groundst_socc_values, \
    get_groundst_orbital_momentum, get_bar_chart

def mapping_between_states(file_msnull, file_msnotnull, states_msnull, states_msnotnull, totalstates_msnull, totalstates_msnotnull):
    """
    Comparing the energies, map states_selected with Ms = 0, that do not have coupling between states_selected triplets,
    and states_selected with Ms not 0, that do have this coupling.
    :return: mapping_list
    """
    def check_mapping(energ_ms_notnull, energ_ms_null, mapping_dic, list):
        """
        Check if states "Ms=0" and "Ms ≠ 0" have the same energy. If so, to include them in the list, check that states
        have not been already included in the list, meaning checking for degeneracies.
        :return: list
        """
        if (energ_ms_notnull == energ_ms_null):
            list.append(mapping_dic)

            for k in range(0, len(list) - 1):
                if (mapping_dic['state ms not null'] == list[k]['state ms not null']) or \
                        (mapping_dic['state ms null'] == list[k]['state ms null']):
                    list.remove(mapping_dic)
        return list

    ms_notnull_energies, ms_notnull_excit_energ = get_eigenenergies(file_msnotnull, totalstates_msnotnull, states_msnotnull)
    ms_null_energies, ms_null_excit_energ = get_eigenenergies(file_msnull, totalstates_msnull, states_msnull)

    mapping_list = []
    mapping_dict = {}

    for i in range(0, len(ms_notnull_energies)):
        for j in range(0, len(ms_null_energies)):
            ener_ms_notnull = np.round(ms_notnull_energies[i], 5)
            ener_ms_null = np.round(ms_null_energies[j], 5)

            mapping_dict = {'state ms not null': i, 'state ms null': j}
            mapping_list = check_mapping(ener_ms_notnull, ener_ms_null, mapping_dict, mapping_list)

    print('Mapping: state Ms not 0 - Ms 0')
    for mapping_dict in mapping_list:
        a = mapping_dict['state ms not null']
        b = mapping_dict['state ms null']
        print(states_msnotnull[a], ' - ', states_msnull[b])
    print('---')
    exit()
    if mapping_list == []:
        raise ValueError("No mapping between states is possible: the states of the two inputs are different.")
    return mapping_dict, mapping_list


def matrix_expansion(sz_small, sz_big, nstates, socs, orbit_moments, spin_moments):
        """
        Expanded the small matrixes with multiplicities "sz_small" to a big matrix with multiplicities
        "sz_big"
        :return: matrix_big
        """
        def general_matrix_expansion(nstates, matrix_small, sz_small, sz_big):
            """
            Expanded the small matrix with multiplicities "sz_small" to a big matrix with multiplicities
            "sz_big"
            :return: matrix_big
            """
            matrix_big = 0

            if (matrix_small.ndim == 2):
                matrix_big = np.zeros((len(nstates) * len(sz_big), len(nstates) * len(sz_big)), dtype=complex)
                for i in range(0, len(nstates)):
                    for j in range(0, len(nstates)):
                        i_position = i * len(sz_small)
                        j_position = j * len(sz_small)

                        i_reordered = i * len(sz_big) + (len(sz_big) - len(sz_small)) // 2
                        j_reordered = j * len(sz_big) + (len(sz_big) - len(sz_small)) // 2

                        for sz_1 in range(0, len(sz_small)):
                            for sz_2 in range(0, len(sz_small)):
                                matrix_big[i_reordered + sz_1, j_reordered + sz_2] = \
                                    matrix_small[i_position + sz_1, j_position + sz_2]

            elif (matrix_small.ndim == 3):
                matrix_big = np.zeros((len(nstates) * len(sz_big), len(nstates) * len(sz_big), 3), dtype=complex)
                for k in range(0, 3):
                    for i in range(0, len(nstates)):
                        for j in range(0, len(nstates)):
                            i_position = i * len(sz_small)
                            j_position = j * len(sz_small)

                            i_reordered = i * len(sz_big) + (len(sz_big) - len(sz_small)) // 2
                            j_reordered = j * len(sz_big) + (len(sz_big) - len(sz_small)) // 2

                            for sz_1 in range(0, len(sz_small)):
                                for sz_2 in range(0, len(sz_small)):
                                    matrix_big[i_reordered + sz_1, j_reordered + sz_2, k] = \
                                        matrix_small[i_position + sz_1, j_position + sz_2, k]
            return matrix_big

        socs = general_matrix_expansion(nstates, socs, sz_small, sz_big)
        orbit_moments = general_matrix_expansion(nstates, orbit_moments, sz_small, sz_big)
        spin_moments = general_matrix_expansion(nstates, spin_moments, sz_small, sz_big)
        return sz_big, socs, orbit_moments, spin_moments


def get_input_values(ras_input, states_ras, selected_states, symmetry_selection, soc_options):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, states_msnull, states_option, symmetry_selection, soc_options
    :return: g-shifts
    """
    totalstates = get_number_of_states(ras_input)

    states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

    selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options)

    spin_matrix, standard_spin_matrix = get_spin_matrices(ras_input, states_ras)

    orbital_matrix = get_orbital_matrices(ras_input, totalstates, states_ras, sz_list)

    return totalstates, states_ras, eigenenergies_ras, selected_socs, sz_list, sz_ground, spin_matrix, \
           standard_spin_matrix, orbital_matrix


def totalstates_mix(states_1, states_2, list_mapping):
        """
        Total number of states is the sum of both total number of states substrating those that are in common
        :param: totalstates_msnull, totalstates_msnotnull, list_mapping
        :return: totalstates
        """
        totalstates = len(states_1) + len(states_2) - len(list_mapping)
        # print('totalstates: ', totalstates)
        # print('---')
        return totalstates


def mixing_lists(list_1, list_2, mapping_list):
        """
        Give the energies of the states in list_1 and list_2 without those states that are
        repeat in both lists (meaning those that are in the mapping list)
        :param: list_1, list_2, mapping_list
        :return: final_list
        """
        elements = list(list_1)

        repeat_states = []
        for mapping_dict in mapping_list:  # States that are in list_1 and list_2
            repeat_states.append(mapping_dict['state ms not null'])

        for state in range(0, len(list_2)):  # States that are in list_2
                if state not in repeat_states:
                    elements.append(list_2[state])

        final_list = np.array(elements)
        return final_list


def include_msnotnull_states(states_msnull, states_msnotnull, mapping_list, soc_msnotnull, total_socs, list_sz):
    """
    SOCs between states Ms ≠ 0 that have not been included in total SOC matrix are now included
    :return:
    """
    def get_mapped_states(list_mapping, states):
        """
        Obtain a list of not mapped states
        :return: no_mapped_states
        """
        mapped_states_index = []
        for mapping_dict in list_mapping:
            mapped_states_index.append(mapping_dict['state ms not null'])

        no_mapped_states = []
        mapped_states = []
        for i in range(0, len(states)):
            if i not in mapped_states_index:
                no_mapped_states.append(i)
            # else:
            #     mapped_states.append(i)
        return no_mapped_states

    def state_arrangement_in_matrix(i, j, states_msnull, no_mapped_states, mapping_list):
        """
        Set the final SOC matrix row or column defined by the state selected.
        :param
        :return: state_i_final, state_j_final
        """
        # Upper left part of the matrix: Begins the writing after the SOCs matrix between states with Ms = 0
        state_i_final = len(states_msnull)
        state_j_final = len(states_msnull)

        if (i in no_mapped_states) and (j in no_mapped_states):
            # Lower right part of the matrix: when both states are not in mapping
            # print('LR')  # "Lower Right"
            state_i_final += no_mapped_states.index(i)
            state_j_final += no_mapped_states.index(j)

        elif (i in no_mapped_states) and (j not in no_mapped_states):
            # Lower left part of the matrix: when rows are in matrix but not the columns
            # print('LL')  # "Lower Left"
            for k in range(0, len(mapping_list)):
                if j == mapping_list[k]['state ms not null']:
                    index = k
                    break
            state_i_final += no_mapped_states.index(i)
            state_j_final = mapping_list[index]['state ms null']

        elif (i not in no_mapped_states) and (j in no_mapped_states):
            # Upper right part of the matrix: when columns are in matrix but not the columns rows
            # print('UR')  # "Upper Right"
            for k in range(0, len(mapping_list)):
                if i == mapping_list[k]['state ms not null']:
                    index = k
                    break
            state_i_final = mapping_list[index]['state ms null']
            state_j_final += no_mapped_states.index(j)

        return state_i_final, state_j_final

    # Obtain list of not mapped states
    no_mapped_states = get_mapped_states(mapping_list, states_msnotnull)

    # Substitution
    for i in range(0, len(states_msnotnull)):
        for j in range(0, len(states_msnotnull)):

            if (i in no_mapped_states) or (j in no_mapped_states):
                state_i_final, state_j_final = state_arrangement_in_matrix(i, j, states_msnull, no_mapped_states, mapping_list)

                # if (state_i_final in mapped_states) or (state_j_final in mapped_states):
                for sz_1 in range(0, len(list_sz)):
                    for sz_2 in range(0, len(list_sz)):
                        soc_row = i * len(list_sz) + sz_1
                        soc_col = j * len(list_sz) + sz_2

                        soc_row_final = state_i_final * len(list_sz) + sz_1
                        soc_col_final = state_j_final * len(list_sz) + sz_2
                        # print(soc_row, soc_col, '-->', soc_row_final, soc_col_final)

                        total_socs[soc_row_final, soc_col_final] = soc_msnotnull[soc_row, soc_col]
                # print('SOC in Ms not null:', soc_msnotnull[i * len(list_sz), j * len(list_sz)],
                #       'SOC in total:', total_socs[state_i_final * len(list_sz), state_j_final * len(list_sz)])
                # exit()
    return total_socs


def socs_mix(states_msnull, states_msnotnull, socs_msnull, socs_msnotnull, mapping_list, sz_list, totalstates):
        """
        Give the SOCS of the states in msnull_ang and soc_msnotnull without those states that are
        repeat in both lists (meaning those that are in the mapping list)
        :param: msnull_ang, msnotnull_ang, mapping_list, list_sz, totalstates
        :return: total_soc
        """
        def exchanging_socs(total_soc, soc_msnotnull, list_mapping, list_sz):
            """
            Put SOCs between states_selected with Ms different than 0 (that are obtained in the output) in the
            SOC matrix of states_selected with Ms 0 (that are not obtained since Clebsh-Gordan coefficient is too small)
            :param: total_socs, soc_msnotnull, list_mapping, list_sz
            :return: total_socs
            """
            # print('SOC of Ms null (not include SOCs between triplets):')
            # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
            #                  for row in np.round((total_socs[:, :] * 219474.63068), 5)]))  # * 219474.63068
            # print()
            for i in list_mapping:
                for j in list_mapping:
                    if i['state ms not null'] != j['state ms not null']:  # If states_selected are not the same (in Ms not null list)

                        for sz_1 in range(0, len(list_sz)):
                            for sz_2 in range(0, len(list_sz)):
                                i_msnotnull = i['state ms not null'] * len(list_sz) + sz_1
                                j_msnotnull = j['state ms not null'] * len(list_sz) + sz_2

                                i_msnull = i['state ms null'] * len(list_sz) + sz_1
                                j_msnull = j['state ms null'] * len(list_sz) + sz_2

                                total_soc[i_msnull, j_msnull] = soc_msnotnull[i_msnotnull, j_msnotnull]

            # print('SOC of Ms null (include SOCs between triplets):')
            # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
            #                  for row in np.round((total_soc[:, :] * 219474.63068), 5)]))  # * 219474.63068
            # print()
            # exit()
            return total_soc

        # 1) Upper left part of the matrix: First block of SOC total matrix is the SOCs of Ms = 0
        total_soc = np.zeros((totalstates * len(sz_list), totalstates * len(sz_list)), dtype=complex)
        for i in range(0, len(socs_msnull)):
            for j in range(0, len(socs_msnull)):
                total_soc[i, j] = socs_msnull[i, j]

        # 2) SOCs in total matrix between states mapped (those in which SOC is not calculated because Clebsch-Gordan coefficient
        # not calculated) are exchanged by those calculated couplings in Ms ≠ 0.
        total_soc = exchanging_socs(total_soc, socs_msnotnull, mapping_list, sz_list)

        # 3) Those SOCs between states Ms ≠ 0 that have not been included in total SOC matrix are now included
        total_soc = include_msnotnull_states(states_msnull, states_msnotnull, mapping_list, socs_msnotnull, total_soc, sz_list)

        # print('msnull_ang:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                 for row in np.round((msnull_ang[:,:]* 219474.63068),5)])) # * 219474.63068
        # print('---')
        # print('msnotnull_ang:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                 for row in np.round((msnotnull_ang[:,:]* 219474.63068),5)])) # * 219474.63068
        # print('---')
        # print('Total SOC:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                 for row in np.round((total_soc[:,:]* 219474.63068),5)])) # * 219474.63068
        # exit()
        return total_soc


def angular_momentums_mix(states_msnull, states_msnotnull, msnull_ang, msnotnull_ang, mapping_list, sz_list, totalstates):
    """
    Give the angular momentum of the states in msnull_ang and soc_msnotnull without those states that are
    repeat in both lists (meaning those that are in the mapping list)
    :param:
    :return:
    """
    # 1) First block of momentum total matrix is the angular moment of Ms = 0
    total_ang = np.zeros((totalstates * len(sz_list), totalstates * len(sz_list), 3), dtype=complex)
    for k in range(0, 3):
        for i in range(0, len(msnull_ang)):
            for j in range(0, len(msnull_ang)):
                total_ang[i, j, k] = msnull_ang[i, j, k]

    # 2) Those momentums between states Ms ≠ 0 that have not been included in total momentum matrix are now included
    for k in range(0, 3):
        total_ang[:,:,k] = include_msnotnull_states(states_msnull, states_msnotnull, mapping_list, msnotnull_ang[:,:,k], total_ang[:,:,k], sz_list)
        hermitian_test(total_ang[:,:,k], sz_list)

    # for k in range(0, 3):
    #     print('msnull_ang:')
    #     print('Dimension: ', k)
    #     print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
    #                      for row in np.round((msnull_ang[:, :, k]), 5)]))
    #     print(" ")
    #     print('---')
    #
    #     print('msnotnull_ang:')
    #     print('Dimension: ', k)
    #     print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
    #                      for row in np.round((msnotnull_ang[:, :, k]), 5)]))
    #     print(" ")
    #     print('---')
    #
    #     print('total_ang:')
    #     print('Dimension: ', k)
    #     print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
    #                      for row in np.round((total_ang[:, :, k]), 5)]))
    #     print(" ")
    #     print('---')
    # exit()
    return total_ang

def count_state_multiplicities(file_msnull, file_ms_notnull, states_ras_msnull, states_ras_msnotnull, list_mapping):
    """
    Count number of triplet and singlet states included.
    :param file_msnull, file_ms_notnull, states_msnull, list_mapping:
    :return: triplets, singlets
    """
    s2_list_1 = s2_from_file(file_msnull, states_ras_msnull)
    s2_list_2 = s2_from_file(file_ms_notnull, states_ras_msnotnull)
    s2_matrix = mixing_lists(s2_list_1, s2_list_2, list_mapping)

    singlets = 0
    triplets = 0
    for i in range(0, len(s2_matrix)):
        if 0 == s2_matrix[i]:
            singlets += 1
        if 2 == s2_matrix[i]:
            triplets += 1
    return singlets, triplets


def print_g_calculation_mixinputs(file, totalstates, states_ras_1,
                        states_ras_2, upper_g_tensor_results_ras, singlets, triplets):

    print("--------------------------------------")
    print("     INPUT SECTION")
    print("--------------------------------------")
    print("File selected: ", file)
    print("Number of states: ", totalstates, " (" , singlets, "singlets and", triplets, "triplets)")
    print("Selected states selected in Ms = 0: ", states_ras_1)
    print("Selected states selected in Ms ≠ 0: ", states_ras_2)

    print(" ")
    print("------------------------")
    print(" RAS-CI RESULTS")
    print("------------------------")
    print('g-factor (x y z dimensions):')
    print(np.round(upper_g_tensor_results_ras.real[0], 3), np.round(upper_g_tensor_results_ras.real[1], 3),
          np.round(upper_g_tensor_results_ras.real[2], 3))
    print('')


def gfactor_presentation_mixinputs(file_msnull, file_ms_notnull, selection_states, states_msnull, states_msnotnull, ppms):
    """
    From both files with Ms = 0 and Ms ≠0 state, obtain the presentation list with the g-values obtained.
    """
    # Obtain all the data from both inputs
    totalstates_1, states_ras_1, eigenenergies_ras_1, selected_socs_1, sz_list_1, sz_ground_1, \
    spin_matrix_1, standard_spin_matrix_1, orbital_matrix_1 = \
        get_input_values(file_msnull, states_msnull, selection_states, symmetry_selection='None', soc_options=0)

    totalstates_2, states_ras_2, eigenenergies_ras_2, selected_socs_2, sz_list_2, sz_ground_2, \
    spin_matrix_2, standard_spin_matrix_2, orbital_matrix_2 = \
        get_input_values(file_ms_notnull, states_msnotnull, selection_states, symmetry_selection='None', soc_options=0)

    # In case both maximum multiplicities differ, expand matrices of the shortest to the largest multiplicity
    sz_list = sz_list_1
    sz_ground = sz_ground_1
    standard_spin_matrix = standard_spin_matrix_1

    if len(sz_list_1) < len(sz_list_2):
        sz_ground = sz_ground_2
        standard_spin_matrix = standard_spin_matrix_2
        sz_list, selected_socs_1, orbital_matrix_1, spin_matrix_1 = \
            matrix_expansion(sz_list_1, sz_list_2, states_ras_1, selected_socs_1, orbital_matrix_1, spin_matrix_1)
    elif len(sz_list_1) > len(sz_list_2):
        sz_ground = sz_ground_1
        standard_spin_matrix = standard_spin_matrix_1
        sz_list, selected_socs_2, orbital_matrix_2, spin_matrix_2 = \
            matrix_expansion(sz_list_2, sz_list_1, states_ras_2, selected_socs_2, orbital_matrix_2, spin_matrix_2)

    # Make a list with the mapping between states of Ms = 0 and Ms ≠0
    dict_mapping, list_mapping = mapping_between_states(file_msnull, file_ms_notnull, states_ras_1, states_ras_2, totalstates_1,totalstates_2)

    totalstates = totalstates_mix(states_ras_1, states_ras_2, list_mapping)

    eigenenergy = mixing_lists(eigenenergies_ras_1, eigenenergies_ras_2, list_mapping)

    socs = socs_mix(states_ras_1, states_ras_2 , selected_socs_1, selected_socs_2, list_mapping, sz_list, totalstates)

    hamiltonian = get_hamiltonian_construction(eigenenergy, socs, sz_list)
    # print('Hamiltonian:')
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                 for row in np.round((hamiltonian[:,:]),5)* 219474.63068]))  # * 219474.63068
    # print()
    # exit()

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

    orbital_matrix = angular_momentums_mix(states_msnull, states_msnotnull, orbital_matrix_1, orbital_matrix_2, list_mapping, sz_list, totalstates)

    spin_matrix = angular_momentums_mix(states_msnull, states_msnotnull, spin_matrix_1, spin_matrix_2, list_mapping, sz_list, totalstates)

    combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, sz_list)

    combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, sz_list)

    g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                                   sz_list, sz_ground, ppms)

    singlet, triplet = count_state_multiplicities(file_msnull, file_ms_notnull, states_ras_1, states_ras_2, list_mapping)

    print_g_calculation_mixinputs(file_msnull, totalstates, states_ras_1, states_ras_2, g_shift, singlet, triplet)


def get_input_data_excited_states(file, state_selections, states_ras):
    """
    Obtaining a matrix with several data for each excited state. The cut-off determines the fraction of the amplitude
    of the 1st configuration that need to have the other configurations to be shown in each state.
    :param: file_ms_notnull, cutoff
    :return: excited_states_presentation_matrix
    """
    totalstates = get_number_of_states(file)

    states_ras = get_selected_states(file, totalstates, states_ras, state_selections, symmetry_selection='None')

    state_symmetries_all, ordered_state_symmetries_all = get_symmetry_states(file, totalstates)

    state_symmetries = []
    ordered_state_symmetries = []
    for i in states_ras:
        state_symmetries.append(state_symmetries_all[i-1])
        ordered_state_symmetries.append(ordered_state_symmetries_all[i-1])

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
    excitation_energies_ras[:] = (excitation_energies_ras[:] - excitation_energies_ras[0]) * 27.211399

    s2_list = s2_from_file(file, states_ras)

    hole_contributions, part_contributions = get_hole_part_contributions(file, totalstates, states_ras)

    socc_values = get_groundst_socc_values(file, totalstates, states_ras)

    orbital_momentum = get_groundst_orbital_momentum(file, totalstates, states_ras)

    return states_ras, totalstates, state_symmetries, eigenenergies_ras, s2_list, hole_contributions, part_contributions, \
        socc_values, orbital_momentum


def excited_states_analysis_mixinputs(file_msnull, file_ms_notnull, states_ras, states_option, plots, save_pict):
    """
    Obtaining a matrix with several data for each excited state. The cut-off determines the fraction of the amplitude
    of the 1st configuration that need to have the other configurations to be shown in each state.
    :param: file_ms_notnull, cutoff
    :return: excited_states_presentation_matrix
    """
    # Data obtaining from both inputs
    states_ras_1, totalstates_1, ordered_state_symmetries_1, eigenenergies_ras_1, s2_list_1, hole_contributions_1, part_contributions_1, \
        socc_values_1, orbital_momentum_1 =  get_input_data_excited_states(file_msnull, states_option, states_ras)

    states_ras_2, totalstates_2, ordered_state_symmetries_2, eigenenergies_ras_2, s2_list_2, hole_contributions_2, part_contributions_2, \
        socc_values_2, orbital_momentum_2 =  get_input_data_excited_states(file_ms_notnull, states_option, states_ras)

    # Mapping between states
    dict_mapping, list_mapping = mapping_between_states(file_msnull, file_ms_notnull, states_ras_1, states_ras_2, totalstates_1, totalstates_2)

    # Mapping data between states
    totalstates = totalstates_mix(states_ras_1, states_ras_2, list_mapping)

    ordered_state_symmetries = mixing_lists(ordered_state_symmetries_1, ordered_state_symmetries_2, list_mapping)

    eigenenergy = mixing_lists(eigenenergies_ras_1, eigenenergies_ras_2, list_mapping)
    eigenenergy[:] = (eigenenergy[:] - eigenenergy[0]) * 27.211399

    s2_list = mixing_lists(s2_list_1, s2_list_2, list_mapping)

    hole_contributions = mixing_lists(hole_contributions_1, hole_contributions_2, list_mapping)

    part_contributions = mixing_lists(part_contributions_1, part_contributions_2, list_mapping)

    socc_values = mixing_lists(socc_values_1, socc_values_2, list_mapping)

    orbital_momentum = mixing_lists(orbital_momentum_1, orbital_momentum_2, list_mapping)

    file_string = file_msnull

    excited_states_presentation_list = [['State', 'Symmetry', 'Hole', 'Part',
                                         'Excitation energy(eV)', 'SOCC(cm-1)',
                                         'Orbital momentum(máx)', 'S^2']]

    for state_index in range(0, totalstates):
        symmetry = ordered_state_symmetries[state_index]

        hole = np.around(float(hole_contributions[state_index]), 2)
        part = np.around(float(part_contributions[state_index]), 2)

        excit_energy = np.round(float(eigenenergy[state_index]), 3)
        soc = np.round(float(socc_values[state_index]), 0)

        orbital_ground_state = np.round(float(orbital_momentum[state_index]), 3)
        s2 = s2_list[state_index]

        excited_states_presentation_list.append([state_index+1, symmetry, hole, part, excit_energy,
                                                 soc, orbital_ground_state, s2])

    excited_states_presentation_matrix = np.array(excited_states_presentation_list, dtype=object)

    print("------------------------")
    print(" EXCITED STATES ANALYSIS")
    print("------------------------")
    print('Most important settings for each state (amplitude_cutoff: 0.3) :')
    print('\n'.join(''.join('{:^20}'.format(item) for item in row)
                    for row in excited_states_presentation_matrix[:, :]))
    print(" ")

    states = list(range(1, totalstates+1))
    if plots == 1:
        get_bar_chart(file_string[:-4], states, eigenenergy, 'Electronic State',
                      'Excitation energy (eV)', 'energ_analysis', save_pict)
        get_bar_chart(file_string[:-4], states, orbital_momentum, 'Electronic State',
                      'Orbital angular momentum', 'orbitmoment_analysis', save_pict)
        get_bar_chart(file_string[:-4], states, socc_values, 'Electronic State',
                      'SOCC (cm-1)', 'socc_analysis', save_pict)


def gfactor_sos_analysis(file_msnull, file_ms_notnull, states_ras, states_option):
    """
    From both files with Ms = 0 and Ms ≠0 state, obtain the presentation list with the g-values obtained.
    :return: analysis
    """
    dict_mapping, list_mapping = mapping_between_states(file_msnull, file_ms_notnull, states_ras, states_option)

    totalstates_1, states_ras_1, eigenenergies_ras_1, selected_socs_1, sz_list_1, sz_ground_1, \
    spin_matrix_1, standard_spin_matrix_1, orbital_matrix_1 = \
        get_input_values(file_msnull, states_ras, states_option, soc_options=0)

    totalstates_2, states_ras_2, eigenenergies_ras_2, selected_socs_2, sz_list_2, sz_ground_2, \
    spin_matrix_2, standard_spin_matrix_2, orbital_matrix_2 = \
        get_input_values(file_ms_notnull, states_ras, states_option, soc_options=0)

    # In case both maximum multiplicities differ, expand one of the matrices to the multiplicity of the other one
    if len(sz_list_1) < len(sz_list_2):
        sz_list, selected_socs_1, orbital_matrix_1, spin_matrix_1 = \
            matrix_expansion(sz_list_1, sz_list_2, states_ras_1, selected_socs_1, orbital_matrix_1, spin_matrix_1)
    elif len(sz_list_1) > len(sz_list_2):
        sz_list, selected_socs_2, orbital_matrix_2, spin_matrix_2 = \
            matrix_expansion(sz_list_2, sz_list_1, states_ras_2, selected_socs_2, orbital_matrix_2, spin_matrix_2)
    elif len(sz_list_1) == len(sz_list_2):
        sz_list = sz_list_1

    totalstates = totalstates_mix(states_ras_1, states_ras_1, list_mapping)

    eigenenergy = mixing_lists(eigenenergies_ras_1, eigenenergies_ras_2, list_mapping)

    socs = socs_mix(selected_socs_1, selected_socs_2, list_mapping, sz_list_1, totalstates)

    hamiltonian = get_hamiltonian_construction(eigenenergy, socs, sz_list_1)

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

    spin_matrix = angular_momentums_mix(spin_matrix_1, spin_matrix_2, list_mapping, sz_list_1, totalstates)

    orbital_matrix = angular_momentums_mix(orbital_matrix_1, orbital_matrix_2, list_mapping, sz_list_1, totalstates)

    combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, sz_list_1)

    combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, sz_list_1)

    g_shift = g_factor_calculation(standard_spin_matrix_1, combination_spin_matrix, combination_orbital_matrix,
                                   sz_list_1, sz_ground_1)
    return g_shift


def plot_g_tensor_vs_states(presentation_matrix, x_title, y_title, main_title, save_picture):
    fig, ax = plt.subplots()

    # MAIN FEATURES:
    fuente = 'sans-serif'  # 'serif'
    small_size = 16
    medium_size = 28
    bigger_size = 26
    weight_selected = 'normal'
    line_width = 2
    marker_size = 10

    # MAJOR AND MINOR TICKS:
    # x_tick = int((max(presentation_matrix[:, 0]))) / 4
    # x_tick = int((max(presentation_matrix[:, 0]))) / 4
    # y_tick = int((max(presentation_matrix[:, :]))) / 4
    # x_tick = 20
    # y_tick = 1
    # ax.xaxis.set_major_locator(MultipleLocator(x_tick))
    # ax.yaxis.set_major_locator(MultipleLocator(y_tick))

    # x_tick_min = x_tick / 2
    # y_tick_min = y_tick / 2
    # ax.xaxis.set_minor_locator(MultipleLocator(x_tick_min))
    # ax.yaxis.set_minor_locator(MultipleLocator(y_tick_min))

    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # LIMIT TO AXIS:
    # ax.set_xlim(xmin=0, xmax=20)
    # ax.set_ylim(ymin=-0.5, ymax=9)

    # LINES:
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 1], 'r',
            label=r'$\mathregular{\Delta g_{xx}}$', linewidth=line_width)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 2], 'b',
             label=r'$\mathregular{\Delta g_{yy}}$', linewidth=line_width)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 3], 'k',
             label=r'$\mathregular{\Delta g_{zz}}$', linewidth=line_width)

    # MARKERS: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 1], 'ro', markersize=marker_size)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 2], 'bo', markersize=marker_size,
            markerfacecolor='none', markeredgewidth=1.5)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 3], 'ko', markersize=marker_size)

    # CHANGING THE FONTSIZE OF TICKS
    plt.xticks(fontsize=small_size, weight=weight_selected)
    plt.yticks(fontsize=small_size, weight=weight_selected)
    # axis.set_major_locator(MaxNLocator(integer=True))

    # LABELS:
    # labelpad: change the space between axis umbers and labels
    plt.xlabel(x_title, fontsize=bigger_size, fontfamily=fuente, labelpad=15,
               weight=weight_selected)
    plt.ylabel(y_title, fontsize=bigger_size, fontfamily=fuente, style='italic',
               weight=weight_selected, labelpad=15)
    # x_min = 0
    # x_max =  11
    # y_min = -45
    # y_max =  5
    # plt.xlim([x_min, x_max])  # Limit axis values
    # plt.ylim([y_min, y_max])  # Limit axis values

    # TITLE:
    # y = 1.05 change the space between title and plot
    plt.title(main_title, fontsize=bigger_size, fontfamily=fuente, y=1.05)

    # LEGEND
    legend = plt.legend(fontsize=medium_size, fancybox=True, framealpha=0.5,
                        labelcolor='linecolor', loc='center')
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    # plt.locator_params(nbins=10)
    # plt.grid()

    # Add an horizontal line in y = 0
    # ax.hlines(y=0, xmin=x_min, xmax=x_max, linewidth=line_width, color='k',
    #           linestyle='dotted')
    # dotted, dashed, solid, dashdot

    # Frame of the plot: https://e2eml.school/matplotlib_framing.html#spinestyle
    line_width = line_width - 0.8
    ax.spines["top"].set_linewidth(line_width)
    ax.spines["bottom"].set_linewidth(line_width)
    ax.spines["left"].set_linewidth(line_width)
    ax.spines["right"].set_linewidth(line_width)

    if save_picture == 0:
        plt.show()
        plt.close()
    else:
        figure_name = main_title + '_sos_analysis.png'
        plt.savefig(figure_name)


def sos_analysis_and_plot_mixinputs(file_msnull, file_ms_notnull, states_selected, states_option, ppms):
    """"
    Calculate the g-shifts in the sum-over-states_selected expansion using
    from 2 states_selected to the total number of states_selected shown in the Q-Chem output.
    :param: file_ms_notnull
    :return: no returned value, it prints the plot
    """
    totalstates = get_number_of_states(file_msnull)
    presentation_list = []

    nstates = get_selected_states(file_msnull, totalstates, states_selected, states_option, symmetry_selection=0)

    for i in range(1, len(nstates)+1):
        states_ras = nstates[0:i]

        g_shift = gfactor_sos_analysis(file_msnull, file_ms_notnull, states_ras, states_option, ppms)

        # if order_symmetry == 1:
        #     presentation_list.append([ordered_state_symmetries[i-1], np.round(g_shift.real[0], 3),
        #                           np.round(g_shift.real[1], 3), np.round(g_shift.real[2], 3)])
        # else:
        presentation_list.append([i, np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3),
                                  np.round(g_shift.real[2], 3)])
        print('Done: ', states_ras)
    presentation_matrix = np.array(presentation_list, dtype=object)

    # To presents deviation from previous g-values instead of the total g-values:
    # presentation_matrix_deviation = np.array(presentation_list, dtype=object)
    # for ndim in [1, 2, 3]:
    #     for i in range(1, len(presentation_matrix)):
    #         presentation_matrix_deviation[i, ndim] = (presentation_matrix[i, ndim])

    print("--------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("--------------------------------")
    # print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))

    plot_g_tensor_vs_states(presentation_matrix, x_title='Electronic State',
                            y_title=r'$\Delta g, ppm$', main_title=file_msnull, save_picture=0)
