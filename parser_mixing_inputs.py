"""
Calculation of the g-tensor using Q-Chem output with RAS-CI
"""
import numpy as np
import matplotlib.pyplot as plt

from parser_gtensor import get_number_of_states, get_symmetry_states, get_selected_states, get_eigenenergies, \
    get_spin_orbit_couplings, from_energies_soc_to_g_values, get_spin_matrices, get_orbital_matrices, hermitian_test


def mapping_between_states(file_msnull, file_msnotnull, states_ras, states_option, states_sym):
    """
    Comparing the energies, map states_selected with Ms = 0, that do not have coupling between states_selected triplets,
    and states_selected with Ms not 0, that do have this coupling.
    :return: mapping_list
    """
    def get_all_energies(file, states_ras, states_option, states_sym):
        """
        Get all energies of RAS-CI states.
        :param: file
        :return: eigenenergies_ras (in atomic units)
        """
        totalstates = get_number_of_states(file)
        selected_states = get_selected_states(file, totalstates, states_ras, states_option, states_sym)
        eigenenergies, excitation_energies = get_eigenenergies(file, totalstates, selected_states)
        return eigenenergies, selected_states

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

    ms_notnull_energies, states_msnotnull = get_all_energies(file_msnotnull, states_ras, states_option, states_sym)
    ms_null_energies, states_msnull = get_all_energies(file_msnull, states_ras, states_option, states_sym)

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
        print(states_msnotnull[a], ' - ', states_msnotnull[b])
    print('---')
    # exit()
    return mapping_dict, mapping_list


def get_input_values(ras_input, states_ras, selected_states, symmetry_selection, soc_options):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, states_ras, selected_states, symmetry_selection, soc_options
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
        :param: totalstates_1, totalstates_2, list_mapping
        :return: totalstates
        """
        totalstates = len(states_1) + len(states_2) - len(list_mapping)
        print('totalstates: ', totalstates)
        print('---')
        return totalstates


def eigenenergy_mix(eigenenergies_1, eigenenergies_2, mapping_list):
        """
        Give the energies of the states in eigenenergies_1 and eigenenergies_2 without those states that are
        repeat in both lists (meaning those that are in the mapping list)
        :param: eigenenergies_1, eigenenergies_2, mapping_list
        :return: eigenenergies
        """
        elements = list(eigenenergies_1)

        repeat_states = []
        for mapping_dict in mapping_list:  # States that are in eigenenergies_1 and eigenenergies_2
            repeat_states.append(mapping_dict['state ms null'])

        for state in range(0, len(eigenenergies_2)):  # States that are in eigenenergies_2
                if state not in repeat_states:
                    elements.append(eigenenergies_2[state])

        eigenenergies = np.array(elements)
        print('Eigenenergies: ', eigenenergies)
        print('---')
        return eigenenergies


def socs_mix(socs_msnull, socs_msnotnull, mapping_list, sz_list, totalstates):
        """
        Give the SOCS of the states in socs_msnull and soc_msnotnull without those states that are
        repeat in both lists (meaning those that are in the mapping list)
        :param: socs_msnull, socs_msnotnull, mapping_list, sz_list, totalstates
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
            #                  for row in np.round((total_socs[:, :] * 219474.63068), 5)]))  # * 219474.63068
            # print()
            return total_soc

        def get_no_mapped_states(list_mapping, nstate):
            """
            Obtain a list of not mapped states
            :return: not_mapped_states
            """
            mapped_states = []
            for mapping_dict in list_mapping:
                mapped_states.append(mapping_dict['state ms not null'])

            not_mapped_states = []
            for state in range(0, nstate):
                if state not in mapped_states:
                    not_mapped_states.append(state)
            return not_mapped_states

        def include_msnotnull_states(nstate, not_mapped_states, soc_msnotnull, total_socs):
            """
            SOCs between states Ms ≠ 0 that have not been included in total SOC matrix are now included
            :return:
            """
            for state_i in range(0, nstate):
                for state_j in range(0, nstate):

                    if (state_i in not_mapped_states) or (state_j in not_mapped_states):
                        state_i_total = nstate
                        state_j_total = nstate

                        if (state_i in not_mapped_states) and (state_j not in not_mapped_states):
                            state_i_total += not_mapped_states.index(state_i)
                            state_j_total = state_j
                        elif (state_i not in not_mapped_states) and (state_j in not_mapped_states):
                            state_i_total = state_i
                            state_j_total += not_mapped_states.index(state_j)
                        elif (state_i in not_mapped_states) and (state_j in not_mapped_states):
                            state_i_total += not_mapped_states.index(state_i)
                            state_j_total += not_mapped_states.index(state_j)

                        # print('States:', state_i, state_j, ', states total: ', state_i_total, state_j_total)

                        for sz_1 in range(0, len(sz_list)):
                            for sz_2 in range(0, len(sz_list)):
                                soc_row = state_i * len(sz_list) + sz_1
                                soc_col = state_j * len(sz_list) + sz_2

                                soc_row_final = state_i_total * len(sz_list) + sz_1
                                soc_col_final = state_j_total * len(sz_list) + sz_2
                                # print(soc_row, soc_col, '-->', soc_row_final, soc_col_final)

                                total_socs[soc_row_final, soc_col_final] = soc_msnotnull[soc_row, soc_col]
            return total_socs

        # 1) First block of SOC total matrix is the SOCs of Ms = 0
        total_soc = np.zeros((totalstates * len(sz_list), totalstates * len(sz_list)), dtype=complex)
        for i in range(0, len(socs_msnull)):
            for j in range(0, len(socs_msnull)):
                total_soc[i, j] = socs_msnull[i, j]

        # 2) SOCs in total matrix between states mapped (those in which SOC is not calculated because Clebsch-Gordan coefficient
        # not calculated) are exchanged by those calculated couplings in Ms ≠ 0.
        total_soc = exchanging_socs(total_soc, socs_msnotnull, mapping_list, sz_list)

        # 3) Obtain list of not mapped states
        nstates = int(len(socs_msnotnull) / len(sz_list))
        not_mapped_states_msnotnull = get_no_mapped_states(mapping_list, nstates)

        # 4) Those SOCs between states Ms ≠ 0 that have not been included in total SOC matrix are now included
        total_soc = include_msnotnull_states(nstates, not_mapped_states_msnotnull, socs_msnotnull, total_soc)

        # print('socs_msnull:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                 for row in np.round((socs_msnull[:,:]* 219474.63068),5)])) # * 219474.63068
        # print('---')
        # print('socs_msnotnull:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                 for row in np.round((socs_msnotnull[:,:]* 219474.63068),5)])) # * 219474.63068
        # print('---')
        # print('Total SOC:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                 for row in np.round((total_soc[:,:]* 219474.63068),5)])) # * 219474.63068
        # exit()
        return total_soc


def angular_momentums_mix(angular_matrix_1, angular_matrix_2, sz_list, totalstates):
    """

    :param angular_matrix_1:
    :param angular_matrix_2:
    :param sz_list:
    :param totalstates:
    :return:
    """


def gfactor_exchange_energies_socs(file_ms_notnull, file_ms_null, states_ras, states_option):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, selected_states, selected_states, symmetry_selection, soc_options
    :return: g-shifts
    """
    totalstates_ms_notnull, energies_ms_notnull, selected_socs_ms_notnull, sz_list_ms_notnull, sz_ground_ms_notnull \
        = get_energies_socs(file_ms_notnull, states_ras, states_option)

    totalstates_ms_null, energies_ms_null, selected_socs_ms_null, sz_list_ms_null, sz_ground_ms_null \
        = get_energies_socs(file_ms_null, states_ras, states_option)

    mapping_dict, mapping_list = mapping_between_states(energies_ms_notnull, energies_ms_null)

    selected_socs_ms_null = exchange_coupling(mapping_list, selected_socs_ms_notnull, selected_socs_ms_null, sz_list_ms_notnull)

    return energies_ms_null, selected_socs_ms_null, sz_list_ms_null, sz_ground_ms_null, totalstates_ms_null
    #
    # hamiltonian_ras = get_hamiltonian_construction(states_ras, energies_ms_null, socs_msnotnull, sz_list_ms_null)
    #
    # eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian_ras)
    #
    # spin_matrix, standard_spin_matrix = get_spin_matrices(file_ms_notnull, states_ras)
    #
    # orbital_matrix = get_orbital_matrices(file_ms_notnull, totalstates_ms_null, states_ras, sz_list_ms_null)
    #
    # combination_spin_matrix = angular_matrixes_obtention(eigenvector, spin_matrix, sz_list_ms_null)
    #
    # combination_orbital_matrix = angular_matrixes_obtention(eigenvector, orbital_matrix, sz_list_ms_null)
    #
    # g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
    #                                sz_list_ms_null, sz_ground_ms_null)
    #
    # print_g_calculation(file_ms_notnull, totalstates_ms_null, states_ras, states_ras, g_shift * 1000, symmetry_selection=0)


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


def sos_analysis_and_plot(file_ms_notnull, file_ms_null, nstates, selected_state, order_symmetry):
    """"
    Calculate the g-shifts in the sum-over-states_selected expansion using
    from 2 states_selected to the total number of states_selected shown in the Q-Chem output.
    :param: file_ms_notnull
    :return: no returned value, it prints the plot
    """
    totalstates = get_number_of_states(file_ms_null)
    presentation_list = []

    nstates = get_selected_states(file_ms_null, totalstates, nstates, selected_state, symmetry_selection=0)

    for i in range(1, len(nstates)+1):
        states_ras = nstates[0:i]
        # eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
        # soc_options = 0
        # selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_ras, soc_options)
        excitation_energies_ras, selected_socs, sz_list, ground_sz, totalstates \
            = gfactor_exchange_energies_socs(file_ms_notnull, file_ms_null, states_ras, selected_state)

        g_shift = from_energies_soc_to_g_values(file_ms_null, states_ras,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, ground_sz)

        g_shift = g_shift * 1000
        state_symmetries, ordered_state_symmetries = get_symmetry_states(file_ms_null, nstates)
        if order_symmetry == 1:
            presentation_list.append([ordered_state_symmetries[i-1], np.round(g_shift.real[0], 3),
                                  np.round(g_shift.real[1], 3), np.round(g_shift.real[2], 3)])
        else:
            presentation_list.append([i, np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3),
                                  np.round(g_shift.real[2], 3)])
    presentation_matrix = np.array(presentation_list, dtype=object)

    # To presents deviation from previous g-values instead of the total g-values:
    presentation_matrix_deviation = np.array(presentation_list, dtype=object)
    for ndim in [1, 2, 3]:
        for i in range(1, len(presentation_matrix)):
            presentation_matrix_deviation[i, ndim] = (presentation_matrix[i, ndim])

    print("--------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("--------------------------------")
    # print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))

    plot_g_tensor_vs_states(presentation_matrix_deviation, x_title='Electronic State',
                            y_title=r'$\Delta g, ppm$', main_title=file_ms_null, save_picture=0)
