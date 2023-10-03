"""
Calculation of the g-tensor using Q-Chem output with RAS-CI
"""
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator, MultipleLocator
from numpy import linalg, sqrt
from pyqchem.parsers.parser_rasci import parser_rasci
from parser_gtensor import get_number_of_states, get_symmetry_states, get_selected_states, order_by_states, \
    get_eigenenergies, get_spin_orbit_couplings, hermitian_test, get_hamiltonian_construction, reordering_eigenvectors, \
    diagonalization, get_spin_matrices, get_orbital_matrices, angular_matrixes_obtention, from_energies_soc_to_g_values, \
    print_g_calculation


def gfactor_exchange_energies_socs(file_ms_notnull, file_ms_null, states_ras, states_option):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, nstates, selected_states, symmetry_selection, soc_options
    :return: g-shifts
    """
    def get_energies_socs(file, nstates, states_option):
        """
        Having the selected states_selected, get the energy and SOCs between them
        :param: file, nstates, states_option
        :return: totalstates, eigenenergies_ras, selected_socs, sz_list, sz_ground
        """
        totalstates = get_number_of_states(file)

        symmetry_selections = 'None'
        nstates = get_selected_states(file, totalstates, nstates, states_option, symmetry_selections)

        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, nstates)
        selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(file, totalstates, nstates, soc_option=0)
        return totalstates, eigenenergies_ras, selected_socs, sz_list, sz_ground

    def mapping_between_states(ms_notnull_energies, ms_null_energies):
        """
        Comparing the energies, mapping between states_selected with Ms = 0, that do not have coupling between states_selected triplets,
        and states_selected with Ms not 0, that do have this coupling.
        :return:
        """
        mapping_list = []
        mapping_dict = {}

        for i in range(0, len(ms_notnull_energies)):
            for j in range(0, len(ms_null_energies)):
                ener_ms_notnull = np.round(ms_notnull_energies[i], 5)
                ener_ms_null = np.round(ms_null_energies[j], 5)

                if ener_ms_notnull == ener_ms_null:
                    mapping_dict = {'state ms not null': i, 'state ms null': j}
                    mapping_list.append(mapping_dict)

        # for mapping_dict in mapping_list:
        #     a = mapping_dict['state ms not null']
        #     b = mapping_dict['state ms null']
        #     print('file_ms_notnull', states_ras[a], 'file_ms_null', states_ras[b])
        # exit()
        return mapping_dict, mapping_list

    def exchange_coupling(mapping_list, selected_socs_ms_notnull, selected_socs_ms_null, sz_list):
        """
        Put SOCs between states_selected with Ms different than 0 (that are obtained in the output) in the
        SOC matrix of states_selected with Ms 0 (that are not obtained since Clebsh-Gordan coefficient is too small)
        :return:
        """
        for i in mapping_list:
            for j in mapping_list:
                if i['state ms not null'] != j['state ms not null']:  # If states_selected are not the same (in Ms not null list)
                    i_state_ms_notnull = i['state ms not null']
                    j_state_ms_notnull = j['state ms not null']

                    i_state_ms_null = i['state ms null']
                    j_state_ms_null = j['state ms null']

                    # print('State Ms not null', states_ras[i_state_ms_notnull], '(', i_state_ms_notnull, ')',
                    #       states_ras[j_state_ms_notnull], '(', j_state_ms_notnull, ')', '; ',
                    #       'State Ms null', states_ras[i_state_ms_null], '(', i_state_ms_null, ')',
                    #       states_ras[j_state_ms_null], '(', j_state_ms_null, ')',)

                    for sz_1 in range(0, len(sz_list)):
                        for sz_2 in range(0, len(sz_list)):
                            i_msnotnull = i_state_ms_notnull * len(sz_list) + sz_1
                            j_msnotnull = j_state_ms_notnull * len(sz_list) + sz_2

                            i_msnull = i_state_ms_null * len(sz_list) + sz_1
                            j_msnull = j_state_ms_null * len(sz_list) + sz_2

                            # print(socs_msnotnull[i_msnull, j_msnull], '<--->', socs_msnull[i_msnotnull, j_msnotnull])
                            selected_socs_ms_null[i_msnull, j_msnull] = selected_socs_ms_notnull[i_msnotnull, j_msnotnull]
                    # print('--END LOOP--')
        # print('SOC:')
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
        #                 for row in np.round((socs_msnotnull[:,:]),5)])) # * 219474.63068
        # exit()
        return selected_socs_ms_null

    # print('File ms not null: ', file_ms_notnull)
    # print('File ms null: ', file_ms_null)
    # print('States: ', states_ras)

    totalstates_ms_notnull, energies_ms_notnull, selected_socs_ms_notnull, sz_list_ms_notnull, sz_ground_ms_notnull = get_energies_socs(file_ms_notnull, states_ras, states_option)

    totalstates_ms_null, energies_ms_null, selected_socs_ms_null, sz_list_ms_null, sz_ground_ms_null = get_energies_socs(file_ms_null, states_ras, states_option)

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
