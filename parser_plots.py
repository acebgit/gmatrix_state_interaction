__author__ = 'Antonio Cebreiro-Gallardo'

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator

from parser_gtensor import *


def save_picture(save_options, file, main_title):
    """
    Function that shows the plot (save_options=0) or save it (save_options=1).
    :param save_options, file, main_title.
    :return: plot (saved or shown)
    """
    if save_options == 1:
        plt.plot()
        figure_name = file + '_' + main_title + '.png'
        plt.savefig(figure_name)
        plt.close()
    else:
        plt.plot()
        plt.show()


def get_bar_chart(file, x_list, y_list, x_title, y_title, main_title, save_pict):
    """
    Print Bar plots: y_list vs x_list.
    :param:
    :return:
    """
    # "matplotlib" help: https://aprendeconalf.es/docencia/python/manual/matplotlib/
    # https://chartio.com/resources/tutorials/how-to-save-a-plot-to-a-file-using-matplotlib/
    # https://pythonspot.com/matplotlib-bar-chart/
    fuente = 'serif'  # "Sans"
    medium_size = 16
    big_size = 20

    y_pos = (list(x_list))
    plt.bar(x_list, y_list, align='center', width=0.5, color='r', edgecolor="black")
    plt.xticks(y_pos)

    plt.title(main_title, fontsize=big_size, fontname=fuente)
    plt.xlabel(x_title, fontsize=medium_size, fontfamily=fuente)
    plt.ylabel(y_title, fontsize=medium_size, fontfamily=fuente)
    # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    plt.grid(True)

    save_picture(save_pict, file, main_title)


def plot_g_tensor_vs_states(file, presentation_matrix, x_title, y_title, main_title, save_options):
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

    save_picture(save_options, file, main_title)


def sos_analysis_and_plot(file, nstates, selected_state, ppms, order_symmetry, save_option):
    """"
    Calculate the g-shifts in the sum-over-states_selected expansion using
    from 2 states_selected to the total number of states_selected shown in the Q-Chem output.
    :param: file_ms_notnull
    :return: no returned value, it prints the plot
    """
    totalstates = get_number_of_states(file)
    presentation_list = []

    nstates = get_selected_states(file, totalstates, nstates, selected_state, symmetry_selection=0)

    for i in range(1, len(nstates)+1):
        states_ras = nstates[0:i]
        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
        soc_options = 0
        selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_ras, soc_options)

        g_shift = from_energies_soc_to_g_values(file, states_ras,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, ground_sz)

        g_shift = from_ppt_to_ppm(g_shift, ppms)

        state_symmetries, ordered_state_symmetries = get_symmetry_states(file, nstates)
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
    # print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))

    print("--------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("--------------------------------")
    file = file[:-4]
    x_title = 'Electronic State'
    y_title = r'$\Delta g, ppm$'
    main_title = 'sos_analysis'

    plot_g_tensor_vs_states(file, presentation_matrix_deviation, x_title, y_title,
                            main_title, save_option)


def gfactor_all_states(file, nstates, ppms):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, states_msnull, states_option, symmetry_selection, soc_options
    :return: g-shifts
    """
    def swapPositions(list, pos1, pos2):
        list[pos1], list[pos2] = list[pos2], list[pos1]
        return list

    totalstates = get_number_of_states(file)
    presentation_list = []
    presentation_list.append(['Ground state', 'gxx', 'gyy', 'gzz'])

    for i in range(0, len(nstates)):
        states_ras = swapPositions(nstates, 0, i)

        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)

        selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_ras, soc_option=0)

        g_shift = from_energies_soc_to_g_values(file, states_ras,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, ground_sz)

        g_shift = from_ppt_to_ppm(g_shift, ppms)

        presentation_list.append([nstates[0], np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3),
                                  np.round(g_shift.real[2], 3)])

    presentation_matrix = np.array(presentation_list, dtype=object)
    print("-----------------------------------------------")
    print(" G-TENSOR WITH DIFFERENT GROUND STATES")
    print("-----------------------------------------------")
    print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))

    presentation_matrix_2 = np.delete(presentation_matrix, 0, 0)
    plot_g_tensor_vs_states(file, presentation_matrix_2, x_title='Number of states_selected',
                            y_title=r'$\Delta g, ppm$', main_title=file, save_options=0)
