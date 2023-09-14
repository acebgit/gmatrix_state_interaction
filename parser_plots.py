import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from parser_gtensor import *


def plot_obtention(file, x_data, y_data):
    # "matplotlib" help: https://aprendeconalf.es/docencia/python/manual/matplotlib/
    # https://chartio.com/resources/tutorials/how-to-save-a-plot-to-a-file-using-matplotlib/
    # https://pythonspot.com/matplotlib-bar-chart/

    fuente = "Sans"
    size = 12

    y_pos = np.arange(len(x_data))
    plt.bar(y_pos, y_data, align='center', alpha=0.5, color='red')
    plt.xticks(y_pos, x_data)

    plt.title(file, fontsize=size, fontname=fuente)
    plt.ylabel('Energies (eV)', fontsize=size, fontname=fuente)
    plt.xlabel('number of state', fontsize=size, fontname=fuente)
    plt.axis([min(x_data) - 2, max(x_data), min(y_data), max(y_data) + 0.1 * max(y_data)])
    # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    # plt.grid(True)

    plt.plot()
    figure_name = file + '.png'
    plt.savefig(figure_name)

    # plt.show()
    plt.close()


def get_bar_chart(file, x_list, y_list, x_title, y_title, main_title):
    y_pos = (list(x_list))
    plt.bar(x_list, y_list, align='center', width=0.5, color='r', edgecolor="black")
    plt.xticks(y_pos)

    fuente = 'serif'

    plt.xlabel(x_title, fontsize=14, fontfamily=fuente)
    plt.ylabel(y_title, fontsize=14, fontfamily=fuente)
    plt.title(main_title, fontsize=18, fontfamily=fuente)

    plt.plot()
    figure_name = file + '_' + main_title + '.png'
    plt.savefig(figure_name)
    plt.show()
    plt.close()


def plot_g_tensor_vs_states(presentation_matrix, x_title, y_title, main_title, save_picture):
    fig, ax = plt.subplots()

    # Features to select:
    fuente = 'sans-serif'  # 'serif'
    small_size = 22
    medium_size = 22
    bigger_size = 24

    weight_selected = 'normal'
    line_width = 2
    marker_size = 10

    # x_min = 0
    # x_max =  11
    # y_min = -45
    # y_max =  5
    # x_tick = 1
    # y_tick = 5

    # for i in range(1, len(presentation_matrix[0, :])):
    #     maximum = max(presentation_matrix[:, i])
    #     if maximum > y_max:
    #         y_max = maximum + 10
    #
    # for i in range(1, len(presentation_matrix[0, :])):
    #     minimum = min(presentation_matrix[:, i])
    #     if minimum < y_min:
    #         y_min = minimum - 10

    # Major and minor ticks:
    # x_tick = int((max(presentation_matrix[:, 0]))) / 4
    # y_tick = 40
    # x_tick_min = x_tick / 2
    # y_tick_min = y_tick / 2

    # ax.xaxis.set_major_locator(MultipleLocator(x_tick))
    # ax.xaxis.set_minor_locator(MultipleLocator(x_tick_min))
    # ax.yaxis.set_major_locator(MultipleLocator(y_tick))
    # ax.yaxis.set_minor_locator(MultipleLocator(y_tick_min))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Lines:
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 1], 'r',
            label='$\\mathregular{\Delta g_{xx}}$', linewidth=line_width)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 2], 'b',
             label='$\mathregular{\Delta g_{yy}}$', linewidth=line_width)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 3], 'k',
             label='$\mathregular{\Delta g_{zz}}$', linewidth=line_width)

    # Markers: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 1], 'ro', markersize=marker_size)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 2], 'bo', markersize=marker_size,
            markerfacecolor='none', markeredgewidth=1.5)
    ax.plot(presentation_matrix[:, 0], presentation_matrix[:, 3], 'ko', markersize=marker_size)

    # changing the fontsize of yticks
    plt.xticks(fontsize=small_size, weight=weight_selected)
    plt.yticks(fontsize=small_size, weight=weight_selected)
    # axis.set_major_locator(MaxNLocator(integer=True))

    # Labels:
    # labelpad: change the space between axis umbers and labels
    plt.xlabel(x_title, fontsize=bigger_size, fontfamily=fuente, labelpad=15,
               weight=weight_selected)
    plt.ylabel(y_title, fontsize=bigger_size, fontfamily=fuente, style='italic',
               weight=weight_selected, labelpad=15)
    # plt.xlim([x_min, x_max])  # Limit axis values
    # plt.ylim([y_min, y_max])  # Limit axis values

    # Title:
    # y=1.05 change the space between title and plot
    plt.title(main_title, fontsize=bigger_size, fontfamily=fuente, y=1.05)

    # Legend
    legend = plt.legend(fontsize=medium_size, fancybox=True, framealpha=0.5,
                        labelcolor='linecolor', loc='center right')
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


def sos_analysis_and_plot(file, nstates):
    """"
    Calculate the g-shifts in the sum-over-states expansion using
    from 2 states to the total number of states shown in the Q-Chem output.
    :param: file
    :return: no returned value, it prints the plot
    """
    totalstates = get_number_of_states(file)
    presentation_list = []

    for i in range(1, len(nstates)+1):
        states_ras = nstates[0:i]
        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
        soc_options = 0
        selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_ras, soc_options)

        g_shift = from_energies_soc_to_g_values(file, states_ras,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, ground_sz)

        g_shift = g_shift * 1000
        # state_symmetries, ordered_state_symmetries = get_symmetry_states(file, nstates)
        # presentation_list.append([ordered_state_symmetries[i-1], np.round(
        #     ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3), np.round(ras_g_values.real[2], 3)])
        print('States used: ', states_ras)
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
    print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))

    plot_g_tensor_vs_states(presentation_matrix_deviation, x_title='Number of states',
                            y_title='$\Delta g, ppm$', main_title=file, save_picture=0)


def gfactor_all_states(file, nstates):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: ras_input, states_ras, selected_states, symmetry_selection, soc_options
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
        soc_options = 0
        selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_ras, soc_options)

        g_shift = from_energies_soc_to_g_values(file, states_ras,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, ground_sz)
        g_shift = g_shift * 1000

        print(nstates)
        presentation_list.append([nstates[0], np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3),
                                  np.round(g_shift.real[2], 3)])

    presentation_matrix = np.array(presentation_list, dtype=object)
    print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))
    print("----------------------------------------------------------")
    print(" G-TENSOR WITH DIFERENT GROUND STATES")
    print("----------------------------------------------------------")
    print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))

    print()
    presentation_matrix_2 = np.delete(presentation_matrix, 0, 0)
    print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix_2[:, :])]))
    plot_g_tensor_vs_states(presentation_matrix_2, x_title='Number of states',
                            y_title='$\Delta g, ppm$', main_title=file, save_picture=0)
