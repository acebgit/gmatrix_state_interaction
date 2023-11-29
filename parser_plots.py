__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator, MultipleLocator

from parser_gtensor import get_number_of_states, get_selected_states, get_eigenenergies, \
    get_spin_orbit_couplings, from_energies_soc_to_g_values, get_symmetry_states, from_ppt_to_ppm
from parser_excitstates import get_groundst_socc_values, get_groundst_orbital_momentum, \
    gshift_estimation_loop, gshift_calculation_loop


def save_picture(save_options, filee, title_main):
    """
    Function that shows the plot (save_options=0) or save it (save_options=1).
    :param save_options, file, main_title.
    :return: plot (saved or shown)
    """
    if save_options == 1:
        plt.plot()
        figure_name = filee + '_' + title_main + '.png'
        plt.savefig(figure_name)
        plt.close()
    else:
        plt.plot()
        plt.show()


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


def sos_analysis_and_plot(file, nstates, selected_state, ppms, estimation, order_symmetry, save_option):
    """"
    Calculate the g-shifts in the sum-over-states_selected expansion using
    from 2 states_selected to the total number of states_selected shown in the Q-Chem output.
    :param: file_ms_notnull
    :return: no returned value, it prints the plot
    """
    totalstates = get_number_of_states(file)
    nstates = get_selected_states(file, totalstates, nstates, selected_state, symmetry_selection=0)
    state_symmetries, ordered_state_symmetries = get_symmetry_states(file, nstates)

    presentation_list = []

    if estimation == 1:

        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, nstates)
        excitation_energies_ras[:] = (excitation_energies_ras[:] - excitation_energies_ras[0]) * 27.211399
        socc_values = get_groundst_socc_values(file, totalstates, nstates)
        orbitmoment_max, orbitmoment_all = get_groundst_orbital_momentum(file, totalstates, nstates)
        gxx_list, gyy_list, gzz_list = gshift_estimation_loop(nstates, orbitmoment_all, socc_values,
                                                              excitation_energies_ras, ppms)

        if order_symmetry == 1:
            for i in range(0, len(ordered_state_symmetries)):
                presentation_list.append([ordered_state_symmetries[i], np.round(gxx_list[i], 3),
                                          np.round(gyy_list[i], 3), np.round(gzz_list[i], 3)])
        else:
            for i in range(0, len(ordered_state_symmetries)):
                presentation_list.append([i, np.round(gxx_list[i][0], 3), np.round(gxx_list[i], 3),
                                          np.round(gyy_list[i], 3), np.round(gzz_list[i], 3)])

    else:

        for i in range(0, len(nstates)):
            states_sos = nstates[0:i+1]
            print('Calculating....', states_sos)
            eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_sos)
            selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_sos, soc_option=0)
            g_shift = from_energies_soc_to_g_values(file, states_sos,
                                                    totalstates, excitation_energies_ras,
                                                    selected_socs, sz_list, ground_sz, ppms)

            if order_symmetry == 1:
                presentation_list.append([ordered_state_symmetries[i], np.round(g_shift.real[0], 3),
                                          np.round(g_shift.real[1], 3), np.round(g_shift.real[2], 3)])
            else:
                presentation_list.append([i+1, np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3),
                                          np.round(g_shift.real[2], 3)])
    presentation_matrix = np.array(presentation_list, dtype=object)

    print("--------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("--------------------------------")
    file = file[:-4]
    x_title = 'Electronic State'
    y_title = r'$\Delta g, ppt$'
    if ppms == 1:
        y_title = r'$\Delta g, ppm$'
    main_title = 'sos_analysis'
    plot_g_tensor_vs_states(file, presentation_matrix, x_title, y_title,
                            main_title, save_option)


def gfactor_all_states(file, nstates, ppms):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, states_msnull, states_option, symmetry_selection, soc_options
    :return: g-shifts
    """
    def swap_positions(lista, pos1, pos2):
        lista[pos1], lista[pos2] = lista[pos2], lista[pos1]
        return lista

    totalstates = get_number_of_states(file)
    presentation_list = ['Ground state', 'gxx', 'gyy', 'gzz']

    for i in range(0, len(nstates)):
        states_ras = swap_positions(nstates, 0, i)

        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)

        selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_ras, soc_option=0)

        g_shift = from_energies_soc_to_g_values(file, states_ras,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, ground_sz, ppms)

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


def compare_gcalculation_gestimation(file, nstates, selected_state, ppms, plotting):
    """
    Make the comparison betweem g-shift calculated with the SOS procedure and the estimation equation
    (equation 11 or article J. Phys.Chem.A2023, 127, 8459−8472).
    :param file:
    :param nstates:
    :param selected_state:
    :param plotting:
    :return:
    """
    def plot_g_tensor_comparison(filee, present_matrix, title_x, title_y, title_main, save_options):
        """
        Plotting results. 
        :param filee: 
        :param present_matrix: 
        :param title_x: 
        :param title_y: 
        :param title_main: 
        :param save_options: 
        :return: 
        """
        fig, ax = plt.subplots()

        # MAIN FEATURES:
        fuente = 'sans-serif'  # 'serif'
        small_size = 16
        medium_size = 28
        bigger_size = 26
        weight_selected = 'normal'
        line_width = 2
        marker_size = 10

        # LINES:
        ax.plot(present_matrix[:, 0], present_matrix[:, 1], 'r',
                label=r'Calc. $\mathregular{\Delta g_{xx}}$', linewidth=line_width)
        ax.plot(present_matrix[:, 0], present_matrix[:, 2], 'b',
                label=r'Calc. $\mathregular{\Delta g_{yy}}$', linewidth=line_width)
        ax.plot(present_matrix[:, 0], present_matrix[:, 3], 'k',
                label=r'Calc. $\mathregular{\Delta g_{zz}}$', linewidth=line_width)

        ax.plot(present_matrix[:, 0], present_matrix[:, 4], 'r--',
                label=r'Estim. $\mathregular{\Delta g_{xx}}$', linewidth=line_width)
        ax.plot(present_matrix[:, 0], present_matrix[:, 5], 'b--',
                label=r'Estim. $\mathregular{\Delta g_{yy}}$', linewidth=line_width)
        ax.plot(present_matrix[:, 0], present_matrix[:, 6], 'k--',
                label=r'Estim. $\mathregular{\Delta g_{zz}}$', linewidth=line_width)

        # MARKERS: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.markers.MarkerStyle.html
        ax.plot(present_matrix[:, 0], present_matrix[:, 1], 'r^', markersize=marker_size,
                markerfacecolor=None, markeredgewidth=1.5, fillstyle='none')
        ax.plot(present_matrix[:, 0], present_matrix[:, 2], 'b^', markersize=marker_size,
                markerfacecolor=None, markeredgewidth=1.5, fillstyle='none')
        ax.plot(present_matrix[:, 0], present_matrix[:, 3], 'k^', markersize=marker_size,
                markerfacecolor=None, markeredgewidth=1.5, fillstyle='none')

        ax.plot(present_matrix[:, 0], present_matrix[:, 4], 'rs', markersize=marker_size,
                markerfacecolor=None, markeredgewidth=1.5, fillstyle='none')
        ax.plot(present_matrix[:, 0], present_matrix[:, 5], 'bs', markersize=marker_size,
                markerfacecolor=None, markeredgewidth=1.5, fillstyle='none')
        ax.plot(present_matrix[:, 0], present_matrix[:, 6], 'ks', markersize=marker_size,
                markerfacecolor=None, markeredgewidth=1.5, fillstyle='none')

        # CHANGING THE FONTSIZE OF TICKS
        plt.xticks(fontsize=small_size, weight=weight_selected)
        plt.yticks(fontsize=small_size, weight=weight_selected)
        # axis.set_major_locator(MaxNLocator(integer=True))

        # LABELS:
        # labelpad: change the space between axis umbers and labels
        plt.xlabel(title_x, fontsize=bigger_size, fontfamily=fuente, labelpad=15,
                   weight=weight_selected)
        plt.ylabel(title_y, fontsize=bigger_size, fontfamily=fuente, style='italic',
                   weight=weight_selected, labelpad=15)

        # TITLE:
        plt.title(title_main, fontsize=bigger_size, fontfamily=fuente, y=1.05)

        # LEGEND
        legend = plt.legend(fontsize=small_size, fancybox=True, framealpha=0.5,
                            labelcolor='linecolor', loc='upper right')
        frame = legend.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')
        line_width = line_width - 0.8
        ax.spines["top"].set_linewidth(line_width)
        ax.spines["bottom"].set_linewidth(line_width)
        ax.spines["left"].set_linewidth(line_width)
        ax.spines["right"].set_linewidth(line_width)

        save_picture(save_options, filee, title_main)

    totalstates = get_number_of_states(file)
    nstates = get_selected_states(file, totalstates, nstates, selected_state, symmetry_selection=0)
    state_symmetries, ordered_state_symmetries = get_symmetry_states(file, nstates)

    # 1) Make the estimation of the g-shift separately for each state
    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, nstates)
    excitation_energies_ras = [(i - excitation_energies_ras[0])*27.211399 for i in excitation_energies_ras]
    soccs = get_groundst_socc_values(file, totalstates, nstates)
    orbitmoment_max, orbitmoment_all = get_groundst_orbital_momentum(file, totalstates, nstates)

    estim_gxx, estim_gyy, estim_gzz = gshift_estimation_loop(nstates, orbitmoment_all, soccs, excitation_energies_ras, ppms)
    estim_gshift = []
    for i in range(0, len(estim_gxx)):
        estim_gshift.append([estim_gxx[i], estim_gyy[i], estim_gzz[i]])

    # 2) Make the calculation of the g-shift separately for each state
    calc_gxx, calc_gyy, calc_gzz = gshift_calculation_loop(nstates, file, totalstates, ppms=0)
    calc_gshift = []
    for i in range(0, len(calc_gxx)):
        calc_gshift.append([calc_gxx[i], calc_gyy[i], calc_gzz[i]])

    # 3) Make the comparison
    # https://pythonforundergradengineers.com/unicode-characters-in-python.html
    presentation_list = []
    dec = 5
    for i in range(0, len(nstates)):
        presentation_list.append([np.round(calc_gshift[i][0], dec),
                                  np.round(calc_gshift[i][1], dec), np.round(calc_gshift[i][2], dec),
                                  np.round(estim_gshift[i][0], dec), np.round(estim_gshift[i][1], dec),
                                  np.round(estim_gshift[i][2], dec)])

    print("--------------------------------")
    print("Comparison calculated-estimated g-shift")
    print("--------------------------------")
    # pd.set_option('display.max_colwidth', 100)
    pd.set_option('display.width', 400)
    pd.set_option('display.max_columns', 10)
    df = pd.DataFrame(presentation_list, index=ordered_state_symmetries,
                      columns=['Calc. \u0394gxx', 'Calc. \u0394gyy', 'Calc. \u0394gzz',
                               'Estim. \u0394gxx', 'Estim. \u0394gyy', 'Estim. \u0394gzz'])
    print(df)

    if plotting == 1:
        file_string = file[:-4]
        x_title = 'Electronic State'
        y_title = r'$\Delta g, ppt$'
        main_title = 'comparison_calc_estim'

        presentation_list = []
        dec = 5
        for i in range(0, len(nstates)):
            presentation_list.append([ordered_state_symmetries[i], np.round(calc_gshift[i][0], dec),
                                      np.round(calc_gshift[i][1], dec), np.round(calc_gshift[i][2], dec),
                                      np.round(estim_gshift[i][0], dec), np.round(estim_gshift[i][1], dec),
                                      np.round(estim_gshift[i][2], dec)])

        presentation_matrix = np.array(presentation_list, dtype=object)
        presentation_matrix_2 = np.delete(presentation_matrix, 0, 0)
        plot_g_tensor_comparison(file_string, presentation_matrix_2, x_title,
                                 y_title, main_title, save_options=0)
