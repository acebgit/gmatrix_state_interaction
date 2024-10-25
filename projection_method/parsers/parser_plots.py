__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator, MultipleLocator

from projection_method.parsers.parser_gtensor import get_number_of_states, get_selected_states, get_eigenenergies, \
    get_spin_orbit_couplings, from_energies_soc_to_g_values, get_symmetry_states
from projection_method.parsers.parser_excitstates import get_groundst_socc_values, get_groundst_orbital_momentum, \
    gshift_estimation_loop, gshift_calculation_loop

# Get the absolute path to local folder 'PyQChem'
import sys
import os
parent_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), 'C:\\Users\\HP.LAPTOP-F127N3L3\\Desktop\\my_programs\\PyQChem'))
if parent_directory not in sys.path: # Add 'folder2' to the Python path
    sys.path.append(parent_directory)
from pyqchem.parsers.parser_rasci import parser_rasci


def save_picture(save_options, filee, title_main):
    """
    Function that shows the plot (save_options=0) or save it (save_options=1).
    :param save_options, file, main_title.
    :return: plot (saved or shown)
    """
    if save_options == 1:
        plt.tight_layout()  # Adjust layout
        plt.plot()
        figure_name = filee + '_' + title_main + '.png'
        plt.savefig(figure_name)
        plt.show()
        plt.close()
    else:
        plt.tight_layout()  # Adjust layout
        plt.plot()
        plt.show()


def plot_g_tensor_vs_states(file, presentation_matrix, x_title, y_title, main_title, save_options):
    fig, ax = plt.subplots()
    plot_type = 1 # 0: plot, 1: bars

    # MAIN FEATURES:
    fuente = 'sans-serif'  # 'serif'
    small_size = 16
    medium_size = 16
    bigger_size = 16
    weight_selected = 'normal'
    line_width = 2
    marker_size = 10

    x = presentation_matrix[:, 0]  # First column for x-axis
    y1 = presentation_matrix[:, 1]  # Second column for the first category
    y2 = presentation_matrix[:, 2]  # Third column for the second category
    y3 = presentation_matrix[:, 3]  # Fourth column for the third category

    #################################
    ###   PLOT TYPE
    #################################
    if plot_type == 0:
        # MAJOR AND MINOR TICKS:
        # x_tick = int((max(x))) / 4
        # x_tick = int((max(x))) / 4
        # y_tick = int((max(x))) / 4
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
        # ax.plot(x, y1, 'r',
        #         label=r'$\mathregular{\Delta g_{xx}}$', linewidth=line_width)
        # ax.plot(x, y2, 'b', 
        #         label=r'$\mathregular{\Delta g_{yy}}$', linewidth=line_width)
        # ax.plot(x, y3, 'k',
        #         label=r'$\mathregular{\Delta g_{zz}}$', linewidth=line_width)

        # MARKERS: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
        ax.plot(x, y1, 'ro', markersize=marker_size)
        ax.plot(x, y2, 'bv', markersize=marker_size,
                markerfacecolor='none', markeredgewidth=1.5)
        ax.plot(x, y3, 'ks', markersize=marker_size)

        # CHANGING THE FONTSIZE OF TICKS
        # plt.xticks(fontsize=small_size, weight=weight_selected)
        # plt.yticks(fontsize=small_size, weight=weight_selected)
        # axis.set_major_locator(MaxNLocator(integer=True))

        # Enable grid lines for both x and y axes
        plt.grid(axis='both', linestyle='--', linewidth=0.7, alpha=0.7)

    #################################
    ###   BAR PLOTS
    #################################
    elif plot_type == 1:
        # Set width of the bars
        bar_width = 0.25

        # Set the positions of the bars on the x-axis
        r1 = np.arange(len(x))
        r2 = [x + bar_width for x in r1]
        r3 = [x + bar_width * 2 for x in r1]

        # Create the bar plot
        plt.bar(r1, y1, width=bar_width, color='blue', edgecolor='blue', label=r'$\mathregular{\Delta g_{xx}}$')
        plt.bar(r2, y2, width=bar_width, color='orange', edgecolor='orange', label=r'$\mathregular{\Delta g_{yy}}$')
        plt.bar(r3, y3, width=bar_width, color='green', edgecolor='green', label=r'$\mathregular{\Delta g_{zz}}$')

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
                        labelcolor='linecolor', loc='best')
    frame = legend.get_frame()
    frame.set_facecolor('white')
    # frame.set_edgecolor('black')

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


def sos_analysis_and_plot(file, nstates, selected_state, analysis, ppms, order_symmetry, save_option):
    """"
    Calculate the g-shifts in the sum-over-states_selected expansion using pair of states: 
    Ground state - Excited states with all the excited states of the output.
    g-shift in absolute value to make the estimation.
    Two options:
    i) "analysis = 0": Effective Hamiltonian 
    ii) "analysis = 1": g-shift calculated with the SOS procedure and the estimation equation 
    (eq. (11) or article J. Phys.Chem.A2023, 127, 8459-8472). Estimation is done with 
    (L SOC / \Delta E), since only the absolute value gives "real" information. 
    :param: file
    :return: prints the plot
    """
    with open(file, encoding="utf8") as f:
        output = f.read()
    output_parsered = parser_rasci(output)

    totalstates = get_number_of_states(file)
    nstates = get_selected_states(file, totalstates, nstates, selected_state, symmetry_selection=0)
    state_symmetries, ordered_state_symmetries = get_symmetry_states(file, totalstates)
    states_selected_symmetries = [ordered_state_symmetries[i-1] for i in nstates]
    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, nstates)
    selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(output_parsered, totalstates, nstates, soc_option=0, soc_order=0)

    print(len(selected_socs))
    print(len(sz_list))
    print(ground_sz)
    print(len(nstates))
    exit()

    presentation_list = []
    gxx_list = []
    gyy_list = []
    gzz_list = []

    if analysis == 0:
        for i in range(0, len(nstates)):
            print('SOS: state ', i)
            states_sos = [nstates[0], nstates[i]]
            excitener_pairstates = [excitation_energies_ras[0], excitation_energies_ras[i]]
            # st_selected_socs, st_sz_list, st_ground_sz = 
            g_shift = from_energies_soc_to_g_values(file, states_sos,
                                                    totalstates, excitener_pairstates,
                                                    st_selected_socs, st_sz_list, st_ground_sz, ppms)
            gxx_list.append(abs(g_shift[0].real))
            gyy_list.append(abs(g_shift[1].real))
            gzz_list.append(abs(g_shift[2].real))

    elif analysis == 1:
        # Energy: from Hartree to a.u.
        excitation_energies_ras[:] = (excitation_energies_ras[:] - excitation_energies_ras[0]) * 27.211399
        socc_values = get_groundst_socc_values(file, totalstates, nstates)
        orbitmoment_max, orbitmoment_all = get_groundst_orbital_momentum(file, totalstates, nstates)
        gxx_list, gyy_list, gzz_list = gshift_estimation_loop(nstates, orbitmoment_all, socc_values,
                                                            excitation_energies_ras, ppms)
    
    if order_symmetry == 1:
        for i in range(0, len(states_selected_symmetries)):
            presentation_list.append([ordered_state_symmetries[i], np.round(gxx_list[i], 3),
                                        np.round(gyy_list[i], 3), np.round(gzz_list[i], 3)])
    else:
        for i in range(0, len(states_selected_symmetries)):
            presentation_list.append([i+1, np.round(gxx_list[i], 3),
                                        np.round(gyy_list[i], 3), np.round(gzz_list[i], 3)])
    presentation_matrix = np.array(presentation_list, dtype=object)

    print("--------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("--------------------------------")

    # Set display options to show all rows and columns
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame([row[0:4] for row in presentation_list],[row[0] for row in presentation_list],columns=['state','gxx','gyy','gzz'])
    print(df.to_string(index=False))
    print()

    file = file[:-4]
    x_title = 'Electronic State'
    y_title = r'$\Delta g, ppt$'
    if ppms == 1:
        y_title = r'$\Delta g, ppm$'
    main_title = 'sos_analysis'
    plot_g_tensor_vs_states(file, presentation_matrix, x_title, y_title,
                            main_title, save_option)


def gfactor_change_ground_state(filee, initial_states, states_option, symmetry_selection, soc_option, ppms):
    """
    Returns the g-shifts for doublet ground state molecules.
    :param: file_ms_notnull, states_msnull, states_option, symmetry_selection, soc_options
    :return: g-shifts
    """
    def swap_positions(lista, pos1, pos2):
        lista[pos1], lista[pos2] = lista[pos2], lista[pos1]
        return lista

    totalstates = get_number_of_states(filee)

    states_selected = get_selected_states(filee, totalstates, initial_states, states_option, symmetry_selection)

    presentation_list = []

    for i in range(0, len(states_selected)):
        istates = swap_positions(states_selected, 0, i)

        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(filee, totalstates, istates)

        selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(filee, totalstates, istates, soc_option)

        g_shift = from_energies_soc_to_g_values(filee, istates,
                                                totalstates, excitation_energies_ras,
                                                selected_socs, sz_list, ground_sz, ppms)

        presentation_list.append([istates[0], np.round(g_shift[0].real, 3), np.round(g_shift[1].real, 3),np.round(g_shift[2].real, 3)])
    
    presentation_matrix = np.array(presentation_list, dtype=object) 
          
    print("-----------------------------------------------")
    print(" G-TENSOR WITH DIFFERENT GROUND STATES")
    print("-----------------------------------------------")
    pd.set_option('display.width', 400)
    pd.set_option('display.max_columns', 10)
    df = pd.DataFrame(presentation_matrix,columns=['Ground state', 'gxx', 'gyy', 'gzz'])
    print(df.to_string(index=False))
    print()


def compare_gcalculation_gestimation(file, nstates, selected_state, ppms, plotting):
    """
    Make the comparison betweem g-shift calculated with the SOS procedure and the estimation equation
    (eq. (11) or article J. Phys.Chem.A2023, 127, 8459−8472).
    Estimation is done with (4 L SOC) / \Delta E, since only the absolute value gives information. 
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
