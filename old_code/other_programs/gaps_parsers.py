import numpy as np
import matplotlib.pyplot as plt


def get_energy_and_spins(files):
    """
    Obtain list of excitation energies and spins of all the RAS states.
    :param: files
    :return: energy_list, spins_list
    """
    energy_list = []
    dft_ener_list = []
    spins_list = []

    word_search = [' Excitation energy (eV) = ', ' <S^2>      :',
                   ' Excitation energy (eV) with DFT-0']

    for line in files:
        if word_search[0] in line:
            line = line.split()
            element = float(line[-1])
            energy_list.append(element)

        elif word_search[1] in line:
            line = line.split()
            element = float(line[-1])
            spins_list.append(element)

        elif word_search[2] in line:
            line = line.split()
            element = float(line[-1])
            dft_ener_list.append(element)
    return energy_list, spins_list, dft_ener_list


def obtain_gaps(spins_list, excit_ener_list):
    """
    Obtain the S0 -> T1 and S0 -> S1 gap.
    :param: spin_list, excit_energy_list
    :return: gap_t_1, gap_s_1
    """
    check_s1 = 0
    check_t1 = 0
    gap_t_1 = 0
    gap_s_1 = 0

    for k in range(1, len(spins_list)):
        if spins_list[0] == 0:  # If singlet is the ground state
            if spins_list[k] == 0 and check_s1 == 0:
                gap_s_1 = np.round(excit_ener_list[k], 3)  # S0 -> T1 gap
                check_s1 += 1
            elif spins_list[k] == 2 and check_t1 == 0:
                gap_t_1 = np.round(excit_ener_list[k], 3)  # S0 -> S1 gap
                check_t1 += 1

        elif spins_list[0] == 2:  # If triplet is the ground state
            if spins_list[k] == 0 and check_s1 == 0:
                gap_t_1 = (-1) * np.round(excit_ener_list[k], 3)  # S0 -> T1 gap
                check_s1 += 1
            elif spins_list[k] == 0 and check_s1 == 1:
                gap_s1_t1 = (-1) * np.round(excit_ener_list[k], 3)  # S1 -> T1 gap
                gap_s_1 = gap_t_1 - gap_s1_t1  # S0 -> S1 gap
                check_s1 += 1
    return gap_t_1, gap_s_1


def ordering_by_actspace(gaps, order):
    """
    Order the obtained gaps by active spaces in lines.
    :param: gaps
    :return: presentation_list
    """
    gaps = sorted(gaps)
    presentation_list = []

    for ii in range(0, len(order)):
        lista = []
        for jj in range(0, len(gaps)):
            if order[ii] in gaps[jj][0]:
                lista.extend(gaps[jj][1:])
        presentation_list.append(lista)

    for ii in range(0, len(order)):
        presentation_list[ii].insert(0, order[ii])
    return presentation_list


def plots(y_matrix, x_title, y_title, main_title):
    fig, ax = plt.subplots()

    # Features to select:
    fuente = 'sans-serif'  # 'serif'
    small_size = 18
    # medium_size = 22
    bigger_size = 18

    weight_selected = 'normal'
    line_width = 2
    marker_size = 10

    # y_min = 0
    # y_max = 0.7

    # Lines:
    ax.plot(y_matrix[:, 0], y_matrix[:, 3], 'r',
            label='bia16', linewidth=line_width)
    ax.plot(y_matrix[:, 0], y_matrix[:, 6], 'b',
            label='bia17', linewidth=line_width)
    ax.plot(y_matrix[:, 0], y_matrix[:, 9], 'k',
            label='bia18', linewidth=line_width)
    ax.plot(y_matrix[:, 0], y_matrix[:, 12], 'g',
            label='bia19', linewidth=line_width)

    # Markers: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    ax.plot(y_matrix[:, 0], y_matrix[:, 3], 'ro', markersize=marker_size)
    ax.plot(y_matrix[:, 0], y_matrix[:, 6], 'bo', markersize=marker_size)
    ax.plot(y_matrix[:, 0], y_matrix[:, 9], 'ko', markersize=marker_size)
    ax.plot(y_matrix[:, 0], y_matrix[:, 12], 'go', markersize=marker_size)

    # changing the fontsize of yticks
    plt.xticks(fontsize=small_size, weight=weight_selected)
    plt.yticks(fontsize=small_size, weight=weight_selected)

    # Labels:
    # labelpad: change the space between axis umbers and labels
    plt.xlabel(x_title, fontsize=bigger_size, fontfamily=fuente, labelpad=15,
               weight=weight_selected)
    plt.ylabel(y_title, fontsize=bigger_size, fontfamily=fuente, style='italic',
               weight=weight_selected, labelpad=15)
    # plt.ylim([y_min, y_max])  # Limit axis values

    # Title:
    # y=1.05 change the space between title and plot
    plt.title(main_title, fontsize=bigger_size, fontfamily=fuente, y=1.05)

    # Legend
    # legend = plt.legend(fontsize=small_size, fancybox=True, framealpha=0.5,
    #                     labelcolor='linecolor', loc='center right')
    legend = plt.legend(fontsize=small_size, fancybox=True, framealpha=0.5,
                        labelcolor='linecolor')
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    # Frame of the plot: https://e2eml.school/matplotlib_framing.html#spinestyle
    line_width = line_width - 0.8
    ax.spines["top"].set_linewidth(line_width)
    ax.spines["bottom"].set_linewidth(line_width)
    ax.spines["left"].set_linewidth(line_width)
    ax.spines["right"].set_linewidth(line_width)

    plt.show()
    plt.close()