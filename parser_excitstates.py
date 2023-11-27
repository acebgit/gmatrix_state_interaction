"""
 Analysis of the excited states_selected. Printing orbitals involved in the configurations
 with the highest amplitudes in the excited states_selected
"""
__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np
import pandas as pd

from parser_plots import get_bar_chart
from parser_gtensor import get_number_of_states, get_symmetry_states, get_selected_states, get_eigenenergies, \
    take_selected_states_values, get_spin_orbit_couplings, from_energies_soc_to_g_values


def s2_from_file(qchem_file, states_selected):
    """
    get s2 of each state from Q-Chem otuput
    :param: file_ms_notnull
    :return: s2
    """
    search = ['  <S^2>      : ']
    elements = []

    with open(qchem_file, encoding="utf8") as qchem_file:
        for line in qchem_file:
            if any(ii in line for ii in search):
                line = line.split()
                element = float(line[2])
                elements.append(np.round(element, 2))

    s2_selected = take_selected_states_values(elements, states_selected)
    s2_each_states = np.array(s2_selected, dtype=float)
    return s2_each_states


def get_hole_part_contributions(file, totalstates, states_selected):
    """
    Take the hole and particle contributions of each state.
    :param: file_ms_notnull, states_option
    :return: hole_contributions, part_contributions
    """
    with open(file, encoding="utf-8") as file:
        data = file.readlines()

    word_search = ['   Hole: ']
    elements = []

    for line in data:
        if any(i in line for i in word_search):
            element = line[39:]
            elements.append(element.split())
        if len(elements) == totalstates:
            break

    hole_elements = take_selected_states_values(elements, states_selected)
    hole_contributions = np.array(hole_elements, dtype=float)

    word_search = ['   Part: ']
    elements = []

    for line in data:
        if any(i in line for i in word_search):
            element = line[39:]
            elements.append(element.split())
        if len(elements) == totalstates:
            break

    part_elements = take_selected_states_values(elements, states_selected)
    part_contributions = np.array(part_elements, dtype=float)
    return hole_contributions, part_contributions


def get_mulliken_spin(file, totalstates, states_selected):
    """
    Get Mulliken charge and spin of the first atom, i.e. a metal in transition atom complexes.
    :param: file_ms_notnull, states_option, states_selected
    :return: charge_mulliken, spin_mulliken
    """
    word_search = '    Mulliken population analysis '
    element_charge = []
    elements_spin = []

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                for i in range(0, 4):  # Take the 5th line
                    next_line = next(file)
                    i += 1

                next_line = next_line.split()

                element = next_line[2]
                element_charge.append(element)

                element = next_line[3]
                elements_spin.append(element)
            if len(element_charge) == totalstates:
                break

    charge_mulliken_selected = take_selected_states_values(element_charge, states_selected)
    spin_mulliken_selected = take_selected_states_values(elements_spin, states_selected)

    charge_mulliken = np.array(charge_mulliken_selected, dtype=float)
    spin_mulliken = np.array(spin_mulliken_selected, dtype=float)
    return charge_mulliken, spin_mulliken


def get_groundst_socc_values(file, totalstates, states_selected):
    """
    Get spin-orbit coupling constant between states_selected
    :param: file_ms_notnull, states_option
    :return: socc
    """
    with open(file, encoding="utf-8") as file:
        data = file.readlines()

    word_search = ['Mean-Field SOCC']
    all_soccs = []
    n_states = 0

    # Search all SOCCs
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[3]
            all_soccs.append(element)
            n_states += 1

    # Put zero's between same states SOCCS
    for i in range(0, totalstates):
        position = i * totalstates + i
        all_soccs.insert(position, '0.000000000')

    # Take SOCCs between the considered ground state and all other states.
    ground_state_soccs = []
    initial_socc = (states_selected[0]-1) * totalstates
    final_socc = (states_selected[0]-1) * totalstates + totalstates
    for i in range(initial_socc, final_socc):
        ground_state_soccs.append(all_soccs[i])

    # Take SOCCs between the considered ground state and the selected states.
    socc_list = take_selected_states_values(ground_state_soccs, states_selected)
    socc = np.array(socc_list, dtype=float)
    return socc


def get_groundst_orbital_momentum(file, totalstates, states_selected):
    """
    Obtaining the orbitals angular momentum between the ground state and all excited states_selected.
    :param: file_ms_notnull, states_option
    :return: orbital_momentum
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    # Search all orbitals angular momentums
    word_search = ['< B | Lx | A >', '< B | Ly | A >', '< B | Lz | A >']
    all_orbit = []
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[8]
            element = element.replace('i', '')
            all_orbit.append(element)

    # Put zero's between same states
    for i in range(0, totalstates):
        position = (i * totalstates + i) * 3
        all_orbit.insert(position, '0.000000000')
        all_orbit.insert(position, '0.000000000')
        all_orbit.insert(position, '0.000000000')

    # Take momentums between the considered ground state and all other states.
    ground_state_orbit = []
    initial_pos = 3 * ((states_selected[0]-1) * totalstates)
    final_pos = 3 * ((states_selected[0]-1) * totalstates + totalstates)
    for i in range(initial_pos, final_pos):
        ground_state_orbit.append(all_orbit[i])

    # Take momentums between the considered ground state and the selected states.
    all_moment = []
    max_moment = []
    nstate = 0
    for i in range(0, len(ground_state_orbit), 3):
        # Select the direction of momentum that is higher in each state
        x_orb = abs(float(ground_state_orbit[i]))
        y_orb = abs(float(ground_state_orbit[i+1]))
        z_orb = abs(float(ground_state_orbit[i+2]))

        comparison = [x_orb, y_orb, z_orb]
        index_max_momentum = comparison.index(max(comparison))
        max_moment.append(float(ground_state_orbit[nstate * 3 + index_max_momentum]))
        nstate += 1

        # Select all the orbital angular momentums
        all_moment.append([float(ground_state_orbit[i]), float(ground_state_orbit[i+1]),
                           float(ground_state_orbit[i+2])])

    lista = take_selected_states_values(max_moment, states_selected)
    max_moment_array = np.array(lista, dtype=float)
    lista = take_selected_states_values(all_moment, states_selected)
    all_moment_array = np.array(lista, dtype=float)

    return max_moment_array, all_moment_array


def get_ras_spaces(qchem_file):
    """
    Get the active space selected in input and the number of orbitals in RAS1.
    :param: qchem_file
    :return: act_orbitals_array, ras_occ
    """
    with open(qchem_file, encoding="utf-8") as file:
        data = file.readlines()

    ras_act = 0
    word_search = ['RAS_ACT  ']
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            ras_act = int(line[1])

    # List with active orbitals
    word_search = ['RAS_ACT_ORB']
    elements = []
    for line in data:
        if any(i in line for i in word_search):
            line = line.replace('[', '')
            line = line.replace(']', '')
            line = line.replace(',', ' ')
            line = line.split()

            try:
                prueba = float(line[1])
            except ValueError:
                raise ValueError("'RAS_ACT_ORB' has not been manually selected.")

            for i in range(1, ras_act+1):
                elements.append(line[i])
            break
    act_orbitals_array = np.array(elements, dtype=int)

    # Number of occupied orbitals
    word_search = ['RAS_OCC']
    ras_occ = 0
    for line in data:
        if any(i in line for i in word_search):
            split_line = line.split()
            ras_occ = int(split_line[1])
            break
    return act_orbitals_array, ras_occ


def get_alpha_beta(qchem_file):
    """
    Get the alpha and beta electrons.
    :param: qchem_file
    :return: alpha, beta
    """
    with open(qchem_file, encoding="utf-8") as file:
        data = file.readlines()

    # Number of alpha and beta electrons
    word_search = ['beta electrons']
    alpha = 0
    beta = 0
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            alpha = int(line[2])
            beta = int(line[5])
            break

    return alpha, beta


def get_highest_amplitudes(file, amplitude_cutoff):
    """
    Get:
    - index_max_amplitudes: list of the indexes of the configurations with relevant amplitudes (higher than a cut-off)
    - configurations_orbitals: "| HOLE  | ALPHA | BETA  | PART" information of the configurations
    with an amplitude higher than a cutoff
    :param: file_ms_notnull
    :return: index_max_amplitudes, configurations_orbitals
    """

    def get_configurations_information(line):
        """
        From highest to lowest amplitudes, get a list with the amplitudes ('amplitude') and
        the orbitals ('orbitals_information') of all the configurations.
        :param: line
        :return: amplitude, orbitals_information
        """
        elements_amplitude = []
        orbitals_information = []

        while line.startswith(" |"):  # All lines with the configurations start with " |"
            line = line.replace('|       |', '-1')  # "-1" to identify where there is no hole-part orbitals
            line = line.replace('|', '')

            element = line.split()[-1]
            elements_amplitude.append(element)  # Take the amplitudes value

            element = line.split()[0:4]
            orbitals_information.append(element)  # Take the "HOLE  | ALPHA | BETA  | PART" sections

            line = next(file)  # Go to next line

        # Order the configurations (amplitudes and orbitals) from highest to lowest amplitudes
        amplitude = np.array(elements_amplitude, dtype=float)
        amplitude = abs(amplitude)

        for ii in range(0, len(amplitude)):
            for jj in range(ii, len(amplitude)):
                if amplitude[jj] > amplitude[ii]:
                    element = amplitude[jj]
                    amplitude[jj] = amplitude[ii]
                    amplitude[ii] = element

                    change_list = orbitals_information[jj]
                    orbitals_information[jj] = orbitals_information[ii]
                    orbitals_information[ii] = change_list
        return amplitude, orbitals_information

    def other_important_orbitals(amplitude, amplitude_cut_off):
        """
        Get indexes of configurations that have the amplitude higher than a cut-off times the
        amplitude of the first configuration (that is the one with the highest amplitude).
        :param: amplitude, amplitude_cutoff
        :return: indexes
        """
        cut_off = amplitude_cut_off * amplitude[0]

        indexes_list = ['0']
        for k in range(1, len(amplitude)):  # Value in 0 is ignored since amplitudes are sorted
            if amplitude[k] > cut_off:
                # If amplitude is higher than cut-off percentage of first configuration, include it
                indexes_list.append(k)
            else:
                break

        indexes = np.array(indexes_list, dtype=int)
        return indexes

    next_line = 0
    for i in range(0, 2):  # Go to the line of the first configuration
        next_line = next(file)
        i += 1

    amplitudes, configurations_orbitals = get_configurations_information(next_line)
    index_max_amplitudes = other_important_orbitals(amplitudes, amplitude_cutoff)
    return index_max_amplitudes, configurations_orbitals


def get_orbital(homo_orbital, configuration_data, initial_active_orbitals):
    """
    Get orbital_list with unpaired electrons in the hole, active or particle configurations
    :param: homo_orbital, configuration_data, initial_active_orbitals
    :return: orbitals (string)
     """
    def from_ras_to_scf_order(homo_orbit, initial_scf_space, ras_orbitals):
        """
        Change from RAS order to SCF order.
        Example: if in (1,2,3,4) orbital_list the active space selected is [1,3],
        RAS_order is (1,2,3,4) while SCF_order is (2,1,3,4), since orbital_list are
        order by RAS1-RAS2-RAS3.
        :param: homo_orbit, initial_scf_space, ras_orbitals
        :return: final_SCF_space
        """
        initial_scf_space = initial_scf_space.tolist()
        scf_orbitals_order = []

        # scf_orbitals_order: orbital_list in HF energetic order
        for orbital in range(1, homo_orbit + 1):  # RAS1
            if orbital not in initial_scf_space:
                scf_orbitals_order.append(orbital)

        for orbital in initial_scf_space:  # RAS2
            scf_orbitals_order.append(orbital)

        for orbital in range(homo_orbit + 1, homo_orbit + 50):  # RAS3
            if orbital not in initial_scf_space:
                scf_orbitals_order.append(orbital)

        # ras_orbitals_order: orbital_list in RAS energetic order, meaning RAS1 - RAS2 - RAS3
        ras_orbitals_order = list(range(1, homo_orbit + 50))

        # print('ras_orbitals_order', 'scf_orbitals_order')
        # for i in range(0, len(ras_orbitals_order)):
        #     print(ras_orbitals_order[i], '--', scf_orbitals_order[i])
        # exit()

        indx = ras_orbitals_order.index(ras_orbitals)
        final_scf_orbital = scf_orbitals_order[indx]
        return final_scf_orbital

    orbital_list = []

    # If orbitals is in hole
    if configuration_data[0] != '-1':
        ras_orbital = int(configuration_data[0])
        new_orbital = from_ras_to_scf_order(homo_orbital, initial_active_orbitals, ras_orbital)
        new_orbital = int(new_orbital)
        orbital_list.append(new_orbital)

    # If orbitals is in active space
    alpha_config = str(configuration_data[1])
    beta_config = str(configuration_data[2])

    for active_orbital_index in range(0, len(alpha_config)):
        if alpha_config[active_orbital_index] != beta_config[active_orbital_index]:
            new_orbital = initial_active_orbitals[active_orbital_index]
            new_orbital = int(new_orbital)
            orbital_list.append(new_orbital)

    # If orbitals is in particle
    if configuration_data[3] != '-1':
        ras_orbital = int(configuration_data[3])
        new_orbital = from_ras_to_scf_order(homo_orbital, initial_active_orbitals, ras_orbital)
        new_orbital = int(new_orbital)
        orbital_list.append(new_orbital)

    if not orbital_list:  # If list is empty, append a 0
        orbital_list.append('-')

    # Put orbitals as a string
    if len(orbital_list) > 1:
        orbitals = ','.join([str(n) for n in orbital_list])
    else:
        orbitals = orbital_list[0]
    return orbitals


def get_configurations_unpaired_orbitals(ras_input, states_selected, cutoff):
    """
    Obtain the orbitals with the unpaired electrons in each relevant configuration of each state.
    :return:
    """
    initial_active_orbitals, ras_occ = get_ras_spaces(ras_input)
    elec_alpha, elec_beta = get_alpha_beta(ras_input)

    word_search = ' | HOLE  | '
    orbit_list = []
    n_state = 0

    with open(ras_input, encoding="utf-8") as ras_input:
        for line in ras_input:
            if word_search in line:  # Go to configurations line

                index_relevant_amplit, state_orbitals = get_highest_amplitudes(ras_input, cutoff)
                n_configuration = 0

                # Include all those configurations with relevant amplitudes in the final list
                for i in index_relevant_amplit:
                    new_orbitals = get_orbital(elec_beta, state_orbitals[i], initial_active_orbitals)
                    orbit_dict = {'State': n_state+1, 'Configuration': n_configuration+1, 'SOMO orbitals': new_orbitals}
                    orbit_list.append(orbit_dict)

                    n_configuration += 1
                n_state += 1

    orbit_list_selected = []
    for state in states_selected:
        for i in range(0, len(orbit_list)):
            if orbit_list[i]["State"] == state:
                orbit_list_selected.append(orbit_list[i])
    initial_active_orbitals_list = list(initial_active_orbitals)
    return orbit_list_selected, initial_active_orbitals_list

def count_singlet_triplets(states, s2_lists):
        singlets = []
        doublets = []
        triplets = []
        quartets = []
        quintets = []
        for ii in range(0, len(states)):
            if s2_lists[ii] == 0:
                singlets.append(ii+1)
            elif s2_lists[ii] == 0.75:
                doublets.append(ii+1)
            elif s2_lists[ii] == 2:
                triplets.append(ii+1)
            elif s2_lists[ii] == 3.75:
                quartets.append(ii+1)
            elif s2_lists[ii] == 6:
                quintets.append(ii + 1)
        print(len(singlets), 'Singlet states: ', singlets)
        print(len(doublets), 'Doublet states: ', doublets)
        print(len(triplets), 'Triplet states: ', triplets)
        print(len(quartets), 'Quartet states: ', quartets)
        print(len(quintets), 'Quintet states: ', quintets)

def add_orbitals_active_space(orbitals, active_space):
    """
    Return the orbitals to be added to the active space list
    :param orbitals:
    :param active_space:
    :return: active space
    """
    if type(orbitals) == int:  # It is a doublet, add SOMO
        active_space.append(int(orbitals))

    elif type(orbitals) == str:  # It is higher multiplicities, add all SOMOs
        orbitals = orbitals.split(',')
        for j in range(0, len(orbitals)):
            active_space.append(int(orbitals[j]))
    return active_space

def get_excited_states_analysis(file, state_selections, states_ras, symmetry_selected, cut_off, cut_soc,
                                cut_ang, cut_gvalue, ppms, plots, save_pict):
    """
    Obtaining a matrix with several data for each excited state. The cut-off determines the fraction of the amplitude
    of the 1st configuration that need to have the other configurations to be shown in each state.
    :param: file_ms_notnull, cutoff
    :return: excited_states_presentation_matrix
    """
    def compare_gvalues_and_estimation(excited_states_presentation_matrix):
        """
        Print a list with 3 columns: state, g-shift and the relation L*SOCC/(E**2)
        :param excited_states_presentation_matrix
        :return:
        """
        largest_gvalues = []

        if cut_gvalue != 0:
            for i in range(2, len(excited_states_presentation_matrix)):
                # Take largest g-shift and put it in a list with each state
                gxx = abs(excited_states_presentation_matrix[i,10] - excited_states_presentation_matrix[i-1,10])
                gyy = abs(excited_states_presentation_matrix[i,11] - excited_states_presentation_matrix[i-1,11])
                gzz = abs(excited_states_presentation_matrix[i,12] - excited_states_presentation_matrix[i-1,12])
                compare_list = [gxx, gyy, gzz]
                compare_list.sort(reverse=True)

                # Take largest relation and put it in a list with each state
                ener_compare = excited_states_presentation_matrix[i,5]
                socc_compare = excited_states_presentation_matrix[i,7]
                orbit_compare = excited_states_presentation_matrix[i,8]
                relation = abs(orbit_compare * socc_compare / ener_compare**2)

                largest_gvalues.append([i, compare_list[0], relation])

            df = pd.DataFrame(largest_gvalues, columns=['State', 'g-shift', 'L*SOCC/E**2'])
            order_states_gvalues = df.sort_values('g-shift')
            print(order_states_gvalues)

    def gshift_calculation_loop(states_ras, file, totalstates, ppms):
        """
        Calculate the g-values in all the selected states.
        :return:
        """
        for i in range(1, len(states_ras) + 1):
            states_sos = states_ras[0:i]
            print('Calculating....    states', states_sos)
            eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_sos)
            selected_socs, sz_list, ground_sz = get_spin_orbit_couplings(file, totalstates, states_sos, soc_option=0)

            g_shift = from_energies_soc_to_g_values(file, states_sos,
                                                    totalstates, excitation_energies_ras,
                                                    selected_socs, sz_list, ground_sz, ppms)
            gxx_list.append(g_shift[0].real)
            gyy_list.append(g_shift[1].real)
            gzz_list.append(g_shift[2].real)
        return gxx_list, gyy_list, gzz_list

    def gshift_estimation_loop(nstates, orbit_moments, soccs, energies):
        """
        Calculate the g-values in all the selected states.
        :return:
        """
        g_xx = [0]
        g_yy = [0]
        g_zz = [0]

        for i in range(1, len(nstates)):
            element_xx = -4 * orbit_moments[i, 0] * soccs[i] / energies[i]
            element_yy = -4 * orbit_moments[i, 1] * soccs[i] / energies[i]
            element_zz = -4 * orbit_moments[i, 2] * soccs[i] / energies[i]

            g_xx.append(element_xx)
            g_yy.append(element_yy)
            g_zz.append(element_zz)
        return g_xx, g_yy, g_zz

    file_string = file

    totalstates = get_number_of_states(file)

    state_symmetries, ordered_state_symmetries = get_symmetry_states(file, totalstates)

    states_ras = get_selected_states(file, totalstates, states_ras, state_selections, symmetry_selected)

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
    excitation_energies_ras[:] = (excitation_energies_ras[:] - excitation_energies_ras[0]) * 27.211399

    s2_list = s2_from_file(file, states_ras)

    hole_contributions, part_contributions = get_hole_part_contributions(file, totalstates, states_ras)

    socc_values = get_groundst_socc_values(file, totalstates, states_ras)

    orbitmoment_max, orbitmoment_all = get_groundst_orbital_momentum(file, totalstates, states_ras)

    elec_alpha, elec_beta = get_alpha_beta(file)

    configuration_orbitals, initial_active_orbitals = get_configurations_unpaired_orbitals(file, states_ras, cut_off)

    if cut_gvalue != 0:
        excited_states_presentation_list = [['State', 'Config.', 'Sym.', 'Hole', 'Part',
                                             'ΔE (eV)', 'Unpaired orb.', 'SOCC',
                                             'Máx Lk', 'S^2', 'gxx', 'gyy', 'gzz']]
        gxx_list = []
        gyy_list = []
        gzz_list = []

        # gxx_list, gyy_list, gzz_list = gshift_calculation_loop(states_ras, file, totalstates, ppms)
        gxx_list, gyy_list, gzz_list = gshift_estimation_loop(states_ras, orbitmoment_all, socc_values, excitation_energies_ras)

    else:
        excited_states_presentation_list = [['State', 'Config.', 'Sym.', 'Hole', 'Part',
                                         'ΔE (eV)', 'Unpaired orb.', 'SOCC',
                                         'Máx Lk', 'S^2']]

    final_active_orbitals = []
    for i in range(0, len(configuration_orbitals)):
        state = configuration_orbitals[i]["State"]
        configuration = configuration_orbitals[i]["Configuration"]

        state_index = states_ras.index(state)
        symmetry = ordered_state_symmetries[state-1]
        hole = np.around(float(hole_contributions[state_index]), 2)
        part = np.around(float(part_contributions[state_index]), 2)
        excit_energy = np.round(float(excitation_energies_ras[state_index]), 3)
        soc = np.round(float(socc_values[state_index]), 3)
        orbital_ground_state = np.round(float(orbitmoment_max[state_index]), 3)
        s2 = s2_list[state_index]

        orbital = configuration_orbitals[i]["SOMO orbitals"]
        if orbital == '-':  # It is a singlet, there is no SOMO so add HOMO
            final_active_orbitals.append(elec_alpha)

        if cut_gvalue != 0:
            g_xx = np.round(gxx_list[state_index], 10)
            g_yy = np.round(gyy_list[state_index], 10)
            g_zz = np.round(gzz_list[state_index], 10)

            g_xx_diff = abs(gxx_list[state_index] - gxx_list[state_index-1])
            g_yy_diff = abs(gyy_list[state_index] - gyy_list[state_index-1])
            g_zz_diff = abs(gzz_list[state_index] - gzz_list[state_index-1])

            if g_xx_diff >= cut_gvalue or g_yy_diff >= cut_gvalue or g_zz_diff >= cut_gvalue:
                excited_states_presentation_list.append([state, configuration, symmetry,
                                                         hole, part, excit_energy, orbital, soc,
                                                         orbital_ground_state, s2, g_xx, g_yy, g_zz])
                final_active_orbitals = add_orbitals_active_space(orbital, final_active_orbitals)

        else:
            if abs(orbital_ground_state) >= cut_ang and abs(soc) >= cut_soc:
                excited_states_presentation_list.append([state, configuration, symmetry,
                                                     hole, part, excit_energy, orbital, soc,
                                                     orbital_ground_state, s2])
                final_active_orbitals = add_orbitals_active_space(orbital, final_active_orbitals)

    excited_states_presentation_matrix = np.array(excited_states_presentation_list, dtype=object)
    # compare_gvalues_and_estimation(excited_states_presentation_matrix)

    orbital_set = set(final_active_orbitals)
    final_active_orbitals = list(orbital_set)
    final_active_orbitals.sort()

    print("------------------------")
    print(" EXCITED STATES ANALYSIS")
    print("------------------------")
    print('Most important settings for each state (amplitude_cutoff:', cut_off, ') :')
    print('\n'.join(''.join('{:^15}'.format(item) for item in row)
                    for row in excited_states_presentation_matrix[:, :]))
    count_singlet_triplets(states_ras, s2_list)
    print(" ")

    print("------------------------")
    print(" IMPROVED ACTIVE SPACE ")
    print("------------------------")

    electrons = get_new_active_space_electrons(initial_active_orbitals, elec_alpha, elec_beta)
    print('Initial active space (HOMO =', elec_alpha, '):', '[', electrons, ',', len(initial_active_orbitals),
          '] ;',
          initial_active_orbitals)

    electrons = get_new_active_space_electrons(final_active_orbitals, elec_alpha, elec_beta)
    print('Final active space (HOMO =', elec_alpha, '):', '[', electrons, ',', len(final_active_orbitals), '] ;',
          final_active_orbitals)

    if plots == 1:
        states_plot = excited_states_presentation_matrix[1:, 0]
        ener_plot = excited_states_presentation_matrix[1:, 5]
        orbit_plot = excited_states_presentation_matrix[1:, 8]
        soc_plot = excited_states_presentation_matrix[1:, 7]

        if cut_gvalue != 0:
            g_plot = []
            for i in range(1, len(excited_states_presentation_matrix)):
                xx_g = excited_states_presentation_matrix[i, 10]
                yy_g = excited_states_presentation_matrix[i, 11]
                zz_g = excited_states_presentation_matrix[i, 12]
                g_comparison = [abs(xx_g), abs(yy_g), abs(zz_g)]
                g_comparison.sort(reverse=True)
                g_plot.append(g_comparison[0])
            get_bar_chart(file_string[:-4], states_plot, g_plot, 'Electronic State',
                      r'$\mathregular{\Delta g_{xx}}$', 'All variables analysis', save_pict)

        get_bar_chart(file_string[:-4], states_plot, ener_plot, 'Electronic State',
                      'Excitation energy (eV)', 'energ_analysis', save_pict)
        get_bar_chart(file_string[:-4], states_plot, orbit_plot, 'Electronic State',
                      'Orbital angular momentum', 'orbitmoment_analysis', save_pict)
        get_bar_chart(file_string[:-4], states_plot, soc_plot, 'Electronic State',
                      'SOCC (cm-1)', 'socc_analysis', save_pict)


def get_new_active_space_electrons(active_space, alpha, beta):
    """
    Electrons in the new active space considering if they are occupied
    or unoccupied observing the HOMO position
     """
    electrons = 0
    for i in range(0, len(active_space)):
        if active_space[i] <= beta:
            electrons += 2
        elif (active_space[i] > beta) and (active_space[i] <= alpha):
            electrons += 1
        elif active_space[i] > alpha:
            pass
    return electrons


def improved_active_space(file, states_option, selected_states, cut_off, cut_soc, cut_ang):
    """
    Obtain an improved active space of RAS-CI Q-Chem output including orbitals
    with unpaired electrons in relevant hole/particle configurations.
    Cut_off is the amplitude at which the second-highest amplitude
    configuration is included.
    :param: file
    :param: states_option, selected_states, cut_off, cut_soc, cut_ang
    """
    totalstates = get_number_of_states(file)

    states_ras = get_selected_states(file, totalstates, selected_states, states_option, symmetry_selection='None')

    elec_alpha, elec_beta = get_alpha_beta(file)

    configuration_orbitals, initial_active_orbitals = get_configurations_unpaired_orbitals(file, states_ras, cut_off)

    socc_values = get_groundst_socc_values(file, totalstates, states_ras)

    orbitmoment_max, orbitmoment_all= get_groundst_orbital_momentum(file, totalstates, states_ras)

    final_active_orbitals = []
    for i in range(0, len(configuration_orbitals)):
        state = configuration_orbitals[i]['State'] - 1
        orbitals = configuration_orbitals[i]['SOMO orbitals']

        if orbitals == '-':  # It is a singlet, there is no SOMO so add HOMO
            final_active_orbitals.append(elec_alpha)

        elif type(orbitals) == int:  # It is a doublet, add SOMO
            if abs(socc_values[state]) >= cut_soc and abs(orbitmoment_max[state]) >= cut_ang:
                final_active_orbitals.append(int(orbitals))

        elif type(orbitals) == str:  # It is higher multiplicities, add all SOMOs
            orbitals = orbitals.split(',')
            for j in range(0, len(orbitals)):
                if abs(socc_values[state]) >= cut_soc and abs(orbitmoment_max[state]) >= cut_ang:
                    final_active_orbitals.append(int(orbitals[j]))

    orbital_set = set(final_active_orbitals)
    final_active_orbitals = list(orbital_set)
    final_active_orbitals.sort()

    print("------------------------")
    print(" IMPROVED ACTIVE SPACE ")
    print("------------------------")

    electrons = get_new_active_space_electrons(initial_active_orbitals, elec_alpha, elec_beta)
    print('Initial active space (HOMO =', elec_alpha, '):', '[', electrons, ',', len(initial_active_orbitals), '] ;',
          initial_active_orbitals)

    electrons = get_new_active_space_electrons(final_active_orbitals, elec_alpha, elec_beta)
    print('Final active space (HOMO =', elec_alpha, '):', '[', electrons, ',', len(final_active_orbitals), '] ;',
          final_active_orbitals)
