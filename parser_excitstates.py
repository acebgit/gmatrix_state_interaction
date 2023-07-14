"""
    EXCITED STATES MOST IMPORTANT CONFIGURATIONS
 Analysis of the excited states and print of those
 orbitals involved in the configurations with the
 highest amplitudes in the excited states
"""
import numpy as np

from parser_init import get_number_of_states, get_eigenenergies, get_selected_states, \
    get_symmetry_states, get_hole_part_contributions, get_socc_values, get_ground_state_orbital_momentum, \
    get_mulliken_spin


def get_ras_spaces(qchem_file):
    """
    Get the active space selected in input and the number of orbitals in RAS1.
    :param: qchem_file
    :return: ras_act_orb, ras_occ
    """
    with open(qchem_file, encoding="utf-8") as file:
        data = file.readlines()

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

            for i in range(1, ras_act+1):
                elements.append(line[i])
            break
    ras_act_orb = np.array(elements, dtype=int)

    # Number of occupied orbitals
    word_search = ['RAS_OCC']
    for line in data:
        if any(i in line for i in word_search):
            split_line = line.split()
            ras_occ = int(split_line[1])
            break

    return ras_act_orb, ras_occ


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
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            alpha = int(line[2])
            beta = int(line[5])
            break

    return alpha, beta


def get_highest_amplitudes(file):
    """
    Get:
    - index_max_amplitudes: list of the indexes of the configurations with relevant amplitudes (higher than a cut-off)
    - configurations_orbitals: "| HOLE  | ALPHA | BETA  | PART" information of the configurations with an amplitude higher
    than a cutoff
    :param: file
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
            line = line.replace('|       |', '-1')  # "-1" to identify where there is no hole-part orbital
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

    def other_important_orbitals(amplitude, amplitude_cutoff):
        """
        Get indexes of configurations that have the amplitude higher than a cut-off times the
        amplitude of the first configuration (that is the one with the highest amplitude).
        :param: amplitude, amplitude_cutoff
        :return: indexes
        """
        cut_off = amplitude_cutoff * amplitude[0]

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
    index_max_amplitudes = other_important_orbitals(amplitudes, amplitude_cutoff=0.7)
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

    # If orbital is in hole
    if (configuration_data[0] != '-1'):
        ras_orbital = int(configuration_data[0])
        new_orbital = from_ras_to_scf_order(homo_orbital, initial_active_orbitals, ras_orbital)
        new_orbital = int(new_orbital)
        orbital_list.append(new_orbital)

    # If orbital is in active space
    alpha_config = str(configuration_data[1])
    beta_config = str(configuration_data[2])

    for active_orbital_index in range(0, len(alpha_config)):
        if alpha_config[active_orbital_index] != beta_config[active_orbital_index]:
            new_orbital = initial_active_orbitals[active_orbital_index]
            new_orbital = int(new_orbital)
            orbital_list.append(new_orbital)

    # If orbital is in particle
    if (configuration_data[3] != '-1'):
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


def print_excited_states(presentation_list, n_states, hole_contributions,
                         part_contributions, socc_values, excitation_energies_ev,
                         state_symmetries, new_orbital, orbital_momentum, mulliken_spin):
    """
    Prepare tHe presentation list with each excited state values
    :param: presentation_list, n_states, hole_contributions,
                         part_contributions, socc_values, excitation_energies_ev,
                         state_symmetries, new_orbital, orbital_momentum, mulliken_spin
    :return: presentation_list, SOC
     """
    state = np.round(int(n_states), 0) + 1
    symmetry = state_symmetries[n_states]

    hole = np.round(float(hole_contributions[n_states]), 2)
    part = np.round(float(part_contributions[n_states]), 2)

    soc = np.round(float(socc_values[n_states]), 0)
    excit_energy = np.round(float(excitation_energies_ev[n_states]), 3)

    orbital_ground_state = np.round(float(orbital_momentum[n_states]), 3)

    spin = np.round(float(mulliken_spin[n_states]), 3)

    presentation_list.append([state, symmetry, hole, part, excit_energy, new_orbital, soc, orbital_ground_state, spin])
    return presentation_list, soc


def get_new_active_space_electrons(new_active_space, homo_orbital):
    """
    Electrons in the new active space considering if they are occupied
    or unoccupied observing the HOMO position
     """
    electrons = 0
    for i in range(0, len(new_active_space)):
        if new_active_space[i] <= homo_orbital:
            electrons += 2
        else:
            pass

    if electrons <= 0:
        electrons = 0
    else:
        electrons = electrons - 1

    return electrons


def get_excited_states_analysis(file):
    """
    Obtention of a matrix with several data for each excited state. Cut-off determines the amplitude
    difference between the 1st and 2nd (or more) important configurations in each state: if it is lower,
    the state will appear another time with the 2nd configuration.
    :param: input
    :return: excited_states_presentation_matrix
    """
    totalstates = get_number_of_states(file)

    state_symmetries, ordered_state_symmetries = get_symmetry_states(file, totalstates)

    states_ras = get_selected_states(file, totalstates, selected_states=0,
                                     states_option=1, symmetry_selection='None')

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)

    hole_contributions, part_contributions = get_hole_part_contributions(file, totalstates)

    mulliken_charge, mulliken_spin = get_mulliken_spin(file, totalstates, states_ras)

    socc_values = get_socc_values(file, totalstates)

    orbital_momentum = get_ground_state_orbital_momentum(file, totalstates)

    initial_active_orbitals, ras_occ = get_ras_spaces(file)

    alpha, homo_orbital = get_alpha_beta(file)

    excited_states_presentation_list = [['State', 'Symmetry', 'Hole', 'Part',
                                         'Excitation energy (eV)', 'Orbitals', 'SOCC (cm-1)',
                                         'Orbital momentum', 'Mulliken Spin']]

    word_search = ' | HOLE  | '
    n_states = 0

    with open(file, encoding="utf-8") as file:
        for line in file:
            if word_search in line:  # Go to configurations line

                index_relevant_amplit, state_orbitals = get_highest_amplitudes(file)

                # Include all those configurations with relevant amplitudes in the final list
                for i in index_relevant_amplit:
                    new_orbitals = get_orbital(homo_orbital, state_orbitals[i], initial_active_orbitals)

                    excited_states_presentation_list, soc = print_excited_states(excited_states_presentation_list,
                                                                                 n_states, hole_contributions,
                                                                                 part_contributions, socc_values,
                                                                                 excitation_energies_ras * 27.211399,
                                                                                 ordered_state_symmetries, new_orbitals,
                                                                                 orbital_momentum, mulliken_spin)

                n_states += 1

    excited_states_presentation_matrix = np.array(excited_states_presentation_list, dtype=object)
    print("------------------------")
    print(" EXCITED STATES ANALYSIS")
    print("------------------------")
    print('Most important settings for each state (amplitude_cutoff: 0.3) :')
    print('\n'.join(''.join('{:^20}'.format(item) for item in row)
                    for row in excited_states_presentation_matrix[:, :]))
    print(" ")


def improved_active_space(file):
    from parser_init import get_number_of_states, get_eigenenergies, get_selected_states, get_symmetry_states, \
        get_hole_part_contributions, get_socc_values

    totalstates = get_number_of_states(file)

    state_symmetries, ordered_state_symmetries = get_symmetry_states(file, totalstates)

    states_ras = 0
    selected_states = 1
    states_ras = get_selected_states(file, totalstates, states_ras, selected_states, symmetry_selection='None')

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)

    hole_contributions, part_contributions = get_hole_part_contributions(file, totalstates)

    orbital_momentum = get_ground_state_orbital_momentum(file, totalstates)

    socc_values = get_socc_values(file, totalstates)

    mulliken_charge, mulliken_spin = get_mulliken_spin(file, totalstates, states_ras)

    initial_active_orbitals, ras_occ = get_ras_spaces(file)

    alpha, homo_orbital = get_alpha_beta(file)

    # TAKE STATES WITH HIGHEST CONTRIBUTION
    new_space_list = []
    final_scf_ordered_space = []
    excited_states_presentation_list = ['State', 'Symmetry', 'Hole', 'Part',
                                        'Excitation energy (eV)', 'Orbitals', 'SOCC (cm-1)',
                                        'Orbital momentum']

    word_search = ' | HOLE  | '
    n_states = 0

    with open(file, encoding="utf-8") as file:
        for line in file:
            if word_search in line:  # Go to configurations line

                index_max_amplitudes, state_orbitals = get_highest_amplitudes(file)

                for i in index_max_amplitudes:
                    new_orbital = get_orbital(homo_orbital, state_orbitals[i], initial_active_orbitals)

                    if new_orbital not in new_space_list:
                        new_space_list.append(new_orbital)
                        excited_states_presentation_list, soc = print_excited_states(excited_states_presentation_list,
                                                                                     n_states, hole_contributions,
                                                                                     part_contributions, socc_values,
                                                                                     excitation_energies_ras
                                                                                     * 27.211399,
                                                                                     state_symmetries, new_orbital,
                                                                                     orbital_momentum, mulliken_spin)

                        if soc != 0 or i == 0:
                            final_scf_ordered_space.append(new_orbital)

                n_states += 1

    print("------------------------")
    print(" IMPROVED ACTIVE SPACE ")
    print("------------------------")

    print('Initial active space (HOMO =', homo_orbital, '):')
    initial_active_orbitals_list = initial_active_orbitals.tolist()
    electrons = get_new_active_space_electrons(initial_active_orbitals_list, homo_orbital)
    print('[', electrons, ',', len(initial_active_orbitals_list), '] ;', initial_active_orbitals_list)

    # print('')
    # new_active_space = np.array(new_space_list, dtype=int)
    # print('New active space (SOCC not zero, HOMO singly occupied):')
    # electrons = get_new_active_space_electrons(final_SCF_ordered_space, homo_orbital)
    # final_SCF_ordered_space.sort()
    # print('[',electrons, ',', len(final_SCF_ordered_space), '] ;', final_SCF_ordered_space)

    print('')
    print('Initial active space + New active space:')

    initial = set(initial_active_orbitals)
    final = set(final_scf_ordered_space)
    in_final_but_not_in_initial = final - initial
    initial_plus_final_space = initial_active_orbitals + list(in_final_but_not_in_initial)

    electrons = get_new_active_space_electrons(initial_plus_final_space, homo_orbital)
    initial_plus_final_space.sort()
    print('[', electrons, ',', len(initial_plus_final_space), '] ;', initial_plus_final_space)
