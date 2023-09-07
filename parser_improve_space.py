"""
Obtain an improved active space of RAS-CI Q-Chem output including orbitals
with unpaired electrons in relevant hole/particle configurations.
"""
import numpy as np
import sys


def get_number_of_states(file):
    """
    Obtain the total number of states in ras
    :param: file
    :return: nstates
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    element = 0
    word_search = ['Requested states: ']

    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[3]
            break
    totalstates = int(element)
    return totalstates


def get_symmetry_states(file, totalstates):
    """
    Create two lists: first with the symmetry of each state (A1,A2,A3...) and second with
    the order of these symmetries (1A1,2A2,3A3...)
    :param: file_ras, nstates
    :return: all_state_symmetries, ordered_state_symmetries
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    word_search = ['State symmetry: ']
    elements = []
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[2]
            element = element.replace('*', '')  # '*' means mix of symmetries
            elements.append(element)

        if len(elements) == totalstates:
            break  # All symmetries appended
    all_state_symmetries = np.array(elements)

    ordered_state_symmetries = []  # State symmetries with the order (1,2,3..)
    for i in range(0, len(all_state_symmetries)):
        number = 1
        for j in range(0, i):
            if all_state_symmetries[i] == all_state_symmetries[j]:
                number += 1
        element = str(number) + all_state_symmetries[i]
        ordered_state_symmetries.append(element)
    return all_state_symmetries, ordered_state_symmetries


def get_selected_states(file, totalstates, symmetry_selection):
    """
    Select the states used depending on "states_option" value:
    0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
    :param: file, nstates, selected_states, states_option, symmetry_selection
    :return: selected_states
    """
    states_option = 1
    selected_states = [0]
    all_symmetries, ordered_symmetries = get_symmetry_states(file, totalstates)

    if states_option == 0:  # Se
        for i in selected_states:
            if i <= 0 or i > totalstates:
                print("The number of states selected must be among the total number of states calculated in QChem.")
                print("Select a different number of states")
                sys.exit()

    elif states_option == 1:
        selected_states = list(range(1, totalstates + 1))

    elif states_option == 2:
        states_selected_by_symmetry = [1]  # GS is always added

        for nstate in range(1, totalstates + 1):
            if (all_symmetries[nstate - 1] == symmetry_selection) and (nstate != 1):
                states_selected_by_symmetry.append(nstate)
        selected_states = states_selected_by_symmetry

        if states_selected_by_symmetry == [1]:
            print('There is not this symmetry.')
            print('Change the symmetry selection:')
            for nstate in range(0, totalstates):
                print('- State', nstate + 1, ', symmetry', all_symmetries[nstate])
            sys.exit()
    return selected_states


def get_eigenenergies(file, totalstates, selected_states):
    """
    Get energies of the RAS-CI states.
    :param: file, nstates, selected_states
    :return: eigenenergies, excitation_energies
    """
    word_search = ' RAS-CI total energy for state  '
    element_energy = []
    elements_excitenergy = []

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                line = line.split()
                element = line[6]
                element_energy.append(element)

                next_line = next(file)
                next_line = next_line.split()
                element = next_line[4]
                elements_excitenergy.append(element)
            if len(element_energy) == totalstates:
                break

    energies_selected = []
    excited_energies_selected = []
    for i in selected_states:
        energies_selected.append(element_energy[i - 1])
        excited_energies_selected.append(elements_excitenergy[i - 1])

    eigenenergies = np.array(energies_selected, dtype=float)
    excitation_energies_ev = np.array(excited_energies_selected, dtype=float)
    excitation_energies = excitation_energies_ev / 27.211399  # From eV to a.u.
    return eigenenergies, excitation_energies


def get_hole_part_contributions(file, totalstates):
    """
    Take the hole and particle contributions of each state.
    :param: file, nstates
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

    hole_contributions = np.array(elements, dtype=float)

    word_search = ['   Part: ']
    elements = []

    for line in data:
        if any(i in line for i in word_search):
            element = line[39:]
            elements.append(element.split())
        if len(elements) == totalstates:
            break
    part_contributions = np.array(elements, dtype=float)
    return hole_contributions, part_contributions


def get_mulliken_spin(file, totalstates, states):
    """
    Get Mulliken charge and spin of the first atom, i.e. a metal in transition atom complexes.
    :param: file, nstates, states
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

    charge_mulliken_selected = []
    spin_mulliken_selected = []

    for i in states:
        charge_mulliken_selected.append(element_charge[i - 1])
        spin_mulliken_selected.append(elements_spin[i - 1])

    charge_mulliken = np.array(charge_mulliken_selected, dtype=float)
    spin_mulliken = np.array(spin_mulliken_selected, dtype=float)
    return charge_mulliken, spin_mulliken


def get_socc_values(file, totalstates):
    """
    Get spin-orbit coupling constant between states
    :param: file, nstates
    :return: socc
    """
    with open(file, encoding="utf-8") as file:
        data = file.readlines()

    word_search = ['Mean-Field SOCC']
    elements = []
    n_states = 0

    elements.append('0.000000')  # SOCC between state 1 and 1 is zero
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[3]
            elements.append(element)
            n_states += 1
        if n_states == totalstates - 1:
            break

    socc = np.array(elements, dtype=float)
    return socc


def get_ground_state_orbital_momentum(file, totalstates):
    """
    Obtaining the orbital angular momentum between the ground state and all excited states.
    :param: file, nstates
    :return: orbital_momentum
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    word_search = ['< B | Lx | A >', '< B | Ly | A >', '< B | Lz | A >']
    elements = ['0.000', '0.000', '0.000']  # L between ground state and it-self does not exist (x,y,z)

    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[8]
            element = element.replace('i', '')
            elements.append(element)

        if len(elements) == totalstates * 3:
            break  # There are all the Lk values with between ground and excited states

    orbital_momentum = []
    nstate = 0

    for i in range(0, len(elements), 3):
        one_state_momentum = []

        x_orb = abs(float(elements[i]))
        one_state_momentum.append(x_orb)

        y_orb = abs(float(elements[i+1]))
        one_state_momentum.append(y_orb)

        z_orb = abs(float(elements[i+2]))
        one_state_momentum.append(z_orb)

        # Select the orbital angular momentum that is higher, to be shown
        index_max_momentum = one_state_momentum.index(max(one_state_momentum))
        orbital_momentum.append(float(elements[nstate*3 + index_max_momentum]))

        nstate += 1
    return orbital_momentum


def get_ras_spaces(qchem_file):
    """
    Get the active space selected in input and the number of orbitals in RAS1.
    :param: qchem_file
    :return: ras_act_orb, ras_occ
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

            for i in range(1, ras_act+1):
                elements.append(line[i])
            break
    ras_act_orb = np.array(elements, dtype=int)

    # Number of occupied orbitals
    word_search = ['RAS_OCC']
    ras_occ = 0
    for line in data:
        if any(i in line for i in word_search):
            split_line = line.split()
            ras_occ = int(split_line[1])
            break

    return ras_act_orb, ras_occ


def s2_from_file(qchem_file):
    """
    get s2 of each state from Q-Chem otuput
    :param: file
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

    s2_each_states = np.array(elements, dtype=float)
    return s2_each_states


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

    # If orbital is in hole
    if configuration_data[0] != '-1':
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


def improved_active_space(file):

    totalstates = get_number_of_states(file)

    state_symmetries, ordered_state_symmetries = get_symmetry_states(file, totalstates)

    states_ras = get_selected_states(file, totalstates, symmetry_selection='None')

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

                index_max_amplitudes, state_orbitals = get_highest_amplitudes(file, cutoff)

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
