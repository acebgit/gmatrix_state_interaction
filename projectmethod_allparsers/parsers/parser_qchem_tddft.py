import json
import numpy as np

state_selection = 1  # 0: use "state_ras" ; 1: use all states_selected
initial_states = [1, 2, 3, 4, 5] # 0 is the ground state
ground_state = '0'
selected_multiplicity = 2

file = '../\
test/qchem_tddft_doublets.out'  # str(sys.argv[1])

#################################
# FUNCTIONS AND CLASSES    ######
#################################
multiplicity_dict = {1: 'Singlet', 2: 'Doublet', 3: 'Triplet', 4: 'Quartet', 5: 'Quintet', 6: 'Sextet'}


def search_word_line(archivo, text):
    """
    Search a text and goes to that line.
    :param archivo
    :param text
    :return: linee
    """
    linee = next(archivo)
    while text not in linee:
        linee = next(archivo)
    return linee


def get_number_of_states(filee):
    """
    Obtain the total number of states selected in TD-DFT.
    :param: filee
    :return: nstates
    """
    nstate = 0
    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'Excited state' in line:
                nstate += 1
    return nstate


def s2_to_s(s2):
    """
    get total spin (s) from s^2
    :param: s2
    :return: total spin (s)
    """
    return 0.5 * (-1 + np.sqrt(1 + 4 * s2))


def s_to_s2(s):
    """
    get s from s^2
    :param: s
    :return: s^2
    """
    return ((2*s + 1)**2 - 1) / 4


def take_states(ff, linne):
        if 'State A: Ground state' in linne:
            state_a = str('0')
        else:
            state_a = str(linne.split()[-1])
        linne = next(ff)
        state_b = str(linne.split()[-1])
        return state_a, state_b


def get_spins(filee, initial_selected_state, states_option, ground_statee, multiplicity_diction, selected_multiplicity):
    """
    Obtain three dictionaries:
    i) SOCs of the selected states
    ii) Approximate spins of the selected states (i.e. s2=0.7625 -> s2=0.75).
    iii) Each state real spin and its approximated value
    :param filee
    :param initial_selected_state
    :param multiplicity_diction
    :return: selected_states_spins, reals2_approxs2_dic
    """
    def change_state_name(state_spin, multip_dictionary, multip_count):
        state_multip = multip_dictionary[int(2 * s2_to_s(float(state_spin)) + 1)]
        multip_count[state_multip] += 1
        state = state_multip[0] + str(multip_count[state_multip])
        return state

    def get_selected_states(option, initial_states, ground_statee, all_states, multiplicity_selected, multiplicity_dict):
        """
        Depending on "states_option" it returns:
        0) "initial_states", ordered by multiplicity (S1, S2... T1, T2...)
        1) all states selected, ordered by multiplicity (S1, S2... T1, T2...)
        :param nstates:
        :param initial_states:
        :param multiplicity_selected:
        :param option:
        :return: final_states
        """
        final_states = []
        if option == 0:
            for i in initial_states:
                if i == 0:
                    final_states.append('0')
                else:
                    named_state = multiplicity_dict[multiplicity_selected][0] + str(i)
                    if named_state in all_states:
                        final_states.append(named_state)
                    else:
                        raise ValueError("The number of states selected selected must be among the total number of states "
                                     "selected calculated in QChem. Select a different number of states selected")

        elif option == 1:
            for i in all_states:
                if i == '0':
                    final_states.append('0')
                elif i[0] == multiplicity_dict[multiplicity_selected][0]:
                    final_states.append(i)

        if ground_statee not in final_states:
            raise ValueError("Ground state is not in selected states. Include it to continue the calculation.")
        else:
            final_states.remove(ground_statee)
            final_states.insert(0, ground_statee)
        return final_states

    all_states_spins = {}
    multiplicity_count = {'Singlet': 0, 'Doublet': 0, 'Triplet': 0, 'Quartet': 0, 'Quintet': 0, 'Sextet': 0, 'Heptet': 0}
    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'State A: Ground state' in line or 'State A: Root 0' in line:
                bra_state, ket_state = take_states(f, line)

                linne = search_word_line(f, 'Ket state')
                state_a_spin = linne.split()[-4]
                linne = search_word_line(f, 'Bra state')
                state_b_spin = linne.split()[-4]

                ket_state = change_state_name(state_b_spin, multiplicity_diction, multiplicity_count)
                all_states_spins.update({bra_state: state_a_spin})
                all_states_spins.update({ket_state: state_b_spin})

    all_states = list(all_states_spins.keys())
    selected_states = get_selected_states(states_option, initial_selected_state, ground_statee, all_states,
                                          selected_multiplicity, multiplicity_diction)

    selected_states_spins = {}
    for i in selected_states:
        if i in list(all_states_spins.keys()):
            selected_states_spins.update({i: all_states_spins[i]})

    if len(selected_states_spins) <= 1:
        raise ValueError("Multiplicity or initial states have not been well selected. ")
    return selected_states_spins, all_states, selected_states


def get_socs(filee, selected_state, ground_state, all_states_named):
    """
    Obtain three dictionaries:
    i) SOCs of the selected states
    ii) Approximate spins of the selected states (i.e. s2=0.7625 -> s2=0.75).
    iii) Each state real spin and its approximated value
    :param filee
    :param selected_state
    :param multiplicity_diction
    :return: selected_socs, selected_states_spins, reals2_approxs2_dictionary
    """
    def take_socs(ff, all_soc, interstate_string):
        """
        Gives i) Dictionary with socs between all states, ii) string "socs_interstates" to be used
        in the next loop in case SOC=0.
        :param ff
        :param all_soc
        :return: all_soc, socs_interstates
        """
        linee = search_word_line(ff, 'Mean-field SO (cm-1)')
        linee = search_word_line(ff, 'Actual matrix elements')
        linee = search_word_line(ff, '<Sz=')

        socs_interstates = []
        while ' <Sz=' in linee:
            linee = linee.replace(',', '|')
            linee = linee.replace('(', '|')
            linee = linee.replace('\n', '')
            linee = linee.replace(')', '')
            linee = linee.split('|')

            socs_sz = []
            for i in range(2, len(linee), 2):
                numero = complex(float(linee[i]), float(linee[i + 1]))
                socs_sz.append(str(numero))
            socs_interstates.append(socs_sz)
            linee = next(ff)
        all_soc.update({interstate_string: socs_interstates})
        return all_soc, socs_interstates

    all_socs = {}
    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'State A:' in line:
                if 'Root' in line:
                    bra_state, ket_state = take_states(f, line)
                    bra_state = all_states_named[int(bra_state)]
                    ket_state = all_states_named[int(ket_state)]
                else:
                    bra_state, ket_state = take_states(f, line)
                interstate = bra_state + "_" + ket_state

                line = search_word_line(f, 'Clebsh-Gordan coefficient')
                clebshgordan_coeff = float(line.split()[-1])
                if clebshgordan_coeff != 0:
                    all_socs, socs_interstate = take_socs(f, all_socs, interstate)
                elif clebshgordan_coeff == 0:
                    # In case no SOC can be calculated, SOCs = 0
                    all_socs.update({interstate: [['0j'] * len(socs_interstate)]})

    selected_socs = {}
    for i in selected_state:
        if i != ground_state:
            if all_states_named.index(i) < all_states_named.index(ground_state):
                inter_states = i + "_" + ground_state
            if all_states_named.index(i) > all_states_named.index(ground_state):
                inter_states = ground_state + "_" + i
            selected_socs.update({inter_states: all_socs[inter_states]})

    if selected_socs == {}:
        raise ValueError("SOC dictionary is empty. Selected states or selected multiplicity have not been well selected.")
    return selected_socs


def get_states_info(filee, selected_state, ground_statee, all_states):
    """
    Obtain:
    i) dictionary with the total and excitation energies of the selected states
    ii) list with all the states with its multiplicities ordered by energy: [T1, S1, T2..]
    iii) dictionary with the transitions of each configuration of each state
    :param filee:
    :param selected_state:
    :param multiplicity_diction:
    :param reals2_approxs2:
    :return: selected_states_energies
    """
    state_index = 1
    multip_count = {'Singlet': 0, 'Doublet': 0, 'Triplet': 0, 'Quartet': 0, 'Quintet': 0, 'Sextet': 0, 'Heptet': 0}

    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'Total energy in the final basis set' in line:
                ground_state_total_energy = float(line.split()[-1])
                all_states_energies = {'0': [ground_state_total_energy, 0.0000]}
            
            if 'beta electrons' in line:
                alpha = int(line.split()[2])
                beta = int(line.split()[5])
                orbital_dict = {}
                for i in range(1, beta+1):
                    orbital_dict.update({'D('+str(i)+')': i})
                for i in range(beta+1, alpha+1):
                    orbital_dict.update({'S('+str(i-beta)+')': i})
                for i in range(alpha+1, alpha+30):
                    orbital_dict.update({'V('+str(i-alpha)+')': i})
                config_list = []
                   
            if 'Excited state' in line:
                nstate = int(line.split()[2].replace(':',''))
                excit_energy = float(line.split()[-1])
                line = next(f)
                total_energy = float(line.split()[-2])

                line = next(f)
                if '<S**2>' in line:  # Word shown in doublets multiplicities in TDDFT
                    state_with_multip = all_states[nstate]
                    all_states_energies.update({state_with_multip: [total_energy, excit_energy]})
                elif 'Multiplicity' in line:  # Word shown in singlets multiplicities in TDDFT
                    multip_count[line.split()[-1]] += 1
                    state_with_multip = line.split()[-1][0] + str(multip_count[line.split()[-1]])
                    all_states_energies.update({state_with_multip: [total_energy, excit_energy]})

                # Go to configurations/transitions lines
                line = next(f)
                line = next(f)
                line = next(f)
                while len(line) > 2: # While line is not empty
                    # Each transition is saved in the dictionary
                    initial_orb = orbital_dict[''.join(line.split()[0:2])]
                    final_orb = orbital_dict[''.join(line.split()[3:5])]
                    config_list.append({'state': state_with_multip, 
                                       'transition/SOMO': str(initial_orb)+'-->'+str(final_orb), 
                                       'amplitude': float(line.split()[7])})
                    line = next(f)

    all_states_energy_order = list(all_states_energies.keys())
    selected_states_energies = {ground_statee: [all_states_energies[ground_statee][0], 0.000]}
    selected_config = []
    selected_config.append([{'state': 0,'transition/SOMO': list(range(beta+1, alpha+1)),'amplitude': 1.000}])

    for i in selected_state:
        if i != ground_statee:
            excit_energy = all_states_energies[i][1] - all_states_energies[ground_statee][1]
            selected_states_energies.update({i: [all_states_energies[i][0], excit_energy]})
            selected_config.append([d for d in config_list if d.get('state') == i])
    return selected_states_energies, all_states_energy_order, selected_config


def get_orbital_angmoment(filee, selected_state, ground_statee, all_states):
    """
    Obtain a dictionary with orbitals angular momentum between ground state and all the rest.
    :param filee:
    :param selected_state:
    :param all_states:
    :return: selected_momentums
    """

    all_momentums = {}
    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'Transition Angular Moments Between' in line:
                for i in range(0, 4):
                    line = next(f)
                while '----------------------------------' not in line:
                    bra_state = all_states[int(line.split()[0])]
                    ket_state = all_states[int(line.split()[1])]
                    interstate = bra_state + "_" + ket_state

                    x = line.split()[2] + "j"
                    y = line.split()[3] + "j"
                    z = line.split()[4] + "j"
                    all_momentums.update({interstate: [x, y, z]})
                    line = next(f)

    selected_momentums = {}
    for i in selected_state:
        if i != ground_statee:
            if all_states.index(i) < all_states.index(ground_state):
                inter_states = i + "_" + ground_state
            if all_states.index(i) > all_states.index(ground_state):
                inter_states = ground_state + "_" + i

            selected_momentums.update({inter_states: all_momentums[inter_states]})

    if selected_momentums == {}:
        raise ValueError("SOC dictionary is empty. Selected states or selected multiplicity have not been well selected.")
    return selected_momentums


def output_json(outpuut):
    outfile_name = outpuut + ".json"
    with open(outfile_name, 'w') as archivo:
        json.dump(output_dict, archivo, separators=(',', ':'), sort_keys=False, indent=4)


#################################
# BEGINNING OF PROGRAM     ######
#################################

# Take number of initial_states
totalstates = get_number_of_states(file)

# Take SOCs of the selected states
# This is done first since in interstate properties is where the approximate spin is defined
# WARNING! TDDFT states by default are going to have spin contamination.
approx_spin_dict, all_states, selected_states = get_spins(file, initial_states, state_selection, ground_state, multiplicity_dict,
                                         selected_multiplicity)

soclist_dict = get_socs(file, selected_states, ground_state, all_states)

# Take energies of the selected states
energylist_dict, all_states_energy_ordered, transitions_dict = get_states_info(file, selected_states, ground_state, all_states)

# Take orbital angular momentum of the selected states
orbitalmomentlist_dict = get_orbital_angmoment(file, selected_states, ground_state, all_states_energy_ordered)

output_dict = {
    "selected_energy_dict": energylist_dict,
    "soclist_dict": soclist_dict,
    "spin_dict": approx_spin_dict,
    "orbitalmomentlist_dict": orbitalmomentlist_dict,
    "transitions_dict": transitions_dict
}

output_json(file)
