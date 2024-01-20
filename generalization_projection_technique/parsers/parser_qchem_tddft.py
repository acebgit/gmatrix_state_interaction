import json
import numpy as np

state_selection = 0  # 0: use "state_ras" ; 1: use all states_selected
initial_states = [0, 1, 2, 3, 4, 5] # 0 is the ground state
ground_state = 'D1'
selected_multiplicity = 2

file = '../../generalization_projection_technique/test/qchem_tddft_doublets.out'  # str(sys.argv[1])

#################################
# FUNCTIONS AND CLASSES    ######
#################################
multiplicity_dict = {1: 'Singlet', 2: 'Doublet', 3: 'Triplet', 4: 'Quartet', 5: 'Quintet', 6: 'Sextet'}


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


def get_selected_states(nstates, states_selected, multiplicity_selected, states_option, ground_statee):
    """
    Depending on "states_option" it returns:
    0) "states_selected", ordered by multiplicity (S1, S2... T1, T2...)
    1) all states selected, ordered by multiplicity (S1, S2... T1, T2...)
    :param nstates:
    :param states_selected:
    :param multiplicity_selected:
    :param states_option:
    :return: states_with_multip
    """
    states_with_multip = []
    if states_option == 0:
        for i in states_selected:
            if i == 0:
                states_with_multip.append('0')
            elif i <= nstates:
                states_with_multip.append(multiplicity_selected[0] + str(i))
            else:
                raise ValueError("The number of states selected selected must be among the total number of states "
                                 "selected calculated in QChem. "
                                 "Select a different number of states selected")

    elif states_option == 1:
        states_with_multip = [(multiplicity_selected[0] + str(i)) for i in range(1, nstates + 1)]

    if ground_statee not in states_with_multip:
        raise ValueError("Ground state is not in selected states. Include it to continue the calculation.")
    else:
        states_with_multip.remove(ground_statee)
        states_with_multip.insert(0, ground_statee)
    return states_with_multip


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


def get_spins(filee, selected_state, multiplicity_diction):
    """
    Obtain three dictionaries:
    i) SOCs of the selected states
    ii) Approximate spins of the selected states (i.e. s2=0.7625 -> s2=0.75).
    iii) Each state real spin and its approximated value
    :param filee
    :param selected_state
    :param multiplicity_diction
    :return: selected_states_spins, reals2_approxs2_dic
    """
    def take_states(ff, linne):
        if 'State A: Ground state' in linne:
            state_a = str('0')
        else:
            state_a = str(linne.split()[-1])
        linne = next(ff)
        state_b = str(linne.split()[-1])
        return state_a, state_b

    def search_word_line(archivo, text):
        """
        Search a text and goes to that line.
        :param archivo
        :param text
        :return: linee
        """
        linee = next(f)
        while text not in linee:
            linee = next(archivo)
        return linee

    def change_state_name(state_spin, multip_dictionary, multip_count):
        state_multip = multip_dictionary[int(2 * s2_to_s(float(state_spin)) + 1)]
        multip_count[state_multip] += 1
        state = state_multip[0] + str(multip_count[state_multip])
        return state

    all_states_spins = {}
    reals2_approxs2_dic = {}
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

                state_b_real_spin = linne.split()[5]
                reals2_approxs2_dic.update({np.round(float(state_b_real_spin), 4): np.round(float(state_b_spin), 4)})

    all_states_named = list(all_states_spins.keys())

    selected_states_spins = {}
    for i in selected_state:
        if i in list(all_states_spins.keys()):
            selected_states_spins.update({i: all_states_spins[i]})

    if len(selected_states_spins) <= 1:
        raise ValueError("Multiplicity or initial states have not been well selected. ")
    return selected_states_spins, reals2_approxs2_dic, all_states_named


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
    def take_states(ff, linne):
        if 'State A: Ground state' in linne:
            state_a = str('0')
        else:
            state_a = str(linne.split()[-1])
        linne = next(f)
        state_b = str(linne.split()[-1])
        return state_a, state_b

    def search_word_line(archivo, text):
        """
        Search a text and goes to that line.
        :param archivo
        :param text
        :return: linee
        """
        linee = next(f)
        while text not in linee:
            linee = next(archivo)
        return linee

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
    all_states_spins = {}
    reals2_approxs2_dic = {}
    multip_count = {'Singlet': 0, 'Doublet': 0, 'Triplet': 0, 'Quartet': 0, 'Quintet': 0, 'Sextet': 0, 'Heptet': 0}

    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'State A:' in line:
                bra_state, ket_state = take_states(f, line)
                bra_state = all_states_named[int(bra_state)]
                ket_state = all_states_named[int(ket_state)]

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


def get_energies(filee, selected_state, multiplicity_diction, reals2_approxs2):
    """
    Obtain:
    i) dictionary with the total and excitation energies of the selected states
    ii) List with all the states with its multiplicities ordered by energy: [T1, S1, T2..]
    :param filee:
    :param selected_state:
    :param multiplicity_diction:
    :param reals2_approxs2:
    :return: selected_states_energies
    """
    with open(filee, encoding="utf8") as f:
        for line in f:
            if ' Total energy in the final basis set' in line:
                ground_state_total_energy = float(line.split()[-1])
                break

    all_states_energies = {}
    all_states = []
    multip_count = {'Singlet': 0, 'Doublet': 0, 'Triplet': 0, 'Quartet': 0, 'Quintet': 0, 'Sextet': 0, 'Heptet': 0}
    with open(filee, encoding="utf8") as f:
        for line in f:

            if 'Excited state' in line:
                excitation_energy_ev = float(line.split()[-1])
                line = next(f)
                state_total_energy = float(line.split()[-2])

                line = next(f)

                if '<S**2>' in line:  # Word shown in doublets multiplicities in TDDFT
                    real_spin = float(line.split()[-1])
                    approx_spin = reals2_approxs2[real_spin]
                    multip = multiplicity_diction[int(2*s2_to_s(approx_spin)+1)]
                    multip_count[multip] += 1
                    state = multip[0] + str(multip_count[multip])

                elif 'Multiplicity' in line:  # Word shown in singlets multiplicities in TDDFT
                    multip_count[line.split()[-1]] += 1
                    state = line.split()[-1][0] + str(multip_count[line.split()[-1]])
                all_states_energies.update({state: [state_total_energy, excitation_energy_ev]})
                all_states.append(state)

    selected_states_energies = {'GS': [ground_state_total_energy, 0.0000]}
    for ket_state in selected_state:
        if ket_state in list(all_states_energies.keys()):
            selected_states_energies.update({ket_state: [all_states_energies[ket_state][0],
                                                         all_states_energies[ket_state][1]]})
    return selected_states_energies, all_states


def get_orbital_angmoment(filee, selected_state, all_states):
    """
    Obtain a dictionary with orbitals angular momentum between ground state and all the rest.
    :param filee:
    :param selected_state:
    :param all_states:
    :return: selected_momentums
    """
    all_momentums = {}
    multip_count = {
        'Singlet': 0,
        'Doublet': 0,
        'Triplet': 0,
        'Quartet': 0,
        'Quintet': 0,
        'Sextet': 0,
        'Heptet': 0,
    }
    with open(filee, encoding="utf8") as f:
        for line in f:

            # One multiplicity in the excited states
            if 'Transition Angular Moments Between Ground and Excited States' in line:
                for i in range(0, 4):
                    line = next(f)
                count = 0
                while '----------------------------------' not in line:
                    interstate = "GS_" + all_states[count]
                    x = line.split()[2] + "j"
                    y = line.split()[3] + "j"
                    z = line.split()[4] + "j"
                    all_momentums.update({interstate: [x, y, z]})
                    count += 1
                    line = next(f)

            # Two or more multiplicities in the excited states
            elif 'Transition Angular Moments Between Ground' in line:
                multip = line.split()[6]
                for i in range(0, 4):
                    line = next(f)
                while '----------------------------------' not in line:
                    multip_count[multip] += 1
                    interstate = "GS_" + multip[0] + str(multip_count[multip])
                    x = line.split()[2] + "j"
                    y = line.split()[3] + "j"
                    z = line.split()[4] + "j"
                    all_momentums.update({interstate: [x, y, z]})
                    line = next(f)

    selected_momentums = {}
    for ket_state in selected_state:
        interstate = "GS_" + str(ket_state)
        if interstate in list(all_momentums.keys()):
            selected_momentums.update({interstate: [all_momentums[interstate][0], all_momentums[interstate][1],
                                                    all_momentums[interstate][2]]})
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

# Take the selected initial_states
selected_states = get_selected_states(totalstates, initial_states, multiplicity_dict[selected_multiplicity],
                                      state_selection, ground_state)

# Take SOCs of the selected states
# This is done first since in interstate properties is where the approximate spin is defined
# WARNING! TDDFT states by default are going to have spin contamination.

approx_spin_dict, reals2_approxs2_dict, all_states_named = get_spins(file, selected_states, multiplicity_dict)

soclist_dict = get_socs(file, selected_states, ground_state, all_states_named)

# Take energies of the selected states
energylist_dict, all_state_list = get_energies(file, selected_states, multiplicity_dict, reals2_approxs2_dict)

# Take orbital angular momentum of the selected states
orbitalmomentlist_dict = get_orbital_angmoment(file, selected_states, all_state_list)

output_dict = {
    "selected_energy_dict": energylist_dict,
    "soclist_dict": soclist_dict,
    "spin_dict": approx_spin_dict,
    "orbitalmomentlist_dict": orbitalmomentlist_dict,
}

output_json(file)
