import json
import numpy as np

state_selection = 0  # 0: use "state_ras" ; 1: use all states_selected
initial_states = [1, 2, 3]
symmetry_selection = 'B1u'  # Symmetry selected states_selected
soc_options = 0

file = '../test/qchem_tddft.out'  # str(sys.argv[1])

#################################
# FUNCTIONS AND CLASSES    ######
#################################


def get_number_of_states(filee):
    """
    Obtain the total number of states selected in TD-DFT.
    :param: file
    :return: nstates
    """
    nstate = 0
    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'Excited state' in line:
                nstate += 1
    return nstate


def get_selected_states(nstates, selected_states, states_option):
    """
    Depending on "states_option" it returns states: 0) in "state_ras" 1) all states selected
    :param: file, nstates, selected_states, states_option, symmetry_selection
    :return: selected_states
    """
    if states_option == 0:
        for i in selected_states:
            if i <= 0 or i > nstates:
                raise ValueError("The number of states selected selected must be among the total number of states "
                                 "selected calculated in QChem. "
                                 "Select a different number of states selected")

    elif states_option == 1:
        selected_states = list(range(1, nstates + 1))
    return selected_states


def s2_to_s(s2):
    """
    get total spin (s) from s^2
    :param: s2_all
    :return: total spin (s)
    """
    return 0.5 * (-1 + np.sqrt(1 + 4 * s2))


def s_to_s2(s):
    """
    get s^2 from total spin (s)
    :param: s2_all
    :return: total spin (s)
    """
    return ( (2*s + 1)**2 - 1 ) / 4

def get_energies_realspins(filee, selected_state):
    """
    Obtain a dictionary with energies and excitation energies of the selected states.
    :param filee:
    :return:
    """
    # def approximate_spins(selected_states_spin):
    #     odd_multip = [i for i in range(1,20,2)]
    #     even_multip = [i for i in range(2,20,2)]
    #
    #     print(selected_states_spin)
    #
    #     selected_states_multip = [2*s2_to_s(float(selected_states_spin[i]))+1 for i in selected_states_spin]
    #     minimum_multiplicity = float(min(selected_states_multip))
    #     s2_values = [s_to_s2(i/4-1) for i in np.arange(minimum_multiplicity,20,2)]
    #     print(selected_states_multip)
    #     print(odd_multip, even_multip)
    #     exit()
    #
    #     for i in range(0, len(selected_states_spin)):
    #         diff_list = []
    #         for s2 in s2_values:
    #             diff = abs(float(selected_states_spin[i]) - s2)
    #             diff_list.append(diff)
    #
    #         min_index = diff_list.index(min(diff_list))
    #         selected_states_spin.update({i: str(s2_values[min_index])})
    #     return selected_states_spin

    with open(filee, encoding="utf8") as f:
        for line in f:
            if '<S^2> =  ' in line:
                ground_state_spin = line.split()[-1]
            if ' Total energy in the final basis set' in line:
                ground_state_total_energy = float(line.split()[-1])
                break

    state = 0
    all_states_energies = {}
    all_states_spins = {}
    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'Excited state' in line:
                state += 1

                # Take excitation energies (eV) and total energies (au)
                excitation_energy_ev = float(line.split()[-1])
                next_line = next(f)
                state_total_energy = float(next_line.split()[-2])
                all_states_energies.update({state: [state_total_energy, excitation_energy_ev]})

                # Take each state spin
                next_line = next(f)
                all_states_spins.update({state: next_line.split()[-1]})

    selected_states_energies = {0: [ground_state_total_energy, 0.0000]}
    selected_states_spins = {0: ground_state_spin}
    for i in selected_state:
        selected_states_energies.update({i: all_states_energies[i]})
        selected_states_spins.update({i: all_states_spins[i]})
    return selected_states_energies, selected_states_spins


def get_socs_approxspins(filee, selected_state):
    def search_word_line(archivo, linee, text):
        linee = next(file)
        while text not in linee:
            linee = next(archivo)
        return linee

    all_socs = {}
    all_states_spins = {}
    with open(filee, encoding="utf8") as file:
        for line in file:
            if 'State A: Root 0' in line:
                state_a = int(line.split()[-1])
                line = next(file)
                state_b = int(line.split()[-1])
                interstate = "0_" + str(state_b)

                line = search_word_line(file, line, 'Ket state')
                state_a_spin = line.split()[-4]
                line = search_word_line(file, line, 'Bra state')
                state_b_spin = line.split()[-4]
                all_states_spins.update({state_a: state_a_spin})
                all_states_spins.update({state_b: state_b_spin})

                line = search_word_line(file, line, 'Mean-field SO (cm-1)')
                line = search_word_line(file, line, 'Actual matrix elements')
                line = next(file)
                line = next(file)

                socs_interstate = []
                while ' <Sz=' in line:
                    line = line.replace(',', '|')
                    line = line.replace('(', '|')
                    line = line.replace('\n', '')
                    line = line.replace(')', '')
                    line = line.split('|')

                    socs_sz = []
                    for i in range(2, len(line), 2):
                        numero = complex(float(line[i]), float(line[i+1]))
                        socs_sz.append(str(numero))
                    socs_interstate.append(socs_sz)
                    line = next(file)
                all_socs.update({interstate: socs_interstate})

    selected_socs = {}
    selected_states_spins = {0: all_states_spins[0]}
    for state_b in selected_state:
        interstate = "0_" + str(state_b)
        selected_socs.update({interstate: all_socs[interstate]})
        selected_states_spins.update({state_b: all_states_spins[state_b]})
    return selected_socs, selected_states_spins


def get_orbital_angmoment(filee, selected_state):
    """
    Obtain a dictionary with orbitals angular momentum between ground state and all the rest, written as strings.
    """
    all_momentums = {}
    with open(filee, encoding="utf8") as file:
        for line in file:
            if 'Transition Angular Moments Between Ground and Excited States' in line:
                for i in range(0, 4): line = next(file)

                while '----------------------------------' not in line:
                    interstate = line.split()[0] + "_" + line.split()[1]
                    x = line.split()[2] + "j"
                    y = line.split()[3] + "j"
                    z = line.split()[4] + "j"
                    all_momentums.update({interstate: [x, y, z]})
                    line = next(file)

    selected_momentums = {}
    for state_b in selected_state:
        interstate = "0_" + str(state_b)
        selected_momentums.update({interstate: all_momentums[interstate]})
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
selected_states = get_selected_states(totalstates, initial_states, state_selection)

# 1) Take energies and spins of the selected states
# WARNING! TDDFT states by default are going to have spin contamination.
energylist_dict, real_spin_dict = get_energies_realspins(file, selected_states)

# 2) Take SOCs of the selected states
soclist_dict, approx_spin_dict = get_socs_approxspins(file, selected_states)

# 4) Take orbital angular momentum of the selected states
orbitalmomentlist_dict = get_orbital_angmoment(file, selected_states)

output_dict = {
    "selected_energy_dict": energylist_dict,
    "soclist_dict": soclist_dict,
    "spin_dict": approx_spin_dict,
    "orbitalmomentlist_dict": orbitalmomentlist_dict,
}

output_json(file)
