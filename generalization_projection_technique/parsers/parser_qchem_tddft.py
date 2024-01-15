import json
# from procedure_projection_technique.parsers.parser_gtensor import *

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


def get_energies_spins(filee, selected_state):
    """
    Obtain a dictionary with energies and excitation energies of the selected states.
    :param filee:
    :return:
    """
    with open(filee, encoding="utf8") as f:
        for line in f:
            if '<S^2> =  ' in line:
                ground_state_spin = line.split()[-1]
            if ' Total energy in the final basis set' in line:
                ground_state_total_energy = line.split()[-1]

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


def get_energies(outpuut, all_states, selected_state):
    """
    Obtain a dictionary with energies and excitation energies of the selected states:
    :param outpuut:
    :param all_states:
    :param selected_state:
    :return: selected_states_energies
    """
    data = outpuut['excited_states']
    all_states_energies = {}
    for i in range(0, len(all_states)):
        lista = [data[i]['total_energy'], data[i]['excitation_energy']]
        all_states_energies.update({all_states[i]: lista})

    selected_states_energies = {}
    for i in selected_state:
        selected_states_energies.update({i: all_states_energies[i]})
    return selected_states_energies


def get_spin_momentum(outpuut, all_state, selected_state):
    """
    Get s2 of each state from Q-Chem output.
    :param outpuut:
    :param all_state:
    :param selected_state:
    :return:
    """
    search = ['  <S^2>      : ']
    s2_all_states = []
    with open(outpuut, encoding="utf8") as outpuut:
        for line in outpuut:
            if any(ii in line for ii in search):
                line = line.split()
                s2_all_states.append(line[2])

    all_spin = {}
    for i in range(0, len(all_state)):
        all_spin.update({all_state[i]: s2_all_states[i]})

    selected_spin = {}
    for i in selected_state:
        selected_spin.update({i: all_spin[i]})
    return selected_spin


def get_socs(outpuut, all_pairstate, selected_pairstate):
    """
    Obtain a dictionary with SOCs between ground state and all the rest, written as strings.
    :param outpuut:
    :param all_pairstate:
    :param selected_pairstate:
    :return:
    """
    data = outpuut['interstate_properties']
    all_socs = {}
    for i in range(2, len(all_pairstate) + 1):
        pair_socs = data[(1, i)]['total_soc_mat']

        soc_list = []
        for row in range(0, len(pair_socs)):
            row_soc = []
            for col in range(0, len(pair_socs[0])):
                row_soc.append(str(pair_socs[row][col]))
            soc_list.append(row_soc)
        all_socs.update({all_pairstate[i - 1]: soc_list})

    selected_socs = {}
    for i in selected_pairstate:
        selected_socs.update({i: all_socs[i]})
    return selected_socs


def get_orbital_angmoment(outpuut, all_pairstate, selected_pairstate):
    """
    Obtain a dictionary with orbitals angular momentum between ground state and all the rest, written as strings.
    :param outpuut:
    :param all_pairstate:
    :param selected_pairstate:
    :return:
    """
    data = outpuut['interstate_properties']
    all_momentums = {}
    for i in range(2, len(all_pairstate) + 1):
        pair_moment = data[(1, i)]['angular_momentum']
        pair_strings = [str(pair_moment[0]), str(pair_moment[1]), str(pair_moment[2])]
        all_momentums.update({all_pairstate[i - 1]: pair_strings})

    selected_momentums = {}
    for i in selected_pairstate:
        selected_momentums.update({i: all_momentums[i]})
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
energylist_dict, spin_dict = get_energies_spins(file, selected_states)

# 2) Take SOCs of the selected states
soclist_dict = get_socs(file)
exit()

# 4) Take orbital angular momentum of the selected states
orbitalmomentlist_dict = get_orbital_angmoment(output, all_pairstates, selected_pairstates)

output_dict = {
    "selected_energy_dict": energylist_dict,
    "soclist_dict": soclist_dict,
    "spin_dict": spin_dict,
    "orbitalmomentlist_dict": orbitalmomentlist_dict,
}

output_json(file)
