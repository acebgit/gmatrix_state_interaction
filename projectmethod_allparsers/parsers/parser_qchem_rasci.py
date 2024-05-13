import json
import numpy as np
import sys 
from pyqchem.parsers.parser_rasci import parser_rasci
from projectmethod.parsers.parser_gtensor import get_number_of_states, get_selected_states, \
    get_symmetry_states
from projectmethod.parsers.parser_excitstates import get_configurations_unpaired_orbitals

state_selection = 1  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
initial_states = [1, 2, 3, 4, 5]
symmetry_selection = 'B1u'  # Symmetry selected states_selected

# # file = str(sys.argv[1])
file = "../../projectmethod_allparsers/test/qchem_rasci.out"

#################################
# FUNCTIONS AND CLASSES    ######
#################################
with open(file, encoding="utf8") as f:
    output = f.read()
output = parser_rasci(output)


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


def get_states_somos(filee, selectedstates):
    """
    Get SOMOs of each of the configuration of each of the states. Output: list of dictionaries with
    each configuration of each state, its SOMOs and its amplitudes. 
    :param: file, selected_states
    """
    configuration_amplitudes, configuration_orbitals, initial_active_orbitals = \
        get_configurations_unpaired_orbitals(filee, selectedstates, cutoff=0)
    transitions = []
    for configuration in configuration_orbitals:
        state = int(configuration['State'])
        amplitude = configuration_amplitudes[configuration_orbitals.index(configuration)]['amplitude']
        transitions.append({'state': state, 
                                'transition/SOMO': configuration['SOMO orbitals'], 
                                'amplitude': np.round(amplitude, 4)})

    state_values = list(set([d['state'] for d in transitions if 'state' in d]))
    transitions_ordered = []
    for state in state_values:
        transitions_ordered.append([d for d in transitions if d['state'] == state])
    return transitions_ordered


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
selected_states = get_selected_states(file, totalstates, initial_states, state_selection, symmetry_selection)

# Take initial_states symmetry
all_symmetry, all_states_symmetry = get_symmetry_states(file, totalstates)

# Take the selected initial_states by symmetry
selected_states_symmetry = [all_states_symmetry[nstate - 1] for nstate in selected_states]

# 1) Take energies of the selected states
energylist_dict = get_energies(output, all_states_symmetry, selected_states_symmetry)

# 2) Take spin of the selected states
spin_dict = get_spin_momentum(file, all_states_symmetry, selected_states_symmetry)

# 3) Take SOCs of the selected states
ground_state = selected_states[0]-1
all_pairstates = [all_states_symmetry[ground_state] + "_" + state_sym for state_sym in all_states_symmetry]

selected_pairstates = [all_states_symmetry[ground_state] + "_" + state_sym for state_sym in selected_states_symmetry
                       if state_sym != all_states_symmetry[ground_state]]

soclist_dict = get_socs(output, all_pairstates, selected_pairstates)

# 4) Take orbital angular momentum of the selected states
orbitalmomentlist_dict = get_orbital_angmoment(output, all_pairstates, selected_pairstates)

# 5) Take SOMOs of each of the configurations of each state
transitions_json = get_states_somos(file, selected_states)

# 6) Put it in the output
output_dict = {
    "selected_energy_dict": energylist_dict,
    "soclist_dict": soclist_dict,
    "spin_dict": spin_dict,
    "orbitalmomentlist_dict": orbitalmomentlist_dict,
    "transitions_dict": transitions_json
}

output_json(file)
