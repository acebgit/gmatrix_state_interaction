import json
import sys

state_selection = 1  # 0: use "state_ras" ; 1: use all states_selected
initial_states = [1,2]
symmetry_selection = 'B1'  # Symmetry selected states_selected

file = '../test/qchem_eomcc.out'

#################################
# FUNCTIONS AND CLASSES    ######
#################################

def get_interstate_properties(filee):
    def search_word_line(archivo, linee, text):
        linee = next(file)
        while text not in linee:
            linee = next(archivo)
        return linee

    all_socs = {}
    all_states_spins = {}
    all_states_excitenergies = {}
    all_momentums = {}

    with open(filee, encoding="utf8") as file:
        for line in file:
            if 'State A' in line:
                # 1) Take state symmetries
                state_a = line.split()[-1]
                line = next(file)
                state_b = line.split()[-1]
                interstate = state_a + "_" + state_b

                # 2) Take state excitation energies
                line = search_word_line(file, line, 'Energy GAP')
                all_states_excitenergies.update({state_a: '0.00000'})
                all_states_excitenergies.update({state_b: line.split()[-2]})

                # 3) Take state orbital angular momentums
                line = search_word_line(file, line, 'Transition angular momentum')
                line = next(file)
                line = line.replace(',', '')
                line = line.replace(')', '')
                line = line.replace('i', 'j')
                all_momentums.update({interstate: [line.split()[2], line.split()[4], line.split()[6]]})

                # 4) Take state's spin
                line = search_word_line(file, line, 'Ket state')
                state_a_spin = line.split()[-4]
                line = search_word_line(file, line, 'Bra state')
                state_b_spin = line.split()[-4]
                all_states_spins.update({state_a: state_a_spin})
                all_states_spins.update({state_b: state_b_spin})

                # 5) Take SOCs
                line = search_word_line(file, line, 'Clebsh-Gordan coefficient')
                clebshgordan_coeff = float(line.split()[-1])
                if clebshgordan_coeff != 0:
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

                elif clebshgordan_coeff == 0:
                    all_socs.update({interstate: [['0j'] * len(socs_interstate)]})
    return all_states_excitenergies, all_momentums, all_states_spins, all_socs


def get_selectedstates_and_totalenergies(filee, states_option, init_states, symmetry_selection):
    """
    Depending on "states_option" it returns states: 0) in "state_ras" 1) all states selected
    :param: file, nstates, selected_states, states_option, symmetry_selection
    :return: selected_states
    """
    nstates = 0
    all_states = []
    total_energies = {}

    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'EOMIP transition' in line:
                state = line.split()[-1]
                all_states.append(state)
                line = next(f)
                total_energies.update({state: line.split()[3]})
                nstates+=1

            if 'State A' in line:
                ground_state = line.split()[-1]
                break

    selected_states = []

    if states_option == 0:
        for i in init_states:
            if i > 0 and i <= nstates:
                check_state = str(i) + "/"
                for state in all_states:
                    if check_state in state:
                        selected_states.append(state)
            else:
                raise ValueError("The number of states selected selected must be among the total number of states "
                                 "selected calculated in QChem. "
                                 "Select a different number of states selected")

    elif states_option == 1:
        selected_states = all_states

    elif states_option == 2:
        for state in all_states:
            if symmetry_selection in state:
                selected_states.append(state)
        if selected_states == []:
            raise ValueError("There is not this symmetry. Change the symmetry selection.")

    if ground_state in selected_states:
        selected_states.remove(ground_state)
    selected_states.insert(0, ground_state)
    return selected_states, total_energies


def get_selected_energies_spins(select_states, all_totalenergies, all_excitenergies, all_spins):
    selected_states_energies = {}
    selected_states_spins = {}
    for i in select_states:
        selected_states_energies.update({i: [float(all_totalenergies[i]), float(all_excitenergies[i])]})
        selected_states_spins.update({i: all_spins[i]})
    return selected_states_energies, selected_states_spins


def get_selected_socs_orbitmoments(select_states, all_socs, all_momentums):
    ground_state = select_states[0]
    select_states.remove(ground_state)

    selected_socs = {}
    selected_momentums = {}
    for state_b in select_states:
        interstate = str(ground_state) + "_" + str(state_b)
        selected_socs.update({interstate: all_socs[interstate]})
        selected_momentums.update({interstate: all_momentums[interstate]})
    return selected_socs, selected_momentums


def output_json(outpuut):
    outfile_name = outpuut + ".json"
    with open(outfile_name, 'w') as archivo:
        json.dump(output_dict, archivo, separators=(',', ':'), sort_keys=False, indent=4)

#################################
# BEGINNING OF PROGRAM     ######
#################################

# Take all the data
allstates_excit_energies, allstates_momentums, allstates_spins, allstates_socs = get_interstate_properties(file)

# Take the selected initial_states
selected_states, allstates_energies = get_selectedstates_and_totalenergies(file, state_selection, initial_states, symmetry_selection)

# Take properties of the selected states
energylist_dict, approx_spin_dict = get_selected_energies_spins(selected_states, allstates_energies, allstates_excit_energies, allstates_spins)

soclist_dict, orbitalmomentlist_dict = get_selected_socs_orbitmoments(selected_states, allstates_socs, allstates_momentums)

output_dict = {
    "selected_energy_dict": energylist_dict,
    "soclist_dict": soclist_dict,
    "spin_dict": approx_spin_dict,
    "orbitalmomentlist_dict": orbitalmomentlist_dict,
}

output_json(file)
