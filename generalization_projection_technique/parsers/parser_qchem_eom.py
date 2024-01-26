import json

state_selection = 1  # 0: use selected states ; 1: use all states_selected
initial_states = [1, 2, 3, 4, 5]
symmetry_selection = 'B1'  # Symmetry selected states_selected

file = '../../\
molecules/eomccsd_outs/benzophenone_singlet_eomsf_sto3g.out'

#################################
# FUNCTIONS AND CLASSES    ######
#################################


def get_interstate_properties(filee):
    """
    Take interstate properties between ground state and all states:
    i) Excitation energy
    ii) Orbital angular momentums
    iii) Spins of the states
    iv) Spin-orbit couplings
    :param filee:
    :return: all_states_excitenergies, all_momentums, all_states_spins, all_socs
    """
    def search_word_line(archivo, text):
        linee = next(f)
        while text not in linee:
            linee = next(archivo)
        return linee

    all_socs = {}
    all_states_spins = {}
    all_states_excitenergies = {}
    all_momentums = {}

    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'State A' in line:
                # 1) Take state symmetries
                state_a = line.split()[-1]
                line = next(f)
                state_b = line.split()[-1]
                interstate = state_a + "_" + state_b

                # 2) Take state excitation energies
                line = search_word_line(f, 'Energy GAP')
                all_states_excitenergies.update({state_a: '0.00000'})
                all_states_excitenergies.update({state_b: line.split()[-2]})

                # 3) Take state orbital angular momentums
                line = search_word_line(f, 'Transition angular momentum')
                line = next(f)
                line = line.replace(',', '')
                line = line.replace(')', '')
                line = line.replace('i', 'j')
                all_momentums.update({interstate: [line.split()[2], line.split()[4], line.split()[6]]})

                # 4) Take state's spin
                line = search_word_line(f, 'Ket state')
                state_a_spin = line.split()[-4]
                line = search_word_line(f, 'Bra state')
                state_b_spin = line.split()[-4]
                all_states_spins.update({state_a: state_a_spin})
                all_states_spins.update({state_b: state_b_spin})

                # 5) Take SOCs
                line = search_word_line(f, 'Clebsh-Gordan coefficient')
                clebshgordan_coeff = float(line.split()[-1])
                if clebshgordan_coeff != 0:
                    line = search_word_line(f, 'Mean-field SO (cm-1)')
                    line = search_word_line(f, 'Actual matrix elements')
                    line = search_word_line(f, '<Sz=')

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
                        line = next(f)
                    all_socs.update({interstate: socs_interstate})

                elif clebshgordan_coeff == 0:
                    all_socs.update({interstate: [['0j'] * len(socs_interstate)]})
    return all_states_excitenergies, all_momentums, all_states_spins, all_socs


def get_selectedstates_and_totalenergies(filee, states_option, init_states, sym_selection):
    """
    Depending on "states_option" it returns states:
    0) initial states
    1) All states in qchem output
    :param filee:
    :param states_option:
    :param init_states:
    :param sym_selection:
    :return: select_states, total_energies
    """
    search = ['EOMIP transition ', 'EOMSF transition ', 'EOMEE transition ', 'EOMEA transition ']
    nstates = 0
    all_states = []
    total_energies = {}

    with open(filee, encoding="utf8") as f:
        for line in f:
            if any(a in line for a in search):
                state = line.split()[-1]
                all_states.append(state)
                line = next(f)
                total_energies.update({state: line.split()[3]})
                nstates += 1

            if 'State A' in line:
                ground_state = line.split()[-1]
                break

    select_states = []
    if states_option == 0:
        for i in init_states:
            if 0 < i <= nstates:
                check_state = str(i) + "/"
                for state in all_states:
                    if check_state in state:
                        select_states.append(state)
            else:
                raise ValueError("The number of states selected selected must be among the total number of states "
                                 "selected calculated in QChem. "
                                 "Select a different number of states selected")

    elif states_option == 1:
        select_states = all_states

    elif states_option == 2:
        for state in all_states:
            if sym_selection in state:
                select_states.append(state)
        if select_states:
            raise ValueError("There is not this symmetry. Change the symmetry selection.")

    if ground_state in select_states:
        select_states.remove(ground_state)
    select_states.insert(0, ground_state)
    return select_states, total_energies


def get_selected_energies_spins(select_states, all_totalenergies, all_excitenergies, all_spins):
    """
    Get selected states energies and spins.
    :param select_states:
    :param all_totalenergies:
    :param all_excitenergies:
    :param all_spins:
    :return: selected_states_energies, selected_states_spins
    """
    selected_states_energies = {}
    selected_states_spins = {}
    for i in select_states:
        selected_states_energies.update({i: [float(all_totalenergies[i]), float(all_excitenergies[i])]})
        selected_states_spins.update({i: all_spins[i]})
    return selected_states_energies, selected_states_spins


def get_selected_socs_orbitmoments(select_states, all_socs, all_momentums):
    """
    Get selected states orbital angular momentums.
    :param select_states:
    :param all_socs:
    :param all_momentums:
    :return: selected_socs, selected_momentums
    """
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
selected_states, allstates_energies = get_selectedstates_and_totalenergies(file, state_selection, initial_states,
                                                                           symmetry_selection)

# Take properties of the selected states
energylist_dict, approx_spin_dict = get_selected_energies_spins(selected_states, allstates_energies,
                                                                allstates_excit_energies, allstates_spins)

soclist_dict, orbitalmomentlist_dict = get_selected_socs_orbitmoments(selected_states, allstates_socs,
                                                                      allstates_momentums)

output_dict = {
    "selected_energy_dict": energylist_dict,
    "soclist_dict": soclist_dict,
    "spin_dict": approx_spin_dict,
    "orbitalmomentlist_dict": orbitalmomentlist_dict,
}

output_json(file)
