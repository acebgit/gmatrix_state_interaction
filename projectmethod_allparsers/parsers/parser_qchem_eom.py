import json
import sys
import numpy as np
import pandas as pd

state_selection = 1  # 0: use selected states ; 1: use all states_selected
initial_states = [1]
symmetry_selection = 'B1'  # Symmetry selected states_selected

# file = str(sys.argv[1])
file = '../../molecules/triangulenes/2Tm_doublet_eomea.out'

#################################
# FUNCTIONS AND CLASSES    ######
#################################

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
                if clebshgordan_coeff != 0: # SOC exists and is taken 
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

                elif clebshgordan_coeff == 0: # SOC does not exist and is 0 
                    all_socs.update({interstate: [['0j'] * len(socs_interstate)]})
    return all_states_excitenergies, all_momentums, all_states_spins, all_socs


def get_selectedstates_and_totalenergies(filee, states_option, init_states, sym_selection):
    """
    It returns the states selection depending on "states_option" value (0: initial states, 
    1: all states in qchem output, 2: symmetry selection) and the total energies of the states. 
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
    states = select_states[1:]

    selected_socs = {}
    selected_momentums = {}
    for state_b in states:
        interstate = str(ground_state) + "_" + str(state_b)
        selected_socs.update({interstate: all_socs[interstate]})
        selected_momentums.update({interstate: all_momentums[interstate]})
    return selected_socs, selected_momentums


def output_json(outpuut):
    outfile_name = outpuut + ".json"
    with open(outfile_name, 'w') as archivo:
        json.dump(output_dict, archivo, separators=(',', ':'), sort_keys=False, indent=4)

#################################
# FUNCTIONS AND CLASSES 2  ######
#################################

def get_eom_type(eom_input):
    """
    Define the type of equation of motions eom_input
    :param: eom_input
    :return: eom_type
    """
    searches = ['Solving for']

    with open(eom_input) as data:
        for line in data:
            if any(i in line for i in searches):
                eom_type = line.split()[2].split('-')[0]
                break
    return eom_type


def get_orbital_symmetries(eom_input):

    lines_list = []
    with open(eom_input, encoding="utf8") as f:
        for line in f:
            line = search_word_line(f, ' Calculation will run')
            line = search_word_line(f, '-- Occupied --')
            while ' Beta MOs, Restricted' not in line:
                lines_list.append(line.split())
                line = next(f)
            break

    index_to_delete = [int(lines_list.index(i)) for i in lines_list if ('Occupied' in i) or ('Virtual' in i)]
    [lines_list.remove(lines_list[i]) for i in index_to_delete]
    symmetry_lines = [lines_list[i] for i in range(1, len(lines_list), 2)]

    orbitals = []
    orbitals += [(row[orbit] + row[orbit+1]) for row in symmetry_lines for orbit in range(0, len(row), 2)]
    return orbitals


def get_scf_energy(eom_input):
    """
    Define the type of equation of motions eom_input
    :param: eom_input
    :return: scf_energy
    """
    searches = ['SCF   energy in the final basis set']

    with open(eom_input, encoding="utf8") as data:
        for line in data:
            if any(i in line for i in searches):
                element = line[39:]
                scf_energy = float(element)
    return scf_energy


def get_energies(eom_input, eom_versioon):
    """
    Take the irreducible representations and the transitions in each irrep.,
    as well as their total and excitation energies
    :param: eom_input
    :return:  irreps, energies, excit_energies
    """
    searches = [eom_versioon + " transition"]
    all_states = []
    energies = {}

    with open(eom_input, encoding="utf8") as f:
        for line in f:
            if any(a in line for a in searches):
                state = line.split()[-1]
                all_states.append(state)
                line = next(f)
                energies.update({state: [float(line.split()[3]), float(line.split()[-2])]})
    return all_states, energies


def get_eom_socc_values(eom_input, optimized_state_index):
    """
    Take the contribution of hole configurations for each state
    :param: lines_eom_input
    :return: hole_contributions
    """
    search_1 = ' Mean-field SO (cm-1)'
    search_2 = 'SOCC '
    socc_list = []

    i = 0

    with open(eom_input, encoding="utf8") as file:
        for line in file:

            if search_1 in line:

                next_line = next(file)

                if search_2 in next_line:
                    if i == optimized_state_index:
                        socc_list.append('0.000000')
                    # SOCC in main state, in which interstate properties are calculated, is zero

                    elements = [next_line.split()]
                    socc_list.append(elements[0][2])

                    i += 1

    return socc_list


def get_maximum_amplitude_orbitals(eom_input, eom_versioon, cut_off):
    """
     Gives orbitals related with maximum amplitudes transitions
     :param: eom_input
     :return: significant_orbitals
     """
    def get_all_transitions(linee):
        def take_symmetry(lista):
            orbitals = []
            if lista == ['infty']:
                orbitals.append('infty')
            else:
                orbitals = [lista[i] + lista[i + 1] for i in range(0, len(lista), 3)]
            return orbitals

        all_transitions = []
        while '->' in linee:
            amplitude = float(linee.split()[0].replace('-', ''))
            orbitals_1 = take_symmetry(linee.split('->')[0].split()[1:])
            orbitals_2 = take_symmetry(linee.split('->')[1].split())
            all_transitions.append({'amplitude': amplitude, 'initial orbital': orbitals_1, 'final orbital': orbitals_2})
            linee = next(data)
        return all_transitions

    def mapping_symmetry_number(linee):
        """
         Used in "get_significant_orbitals". It takes the
         "Transitions between orbitals" lines without the amplitudes
         :param: linee
         :return: transition_orbital
         """
        linee = next(data)
        linee = next(data)
        linee = next(data)

        orbital_list = []
        while linee.isspace() ==  False:
            orbital_list.append({'symmetry': linee.split()[3]+linee.split()[4], 'number': linee.split()[0]})
            linee = next(data)
        return orbital_list

    def from_symmetry_to_number(transitions, mapping, nstate):
        final_list = []
        for i in transitions:
            init_orb = 0
            final_orb = 0
            for map in mapping:
                if i['initial orbital'][0] == map['symmetry']:
                    init_orb = map['number']
                elif i['final orbital'][0] == map['symmetry']:
                    final_orb = map['number']
            if init_orb == 0:
                init_orb = 'infty'
            if final_orb == 0:
                final_orb = 'infty'
            final_list.append(
                {'state': nstate, 
                 'transition/SOMO': init_orb+'->'+final_orb, 
                 'amplitude': i['amplitude']})
        return final_list

    searches = [eom_versioon + " transition"]
    norbitals_transition = []

    with open(eom_input) as data:
        for line in data:
            if any(i in line for i in searches):
                state = line.split()[-1]
                line = search_word_line(data, 'Transitions between orbitals')
                line = next(data)
                symmetry_transitions = get_all_transitions(line)

                selected_transitions = [i for i in symmetry_transitions if i['amplitude'] >= cut_off]
                if selected_transitions == []:
                    max_amplitude = max([i['amplitude'] for i in symmetry_transitions])
                    selected_transitions = [i for i in symmetry_transitions if i['amplitude'] >= max_amplitude]
                
                line = search_word_line(data, 'Summary of significant orbitals')
                mapping_symmetry_norbital = mapping_symmetry_number(line)

                transitions = from_symmetry_to_number(selected_transitions, mapping_symmetry_norbital, state)
                norbitals_transition.append(transitions)
    return norbitals_transition


def second_smallest_number(numbers):
    m1 = m2 = float('inf')
    for x in numbers:
        if x <= m1:
            m1, m2 = x, m1
        elif x < m2:
            m2 = x
    return m2


def eom_results(eom_presentation_list, irreps, total_energies, excitation_energies, eom_soc_constants, orbitals):
    eom_state = 0
    transition = 0
    selected_eom_excitation_energies_list = []
    selected_eom_soc_constants_list = []

    threshold_excitation_energy = 8  # in eV, energy difference with respect to ground state

    for i in range(0, len(irreps)):
        state_irreps = irreps[i][1]
        excit_energy = np.round(float(total_energies[i]), 3)
        excitation_energy = np.round(float(excitation_energies[i] * 27.211399), 3)
        socc = np.round(float(eom_soc_constants[i]), 2)
        # socc = 0

        if excitation_energy <= threshold_excitation_energy:
            selected_eom_excitation_energies_list.append(excitation_energies[i])
            # selected_eom_soc_constants_list.append(eom_soc_constants[i])
            selected_eom_soc_constants_list.append(socc)

            transition += 1

            eom_presentation_list.append([state_irreps, transition, excit_energy,
                                          excitation_energy, socc, orbitals[i * 3], orbitals[i * 3 + 1],
                                          orbitals[i * 3 + 2]])
            eom_state += 1

        if (i < len(irreps) - 1) and (irreps[i][1] != irreps[i + 1][1]):
            eom_presentation_list.append(['---'])

    second_smallest_excit_energy = second_smallest_number(excitation_energies)
    if second_smallest_excit_energy > threshold_excitation_energy:
        print('Increase the threshold_excitation_energy')
        exit()

    selected_excitation_energies_eom = np.array(selected_eom_excitation_energies_list, dtype=float)
    selected_eom_soc_constants = np.array(selected_eom_soc_constants_list, dtype=float)

    return eom_presentation_list, selected_excitation_energies_eom, selected_eom_soc_constants


def get_list_all_transitions(file, transition_cutoff, ref_state): 
    # orbital_sym = get_orbital_symmetries(file)

    eom_version = get_eom_type(file)

    scf_reference_energy = get_scf_energy(file)

    states, total_energies = get_energies(file, eom_version)

    # allstates_excit_energies, allstates_momentums, allstates_spins, allstates_socs = get_interstate_properties(file)

    # eom_socc = get_eom_socc_values(eom_input)

    orbitals_in_transitions = get_maximum_amplitude_orbitals(file, eom_version, transition_cutoff)

    presentation_list = []
    row_list = []
    for transition in list(total_energies.keys()):
        for trans in orbitals_in_transitions:
            if trans[0]['transition'] == transition:
                ener_excit = np.round(total_energies[transition][1] - total_energies[ref_state][1], 4)
                for subtrans in trans:
                    presentation_list.append([total_energies[transition][0], ener_excit, subtrans['orbitals']])
                    row_list.append(transition)
    df = pd.DataFrame(np.array(presentation_list), index=row_list, columns=['Total energ', 'Excit. ener.', 'transition'])

    print("")
    print("------------------------")
    print("    eom-CC ANALYSIS ")
    print("------------------------")
    print("SCF   reference energy: ", scf_reference_energy)
    print()
    print(df.to_string())


def get_transitions_json(filee): 
    eom_version = get_eom_type(filee)

    states, total_energies = get_energies(filee, eom_version)

    orbitals_in_transitions = get_maximum_amplitude_orbitals(filee, eom_version, cut_off=0)

    filtered_transitions = []
    [filtered_transitions.append(d) for state in selected_states for d in orbitals_in_transitions if d[0].get('state') == state]
    return filtered_transitions


#################################
# BEGINNING OF PROGRAM     ######
#################################

# 1) Take all the data
allstates_excit_energies, allstates_momentums, allstates_spins, allstates_socs = get_interstate_properties(file)

# 2) Take the selected initial_states and selected states energies
selected_states, allstates_energies = get_selectedstates_and_totalenergies(file, state_selection, initial_states,
                                                                           symmetry_selection)

# 3) Take properties of the selected states
energylist_json, approx_spin_json = get_selected_energies_spins(selected_states, allstates_energies,
                                                                allstates_excit_energies, allstates_spins)

soclist_json, orbitalmomentlist_json = get_selected_socs_orbitmoments(selected_states, allstates_socs,
                                                                      allstates_momentums)

# 4) Take all the transitions of the selected states
transitions_json = get_transitions_json(file)

output_dict = {
    "selected_energy_dict": energylist_json,
    "soclist_dict": soclist_json,
    "spin_dict": approx_spin_json,
    "orbitalmomentlist_dict": orbitalmomentlist_json,
    "transitions_dict": transitions_json
}

output_json(file)
