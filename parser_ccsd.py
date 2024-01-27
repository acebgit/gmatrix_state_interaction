
archivo = '../\
molecules/ccsd/anthracene_opt_ccsd.out'
cutoff_amp = 0.5
cutoff_soccs = 20
selected_multip = 2


def get_somos(file):
    """
    Obtain the SOMOs in TD-DFT.
    :param: file
    :return: totalstates
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    element = 0
    word_search = ['beta electrons']
    for line in data:
        if any(ii in line for ii in word_search):
            line = line.split()
            alpha = int(line[2])
            beta = int(line[5])
            somo_orbitals = list(range(beta+1, alpha+1))
            break
    return somo_orbitals


def get_number_of_states(file):
    """
    Obtain the total number of states selected in TD-DFT.
    :param: file
    :return: totalstates
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    element = 0
    word_search = ['NRoots =']
    for line in data:
        if any(ii in line for ii in word_search):
            line = line.split()
            element = line[2]
            element = element.replace(',', '')
            break
    totalstate = int(element)
    return totalstate


def get_energy_and_transitions(file):
    """
    Get excitation energies and orbital transitions in TD-DFT.
    :parameter file
    :return energy_list, transitions_list
    """
    word_search = 'Excited state'
    energy_list = []
    trans_moment_list = []
    transitions_list = []

    multip_dict = {'Singlet': 0, 'Doublet': 1,
                   'Triplet': 2, 'Quartet': 3,
                   'Quintet': 4, 'Sextet': 5}

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                # 1) Take the excitation energies
                line = line.split()
                states = int(line[2].replace(':', ''))
                excitation_energy = line[7]

                for jump in range(0, 2):
                    next_line = next(file)
                state_multip = next_line.split()[1]

                mapping_dict = {'state': states, 'multiplicity': multip_dict[state_multip], 'energy': excitation_energy}
                energy_list.append(mapping_dict)

                # 2) Take the transition angular momentums
                next_line = next(file)
                transmom = next_line.split()
                mapping_dict = {'state': states, 'multiplicity': multip_dict[state_multip], 'moment_x': transmom[2], 'moment_y': transmom[4], 'moment_z': transmom[6]}
                trans_moment_list.append(mapping_dict)

                # 3) Take the transitions
                for jump in range(0, 2):
                    next_line = next(file)
                excitation = 0

                while 'amplitude' in next_line:
                    next_line = next_line.split()
                    transmom = ' '.join(next_line[0:5]+next_line[-1:])
                    if 'alpha' in next_line or 'beta' in next_line:
                        amplitude = next_line[-2]  # ' '.join(next_line[-2:])
                    else:
                        amplitude = next_line[-1]  # ' '.join(next_line[-2:])

                    mapping_dict = {'state': states, 'multiplicity': multip_dict[state_multip], 'excitation': excitation,
                                    'transition': transmom, 'amplitude': amplitude}
                    excitation += 1
                    transitions_list.append(mapping_dict)
                    next_line = next(file)
    # [print(i) for i in transitions_list]
    return energy_list, trans_moment_list, transitions_list


def get_soccs(file, energy_dict):
    """
    Get spin-orbit coupling constants between GS and all the excited states in TD-DFT.
    :parameter file
    :return energy_list, transitions_list
    """
    word_search = 'State A: Ground state'
    socc_list = []
    multip_count = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:

                for jump in range(0, 4):  # Take the state multiplicity
                    next_line = next(file)
                state_multip = int(float(next_line.split()[-4]))

                for i in range(multip_count[state_multip], len(energy_dict)):  # Take the state number from energy list
                    if state_multip == energy_dict[i]['multiplicity']:
                        state_in_dict = energy_dict[i]['state']
                        multip_count[state_multip] = state_in_dict
                        break

                for jump in range(0, 30):  # Go to Mean-field SOCC
                    next_line = next(file)
                    if 'Clebsh-Gordon coefficient is too small' in next_line:
                        mapping_dict = {'state': state_in_dict, 'multiplicity': state_multip, 'socc': 0}
                        break

                element = next_line.split()
                if 'SOCC' in element:  # Take the Mean-field SOCC
                    mapping_dict = {'state': state_in_dict, 'multiplicity': state_multip, 'socc': element[2]}
                socc_list.append(mapping_dict)

    sorted_socc = sorted(socc_list, key=lambda x: x['state'])
    # [print(i) for i in sorted_socc]
    # exit()
    return sorted_socc


somos = get_somos(archivo)
totalstates = get_number_of_states(archivo)
energy, transit_moments, transitions = get_energy_and_transitions(archivo)
soccs = get_soccs(archivo, energy)

