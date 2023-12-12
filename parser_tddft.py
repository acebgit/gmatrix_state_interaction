import numpy as np

archivo = '../tddft_outs/naphthalene_tddft.out'
cutoff_amp = 0.25
cutoff_soccs = 10


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

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                # 1) Take the excitation energies
                line = line.split()
                states = line[2].replace(':', '')
                excitation_energy = line[7]
                mapping_dict = {'state': states, 'energy': excitation_energy}
                energy_list.append(mapping_dict)

                # 2) Take the transition angular momentums
                for jump in range(0, 3):
                    next_line = next(file)
                element = next_line.split()
                mapping_dict = {'state': states, 'moment_x': element[2], 'moment_y': element[4], 'moment_z': element[6]}
                trans_moment_list.append(mapping_dict)

                # 3) Take the transitions
                for jump in range(0, 2):
                    next_line = next(file)
                excitation = 0

                while 'amplitude' in next_line:
                    next_line = next_line.split()
                    element = ' '.join(next_line[0:5]+next_line[-1:])
                    amplitude = next_line[-2]  # ' '.join(next_line[-2:])

                    mapping_dict = {'state': states, 'excitation': excitation,
                                    'transition': element, 'amplitude': amplitude}
                    excitation += 1
                    transitions_list.append(mapping_dict)
                    next_line = next(file)
    # [print(i) for i in transitions_list]
    return energy_list, trans_moment_list, transitions_list


def get_soccs(file, totalstate):
    """
    Get spin-orbit coupling constants between GS and all the excited states in TD-DFT.
    :parameter file
    :parameter totalstate
    :return energy_list, transitions_list
    """
    word_search = 'Mean-field SO (cm-1)'
    socc_list = []
    states = 1

    with open(file, encoding="utf8") as file:
        for line in file:
            if word_search in line:
                for jump in range(0, 11):  # Take the transitions
                    next_line = next(file)

                element = next_line.split()
                mapping_dict = {'state': states, 'socc': element[2]}
                socc_list.append(mapping_dict)

                states += 1
                if states > totalstate:
                    break
    # [print(i) for i in socc_list]
    return socc_list


somos = get_somos(archivo)
totalstates = get_number_of_states(archivo)
energy, transit_moments, transitions = get_energy_and_transitions(archivo)
soccs = get_soccs(archivo, totalstates)

# [print(i) for i in transitions]
presentation_list = [['State', 'Transition', 'Amp.', 'ΔE (eV)', 'SOCC', 'TMx', 'TMy', 'TMz']]
for i in range(0, len(transitions)):
    state = int(transitions[i]['state']) - 1
    trans = transitions[i]['transition']
    amp = np.round(float(transitions[i]['amplitude']), 2)

    energ = np.round(float(energy[state]["energy"]), 2)
    soc = np.round(float(soccs[state]["socc"]), 0)
    trans_x = transit_moments[state]["moment_x"]
    trans_y = transit_moments[state]["moment_y"]
    trans_z = transit_moments[state]["moment_z"]

    if abs(amp) >= cutoff_amp and soc >= cutoff_soccs:
        presentation_list.append([state, trans, amp, energ, soc, trans_x, trans_y, trans_z])
presentation_matrix = np.array(presentation_list, dtype=object)

print('cut-off amplitude', cutoff_amp, ', cut-off soccs: ', cutoff_soccs)
print('SOMOs: ', somos)
print('\n'.join(''.join('{:^15}'.format(item) for item in row)
                for row in presentation_matrix[:, :]))
