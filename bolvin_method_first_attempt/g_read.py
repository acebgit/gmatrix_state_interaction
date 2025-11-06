import sys
import numpy as np


def get_number_of_states(filename):
    """
    Data to be obtained from the output:
    - Number of states
    - Eigenenergies
    - spin-orbit couplings
    - spin angular momentum
    - orbital angular momentum
    """

    with open(filename) as f:
        data = f.readlines()

    searches = ['Requested states: ']
    elements = []
    
    for line in data:

        if any(i in line for i in searches):

            element = line[29:33]
            elements.append(element.split())

    total_state = np.array(elements, dtype=float)
    total_states = int(total_state)

    return total_states


def get_states_same_symmetry(filename, ground_state, totalstates, symmetry_selection):
    """
    Search the symmetry of each state and choose those with the selected symmetry
    """
    with open(filename) as f:
        data = f.readlines()

    searches = ['Electronic state symmetry: ']
    elements = []

    for line in data:
        if any(i in line for i in searches):
            element = line[27:30]
            elements.append(element.split())

    state_symmetries = np.array(elements)

    symmetry_selected_states = []
    symmetry_selected_states.append(ground_state)

    for i in range(1, totalstates+1):
        if (state_symmetries[i-1] == symmetry_selection) and (i != ground_state):
            symmetry_selected_states.append(i)

    if symmetry_selected_states == [ground_state]:
        print('There is not this kind of symmetry: change the symmetry selection')
        sys.exit()

    return symmetry_selected_states


def get_eigenenergies(states_ras, filename):

    searches = [' RAS-CI total energy for state  ']  # words to be searched in sentences
    elements = []

    with open(filename) as f:
        data = f.readlines()

    for line in data:

        if any(i in line for i in searches):

            element = line[37:]
            elements.append(element.split())

    energies_selected = []

    for i in states_ras:
        energies_selected.append(elements[i - 1])
        energies_selected.append(elements[i - 1])

    eigen_energies = np.array( energies_selected, dtype=float )

    return eigen_energies


def get_spin_orbit_couplings(totalstates, states_ras, filename):

    searches = 'Total mean-field SOC matrix'
    elements = []

    with open(filename) as file:
        for line in file:

            if searches in line:
                next_line = next(file)
    
                if next_line.startswith("                     |Sz'=-0.50>"):
                    next_line = next(file)
    
                    while next_line.startswith(" <Sz="):
                        if next_line.startswith((" <Sz=-0.50|", " <Sz= 0.50|")):
                            elements.append(next_line[12:35].split())  # add elements |Sz'=-0.50> to a list
                            elements.append(next_line[37:61].split())  # add elements |Sz'=+0.50> to a list
                            next_line = next(file)
                        else:
                            next_line = next(file)
    
                if next_line.startswith("                     |Sz'=-1.50>"):
                    next_line = next(file)
    
                    while next_line.startswith(" <Sz="):
                        if next_line.startswith((" <Sz=-0.50|", " <Sz= 0.50|")):
                            elements.append(next_line[37:61].split()) # add elements |Sz'=-0.50> to a list
                            elements.append(next_line[63:87].split()) # add elements |Sz'=+0.50> to a list
                            next_line = next(file)
                        else:
                            next_line = next(file)

    # how states are selected: < B -1/2 | A -1/2 >, < B -1/2 | A +1/2 >, < B +1/2 | A -1/2 >, < B +1/2 | A +1/2 >
    soc_selected_list = []
    for state_B in states_ras:  # < B |
        for state_A in states_ras:  # | A >

            n_soc_per_state = 0
            while n_soc_per_state < 4:

                if (state_B == state_A):
                    soc_selected_list.append(['0.000000', '0.000000'])

                if (state_B > state_A):
                    val = 4 * (totalstates - 1) * (state_B - 1) + 4 * (state_A - 1) + n_soc_per_state
                    soc_selected_list.append( elements [ val ] )

                if (state_B < state_A):
                    val = 4 * (totalstates - 1) * (state_B - 1) + 4 * (state_A - 2) + n_soc_per_state
                    soc_selected_list.append( elements [ val ] )

                n_soc_per_state = n_soc_per_state + 1

    soc_selected = np.array(soc_selected_list, dtype=float)

    # how states are saved: | A > in rows, | B > in columns. Order of spin: +1/2 , -1/2
    soc = np.zeros( (len(states_ras) * 2 , len(states_ras) * 2), dtype=complex)
    
    n_soc_per_state = 0
    for i in range(0, len(states_ras) * 2, 2):
        for j in range(0, len(states_ras) * 2, 2):

            soc[i, j] = complex( soc_selected[n_soc_per_state, 0] , soc_selected[n_soc_per_state, 1] )
            n_soc_per_state = n_soc_per_state + 1

            soc[i, j + 1] = complex(float(soc_selected[n_soc_per_state, 0]), float(soc_selected[n_soc_per_state, 1]))
            n_soc_per_state = n_soc_per_state + 1

            soc[i + 1, j] = complex(float(soc_selected[n_soc_per_state, 0]), float(soc_selected[n_soc_per_state, 1]))
            n_soc_per_state = n_soc_per_state + 1

            soc[i + 1, j + 1] = complex(float(soc_selected[n_soc_per_state, 0]), float(soc_selected[n_soc_per_state, 1]))
            n_soc_per_state = n_soc_per_state + 1

    return soc


def get_spin_of_each_state(filename):

    search = [' RAS-CI total energy for state ']
    # state_num = []
    elements = []

    with open(filename) as file:
        for line in file:
            if any(i in line for i in search):

                next_line = next(file)
                next_line = next(file)

                element = next_line[16:22]
                elements.append(element.split())

    spin_each_state = np.array(elements, dtype=float)

    return spin_each_state


def get_spin_matrices(filename, states_ras, state_spin):
    """
    spin values are written in spin_matrix with < B(S,Sz) | in rows
    and | A(S',Sz') > in columns, and third dimension is the direction.

    For example, element spin_matrix[1,0,1] is < 1(1/2, +1/2) | Sy | 1(1/2, +1/2) >
    """
    
    elements = []
    search = ['Spin Matrices']
    search_spin_one_time = 0

    with open(filename) as file:
        for line in file:
            if any(i in line for i in search):
                next_line = next(file)
                next_line = next(file)

                if next_line.startswith("                     |Sz'=-0.50>") and (search_spin_one_time == 0):

                    # Save Sx matrix
                    next_line = next(file)
                    elements.append(next_line[13:35].split())
                    elements.append(next_line[37:61].split())
                    next_line = next(file)
                    elements.append(next_line[13:35].split())
                    elements.append(next_line[37:61].split())

                    # Save Sy matrix
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    elements.append(next_line[13:35].split())
                    elements.append(next_line[37:61].split())
                    next_line = next(file)
                    elements.append(next_line[13:35].split())
                    elements.append(next_line[37:61].split())

                    # Save Sz matrix
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    elements.append(next_line[13:35].split())
                    elements.append(next_line[37:61].split())
                    next_line = next(file)
                    elements.append(next_line[13:35].split())
                    elements.append(next_line[37:61].split())

                    search_spin_one_time = search_spin_one_time + 1

                if next_line.startswith("                     |Sz'=-1.50>") and (search_spin_one_time == 1):
                    
                    # Save Sx matrix
                    next_line = next(file)
                    next_line = next(file)
                    elements.append(next_line[37:61].split())
                    elements.append(next_line[63:87].split())
                    next_line = next(file)
                    elements.append(next_line[37:61].split())
                    elements.append(next_line[63:87].split())

                    # Save Sy matrix
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    elements.append(next_line[37:61].split())
                    elements.append(next_line[63:87].split())
                    next_line = next(file)
                    elements.append(next_line[37:61].split())
                    elements.append(next_line[63:87].split())

                    # Save Sz matrix
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    next_line = next(file)
                    elements.append(next_line[37:61].split())
                    elements.append(next_line[63:87].split())
                    next_line = next(file)
                    elements.append(next_line[37:61].split())
                    elements.append(next_line[63:87].split())
                    
                    break

    selected_spin = np.array(elements, dtype=float)
    spin_matrix = np.zeros((len(states_ras) * 2, len(states_ras) * 2, 3), dtype=complex)

    for n_state in states_ras:

        i = states_ras.index(n_state) * 2
        for n_dim in range(0, 3):

            if state_spin[n_state - 1] < 1: elem = 4 * n_dim  # Doublet cases
            elif 3 < state_spin[n_state - 1] < 4: elem = 12 + 4 * n_dim # Quartet cases

            spin_matrix[i, i, n_dim] = complex( selected_spin[elem, 0], selected_spin[elem, 1] )
            spin_matrix[i, i + 1, n_dim] = complex( selected_spin[elem+1, 0], selected_spin[elem+1, 1] )
            spin_matrix[i + 1, i, n_dim] = complex( selected_spin[elem+2, 0], selected_spin[elem+2, 1] )
            spin_matrix[i + 1, i + 1, n_dim] = complex( selected_spin[elem+3, 0], selected_spin[elem+3, 1] )

            # elif 3 < state_spin[n_state - 1] < 4:
            #
            #     spin_matrix[i, i, n_dim] = 0 + 0j
            #     spin_matrix[i, i + 1, n_dim] = 0 + 0j
            #     spin_matrix[i + 1, i, n_dim] = 0 + 0j
            #     spin_matrix[i + 1, i + 1, n_dim] = 0 + 0j

    return spin_matrix


def get_orbital_matrices(filename, totalstates, states_ras):
    """
    orbital angular momentum values are written in orbital_matrix with < B | in rows
    and | A > in columns, and third dimension is the direction.

    For example, element spin_matrix[2,0,0] is < state 2 | Lx | state 1 >
    """

    with open(filename) as f:
        data = f.readlines()

    searches = ['< B | Lx | A >', '< B | Ly | A >', '< B | Lz | A >' ] 
    elements = [] 

    for line in data:

        if any(i in line for i in searches):
            element = line[19:32]
            elements.append(element.split())

    angular_selection_list = []
    
    for j in states_ras:
        for i in states_ras:
            
            ndim = 0
            while ndim < 3:
                
                if (i == j):
                    angular_selection_list.append(['0.000000'])
                    
                elif i > j:
                    angular_position = 3 * (totalstates - 1) * (i - 1) + 3 * (j - 1) + ndim
                    angular_selection_list.append( elements[ angular_position ] )
                    
                elif i < j:
                    angular_position = 3 * (totalstates - 1) * (i - 1) + 3 * (j - 2) + ndim
                    angular_selection_list.append( elements[ angular_position ] )
                ndim = ndim + 1

    angular_selection = np.array(angular_selection_list, dtype=float)
    orbital_matrix = np.zeros( (len(states_ras)*2 ,len(states_ras)*2, 3), dtype=complex)

    t = 0
    for i in range(0,len(states_ras)*2,2):
        for j in range(0,len(states_ras)*2,2):
            for k in range(0,3):

                orbital_matrix[i, j, k] = complex(0, angular_selection[t + k])
                orbital_matrix[i+1, j+1, k] = complex(0, angular_selection[t + k])
            t = t + 3

    return orbital_matrix
