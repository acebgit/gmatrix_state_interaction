"""
    READ DATA OF Q-CHEM OUTPUTS
"""
import sys
import numpy as np

def get_number_of_states(file_ras):
    """
    Obtain the total number of states in RAS
    :param file_ras: output Q-Chem
    :return totalstates: number of states
    """
    with open(file_ras, encoding="utf8") as f:
        data = f.readlines()

    searches = ['Requested states: ']
    for line in data:
        if any(i in line for i in searches):
            line = line.split()
            element = line[3]
            break
    totalstates = int(element)
    return totalstates

def get_selected_states(input_ras, totalstates, symmetry_selection, states_ras, selected_states):
    """
    Select the states used depending on "selected_states" value:
    0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
    :param input_ras, totalstates, selected_states, symmetry_selection, states_ras
    :return states_ras
    """
    def get_states_same_symmetry(file_ras, totalstates, symmetry_selection):
        """
        Search the symmetry of each state and choose those with the selected symmetry
        :param file_ras, totalstates, symmetry_selection: output Q-Chem,
        number of states, symmetry selected
        :return states_selected_by_symmetry: symmetry of each state array
        """
        ground_state = 1  # Ground state to be included in symmetry calculation

        with open(file_ras, encoding="utf8") as f:
            data = f.readlines()

        searches = ['State symmetry: ']
        elements = []
        for line in data:
            if any(i in line for i in searches):
                line = line.split()
                element = line[2]
                element = element.replace('*', '') # '*' means mix of symmetries
                elements.append(element)

            if len(elements) == totalstates: break
        state_symmetries = np.array(elements)

        states_selected_by_symmetry = []
        states_selected_by_symmetry.append(ground_state) # GS is always added
        for nstate in range(1, totalstates + 1):
            if (state_symmetries[nstate - 1] == symmetry_selection) and (nstate != ground_state):
                states_selected_by_symmetry.append(nstate)

        if (states_selected_by_symmetry == [ground_state]):
            print('There is not this symmetry.')
            print('Change the symmetry selection:')
            for nstate in range(0, totalstates):
                print('- State', nstate+1, ', symmetry', state_symmetries[nstate])
            sys.exit()
        return states_selected_by_symmetry

    if (selected_states == 0):
        for i in states_ras:
            if (i <= 0 or i > totalstates):
                print("The number of states selected must be among the total number of states calculated in QChem.")
                print("Select a different number of states")
                sys.exit()

    elif (selected_states == 1):
        states_ras = list( range(1, totalstates+1) )

    elif (selected_states == 2):
        states_ras = get_states_same_symmetry(input_ras, totalstates, symmetry_selection)

    return states_ras

def get_eigenenergies(selected_states, input, theory_level, totalstates):
    """
    Get energies of the states: 0: RAS ; 1: CASSCF ; 2: SF-DFT
    :param selected_states, input_ras, theory_level
    :return: eigenenergies
    """
    with open(input, encoding="utf8") as f:
        data = f.readlines()

    if theory_level == 0:
        searches = [' RAS-CI total energy for state  ']
    elif theory_level == 1:
        searches = [' CAS-SCF            ']
    elif theory_level == 2:
        searches = [' Total energy for state ']

    elements = []

    for line in data:
        if any(i in line for i in searches):

            if theory_level == 0:
                line = line.split()
                element = line[6]
            elif theory_level == 1:
                element = line[29:47]
            elif theory_level == 2:
                element = line[43:59]
            elements.append(element)

        if (len(elements) == totalstates): break

    energies_selected = []
    for i in selected_states:
        energies_selected.append(elements[i - 1])
    eigenenergies = np.array(energies_selected, dtype=float)

    excitation_energies_eV = np.zeros(len(eigenenergies))
    for i in range(0, len(eigenenergies)):
        excitation_energies_eV[i] = (eigenenergies[i] - eigenenergies[0]) * 27.211399

    return eigenenergies, excitation_energies_eV

def get_spin_orbit_couplings(selected_SOC, totalstates, n_states, file_ras):
    """
    Spin-orbit coupling values are written in matrix with 'bra' in rows
    and 'ket' in columns, with spin order -1/2 , +1/2.
    :param selected_SOC, totalstates, n_states, file_ras: type of SOC
    considered, number of states, states selected, output Q-Chem
    :return SOC_matrix: SOC matrix
    """

    def s2_from_file(file, totalstates):
        """
        get s2 of each state from Q-Chem otuput
        :param s2: file
        :return: s2
        """
        search = ['  <S^2>      : ']
        elements = []

        with open(file, encoding="utf8") as file:
            for line in file:
                if any(i in line for i in search):

                    line = line.split()
                    element = float(line[2])
                    elements.append(element)
                if (len(elements) == totalstates): break

        s2_states = np.array(elements, dtype=int)
        return s2_states

    def get_total_number_SOC(s2_each_state):
        """
        get s2 of each state from Q-Chem otuput
        :param s2: file
        :return: s2
        """
        total_number_SOC = 0
        for s2 in s2_each_state:
            s = 0.5 * (-1 + np.sqrt(1 + 4 * s2))
            total_number_SOC += (2 * s + 1) * 2
        print(total_number_SOC)
        exit()
        return total_number_SOC

    if selected_SOC == 0:
        searches = 'Total mean-field SOC matrix'
    elif selected_SOC == 1:
        searches = ' 1-elec SOC matrix '
    elif selected_SOC == 2:
        searches = ' 2-elec mean-field SOC matrix '

    s2_each_state = s2_from_file(file_ras, totalstates)
    print(s2_each_state)
    total_number_SOC = get_total_number_SOC(s2_each_state)
    exit()

    with open(file_ras, encoding="utf8") as file:
        for line in file:
            if searches in line:

                next_line = next(file)

                if next_line.startswith("                     |Sz'=-0.50>"):
                    next_line = next(file)

                    while next_line.startswith(" <Sz="):
                        if next_line.startswith((" <Sz=-0.50|", " <Sz= 0.50|")):
                            # next_line = next_line.replace('***********i', '0.00000')
                            elements.append(next_line[12:35].split())  # add elements |Sz'=-0.50> to a list
                            elements.append(next_line[37:61].split())  # add elements |Sz'=+0.50> to a list
                            next_line = next(file)
                        else:
                            next_line = next(file)

                if next_line.startswith("                     |Sz'=-1.50>"):
                    next_line = next(file)

                    while next_line.startswith(" <Sz="):
                        if next_line.startswith((" <Sz=-0.50|", " <Sz= 0.50|")):
                            # next_line = next_line.replace('***********i', '0.00000')
                            elements.append(next_line[37:61].split()) # add elements |Sz'=-0.50> to a list
                            elements.append(next_line[63:87].split()) # add elements |Sz'=+0.50> to a list
                            next_line = next(file)
                        else:
                            next_line = next(file)

            if (len(elements) == total_number_SOC): break
    #--------------------------------------------------------------------------------
    soc_selected_list = []
    for ket in n_states:  # | A >
        for bra in n_states:  # < B |

            ket_index = 0
            while ket_index < 4:
                if (ket == bra): # SOC between same state values zero
                    soc_selected_list.append(['0.000000', '0.000000'])

                elif (ket < bra): # In Q-Chem, SOCs written when bra > ket
                    ket_position = 0
                    for ket_number in range(1, ket):
                        ket_position = ket_position + 4 * (totalstates - ket_number)
                    bra_position = 4 * (bra - ket - 1)

                    soc_position = ket_position + bra_position + ket_index
                    soc_selected_list.append( elements [ soc_position ] )

                else: pass

                ket_index = ket_index + 1
    soc_selected = np.array(soc_selected_list, dtype=float)
    #--------------------------------------------------------------------------------
    soc = np.zeros( (len(n_states) * 2 , len(n_states) * 2), dtype=complex)
    nel = 0

    for ket in n_states:  # |A, -1/2 >, |A, +1/2 >
        for bra in n_states:  # <B, -1/2|, <B, +1/2|
            ket_index = n_states.index(ket) * 2
            bra_index = n_states.index(bra) * 2

            if (ket == bra):
                soc[bra_index, ket_index] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                nel = nel + 1
                soc[bra_index, ket_index+1] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                nel = nel + 1
                soc[bra_index+1, ket_index] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                nel = nel + 1
                soc[bra_index+1, ket_index+1] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                nel = nel + 1

            if (ket < bra):  # In Q-Chem, SOCs written when ket < bra
                soc[bra_index, ket_index] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                soc[ket_index, bra_index] = np.conj(soc[bra_index, ket_index])
                nel = nel + 1

                soc[bra_index, ket_index+1] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                soc[ket_index+1, bra_index] = np.conj(soc[bra_index, ket_index+1])
                nel = nel + 1

                soc[bra_index+1, ket_index] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                soc[ket_index, bra_index+1] = np.conj(soc[bra_index+1, ket_index])
                nel = nel + 1

                soc[bra_index+1, ket_index+1] = complex(soc_selected[nel, 0], soc_selected[nel, 1])
                soc[ket_index+1, bra_index+1] = np.conj(soc[bra_index+1, ket_index+1])
                nel = nel + 1

    # print('\n'.join([''.join(['{:^20}'.format(item) for item in row])\
    #                  for row in np.round((soc[:,:]),5)]))
    # print(" ")
    # exit()

    soc = soc / 219474.63068  # From cm-1 to a.u.
    return soc

def get_spin_matrices(file, n_states):
    """
    Obtain the 3 dimensions of spin (Sx, Sy, Sz) from s2 of each state and
     put them in a 3-dim array. Spin is written with 'bra' in rows and
     'ket' in columns, with spin order -1/2 , +1/2.
    :param file, n_states: output Q-Chem, number of states
    :return total_spin_matrix: matrix with spins < m' | S | m >  in 3-D
    """
    def s2_from_file(file):
        """
        get s2 of each state from Q-Chem otuput
        :param s2: file
        :return: s2
        """
        search = ['  <S^2>      : ']
        elements = []

        with open(file, encoding="utf8") as file:
            for line in file:
                if any(i in line for i in search):

                    element = line[16:]
                    elements.append(element.split())

        s2_states = np.array(elements, dtype=float)
        return s2_states

    def s2_single_values(states_s2):
        """
        get s2 values from the s2 of all the states, to obtain all values of
        s2 only one time.
        :param s2: s2_states
        :return: s2_values
        """
        s2_values_list = []
        for i in range(0, len(states_s2)):
                if states_s2 [i] not in s2_values_list:
                    s2_values_list.append(states_s2 [i])

        s2_values = np.array(s2_values_list, dtype=float)
        return s2_values

    def s2_to_s(s2):
        """
        get total spin (s) from s2
        :param s2: s2
        :return: total spin (s)
        """
        return 0.5 * (-1 + np.sqrt(1 + 4 * s2))

    def s_to_s2(s):
        """
        get s2 from total spin (s)
        :param s: total spin (s)
        :return: s2
        """
        return s * (s + 1)

    # Abel function to determine the spin matrix:
    # https://github.com/abelcarreras/PyQchem/blob/1b1a0291f2737474955a5045ffbc56a2efd50911/pyqchem/tools/spin.py#L14

    def spin_matrices(s):
        """
        Get spin matrices Sx, Sy, Sz between two spin states (s,m) and (s,m') such that
        sx = < m' | Sx | m >, sy = < m' | Sy | m > and sz = < m' | Sz | m >
        :param s: total spin (s)
        :return: Sx, Sy, Sz
        """

        def are_equal(a, b, thresh=1e-4):
            return abs(a - b) < thresh

        def sz_values(s):
            return np.arange(-s, s + 1)

        # spin-multiplicities
        multiplicity = len(sz_values(s))

        # initialize Sx, Sy, Sz
        Sx = np.zeros((multiplicity, multiplicity), dtype=complex)
        Sy = np.zeros((multiplicity, multiplicity), dtype=complex)
        Sz = np.zeros((multiplicity, multiplicity), dtype=complex)

        # build spin matrices
        for i, sz_bra in enumerate(sz_values(s)):
            for j, sz_ket in enumerate(sz_values(s)):

                if are_equal(sz_bra, sz_ket):
                    Sz[i, j] = sz_ket

                if are_equal(sz_bra, sz_ket + 1):
                    Sx[i, j] = 0.5 * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)
                    Sy[i, j] = -0.5j * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)

                if are_equal(sz_bra, sz_ket - 1):
                    Sx[i, j] = 0.5 * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)
                    Sy[i, j] = 0.5j * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)

        return Sx, Sy, Sz

    s2_states = s2_from_file(file)
    single_s2_values = s2_single_values(s2_states)
    number_of_spins = len(single_s2_values)

    total_spin_matrix = np.zeros((len(n_states) * 2, len(n_states) * 2, 3), dtype=complex)

    for i in range(0, number_of_spins):
        for j in n_states:
            if s2_states[j - 1] == single_s2_values[i]:

                s = s2_to_s(single_s2_values[i])
                Sx, Sy, Sz = spin_matrices(s)

                s_dim = len(Sx) // 2 - 1
                total_dim = n_states.index(j) * 2

                total_spin_matrix[total_dim, total_dim, 0] = Sx[s_dim, s_dim]
                total_spin_matrix[total_dim+1, total_dim, 0] = Sx[s_dim+1, s_dim]
                total_spin_matrix[total_dim, total_dim+1, 0] = Sx[s_dim, s_dim+1]
                total_spin_matrix[total_dim+1, total_dim+1, 0] = Sx[s_dim+1, s_dim+1]
                
                total_spin_matrix[total_dim, total_dim, 1] = Sy[s_dim, s_dim]
                total_spin_matrix[total_dim+1, total_dim, 1] = Sy[s_dim+1, s_dim]
                total_spin_matrix[total_dim, total_dim+1, 1] = Sy[s_dim, s_dim+1]
                total_spin_matrix[total_dim+1, total_dim+1, 1] = Sy[s_dim+1, s_dim+1]
                
                total_spin_matrix[total_dim, total_dim, 2] = Sz[s_dim, s_dim]
                total_spin_matrix[total_dim+1, total_dim, 2] = Sz[s_dim+1, s_dim]
                total_spin_matrix[total_dim, total_dim+1, 2] = Sz[s_dim, s_dim+1]
                total_spin_matrix[total_dim+1, total_dim+1, 2] = Sz[s_dim+1, s_dim+1]
    return total_spin_matrix

def get_orbital_matrices(file_ras, totalstates, n_states):
    """
    Orbital angular momentum values are written in matrix with
    'bra' in rows and 'ket' in columns, with spin order -1/2 , +1/2.
    Third dimension is the direction.
    :param file_ras, totalstates, n_states: output Q-Chem, number of 
    states, states selected
    :return: orbital_matrix
    """
    searches = ['< B | Lx | A >', '< B | Ly | A >', '< B | Lz | A >' ]
    elements = []

    # take L values from output
    with open(file_ras, encoding="utf8") as f:
        data = f.readlines()

    for line in data:
        if any(i in line for i in searches):
            element = line[19:32]
            elements.append(element.split())
    # ----------------------------------------------------------------------------------
    # put orbital angular momentum in the array
    angular_selection_list = []

    for ket in n_states:  # | A >
        for bra in n_states:  # < B |

            n_dim = 0
            while n_dim < 3:

                if (ket == bra):  # momentum between same state values zero
                    angular_selection_list.append(['0.000000'])
                    # print('eq',ket,bra)

                elif (ket < bra): # In Q-Chem, L written when ket < bra
                    # print('diff', ket, bra)
                    ket_position = 0
                    for i in range(1, ket):
                        ket_position = ket_position + 3 * (totalstates - i)
                    bra_position = 3 * (bra - ket - 1)

                    orbital_position = ket_position + bra_position + n_dim
                    angular_selection_list.append(elements[orbital_position])

                n_dim = n_dim + 1

    angular_selection = np.array(angular_selection_list, dtype=float)
    # ----------------------------------------------------------------------------------------------
    orbital_matrix = np.zeros((len(n_states) * 2, len(n_states) * 2, 3), dtype=complex)
    nel = 0

    for ket in n_states:  # |A, -1/2 >, |A, +1/2 >
        for bra in n_states:  # <B, -1/2|, <B, +1/2|
            ket_index = n_states.index(ket) * 2
            bra_index = n_states.index(bra) * 2

            if (ket == bra):
                for k in range(0,3):
                        # print('eq', bra, ket, '-->', nel + k)
                        orbital_matrix[bra_index, ket_index, k] = complex(0, angular_selection[nel + k])
                        orbital_matrix[bra_index+1, ket_index+1, k] = complex(0, angular_selection[nel + k])

                        orbital_matrix[ket_index, bra_index, k] = (-1) * orbital_matrix[bra_index, ket_index, k]
                        orbital_matrix[ket_index+1, bra_index+1, k] = (-1) * orbital_matrix[bra_index+1, ket_index+1, k]
                nel = nel + 3

            if (ket < bra):  # In Q-Chem, SOCs written when ket < bra
                for k in range(0,3):
                        # print('diff', bra, ket, '-->', nel + k)
                        orbital_matrix[bra_index, ket_index, k] = complex(0, angular_selection[nel + k])
                        orbital_matrix[bra_index+1, ket_index+1, k] = complex(0, angular_selection[nel + k])

                        orbital_matrix[ket_index, bra_index, k] = (-1) * orbital_matrix[bra_index, ket_index, k]
                        orbital_matrix[ket_index+1, bra_index+1, k] = (-1) * orbital_matrix[bra_index+1, ket_index+1, k]
                nel = nel + 3

    # print('Angular momentums (x,y,z):')
    # for k in range(0,3):
    #    print('Dimension: ', k)
    #    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                     for row in np.round((orbital_matrix[:,:,k]),5)]))
    #    print(" ")
    return orbital_matrix
