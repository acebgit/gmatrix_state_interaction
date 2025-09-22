import json
import sys 
import numpy as np
from numpy import linalg, sqrt
import pandas as pd
import math 
from scipy import constants

import matplotlib
# matplotlib.use('Agg')  # Set backend to 'Agg' for headless environments
import matplotlib.pyplot as plt

# import timeit
# from functools import partial


def s2_to_s(s2):
    """
    get total spin (s) from s^2
    :param: s2_all
    :return: total spin (s)
    """
    return 0.5 * (-1 + np.sqrt(1 + 4 * s2))


def s_to_s2(spin):
    """
    get s2_all from total spin (s)
    :param: spin: total spin (s)
    :return: s2_all
    """
    return spin * (spin + 1)


def extract_data_from_json(filee):
    """
    Get lists with the information from "json" output filee.
    :param filee:
    :return: total_energy, excitation_energy_list, spin_list, soc_list, orbital_momentum_list
    """
    with open(filee, 'r') as f:
        object_text = f.read()
    outpuut_dict = json.loads(object_text)

    # total_energy_list = []
    # excitation_energy_list = []
    # for i in outpuut_dict['selected_energy_dict']:
    #     total_energy_list.append(float(outpuut_dict['selected_energy_dict'][i][0]))
    #     excit_energy = float(outpuut_dict['selected_energy_dict'][i][1]) 
    #     excitation_energy_list.append(excit_energy)

    # spin_list = [float(outpuut_dict['spin_dict'][i]) for i in outpuut_dict['spin_dict']]
    # soc_list = [outpuut_dict['soclist_dict'][i] for i in outpuut_dict['soclist_dict']]
    # orbital_momentum_list = [[complex(outpuut_dict['orbitalmomentlist_dict'][i][0]),
    #                         complex(outpuut_dict['orbitalmomentlist_dict'][i][1]),
    #                         complex(outpuut_dict['orbitalmomentlist_dict'][i][2])] 
    #                         for i in outpuut_dict['orbitalmomentlist_dict']]
    # transitions_list = [i for i in outpuut_dict['transitions_dict']]
    return outpuut_dict


def get_selected_dict(output__dict, totalstatess, initial__states, states__selection, sym__selected):
    """
    Depending on "states_option" it returns states: 
    0) in "state_ras" 
    1) all states selected 
    2) States selected by selected symmetry "symmetry_selection"
    :param: file, totalstates, selected_states, states_option, symmetry_selection
    :return: selected_states
    """
    selected_dict = {}

    if states__selection == 0:
        selected_states = initial__states

    elif states__selection == 1:
        selected_states = list(range(1, totalstatess+1))

    elif states__selection == 2:
        # Check that selected symmetry is between the symmetries
        mapping_symmetry = {key: index+1 for index, key in enumerate(output__dict["energy_dict"].keys())}
        if not any(sym__selected in key for key in mapping_symmetry.keys()):
            raise ValueError("The symmetry selected is not between the states.")
    
    # Create a dictionary with indices as order of the states
    # This is useful for mapping states symmetry in RASCI with their position
    mapping_symmetry = {index+1: key for index, key in enumerate(output__dict["energy_dict"].keys())}

    for name_dict, sub_dict in output__dict.items():
        # Check that initial states are between 1 and totalstates
        totalstates = len(output__dict["energy_dict"])
        if any(i <= 0 or i > totalstates for i in selected_states):
            raise ValueError("The number of states selected selected must be among the total number of states "
                                "selected calculated in QChem. ")

        data_dict = {}

        # Dictionaries scf_alpha_beta_electrons
        if name_dict == "scf_alpha_beta_electrons":
            data_dict = sub_dict

        # Dictionaries energy_dict, spin_dict, transitions_dict
        if name_dict == "energy_dict" or name_dict == "spin_dict":
            data_dict.update({i: sub_dict[mapping_symmetry[i]] for i in selected_states})
        
        # Dictionaries soc_matrix_dict, socc_dict, angmoment_dict
        if name_dict == "soc_matrix_dict" or name_dict == "socc_dict" or name_dict == "angmoment_dict": 
            for state_a in selected_states:
                for state_b in selected_states:
                    if state_a != state_b:
                        braket = mapping_symmetry[state_a] + "_" + mapping_symmetry[state_b]
                        if braket in sub_dict:
                            data_dict[str(state_a) + "_" + str(state_b)] = sub_dict[braket]

        # Dictionary transitions_dict
        if name_dict == "transitions_dict":
            data_list = []

            for transition in sub_dict:
                for state in selected_states:
                    if state == transition[0]['state']:
                        data_list.append(transition)
            data_dict = data_list

        selected_dict[name_dict] = data_dict
    return selected_dict


def get_soc_matrix(soc_dict, nstates, states_length_sz):
    """
    Construct the SOC matrix from the json filee dictionary. Result in eV. 
    Matrix is formed as: 
    < B, S, Sz | HSOC | A, S, Sz' >   -->   < row, sz_row | HSOC | column, sz_col >, 
    from -S to +S from left to right. 
    :param: soc_dict, nstates, states_length_sz
    :return: socmatrix
    """
    ndim = sum(states_length_sz.values())
    socmatrix = np.zeros((ndim, ndim), dtype=complex)
    
    initial_row = 0 
    for row_state in nstates:  # Iterate over bra < B | states (rows) 
        bra_size = states_length_sz[row_state] # Number of microstates in < B | 

        initial_column = 0 
        for col_state in nstates:  # Iterate over ket | A > states (columns) 
            ket_size = states_length_sz[col_state] # Number of microstates in | A > 

            # < B | larger than | A >, i.e. lower part triangle is calculated
            if row_state > col_state:

                # Create the interstate key as "A_B"
                interstate = str(col_state) + "_" + str(row_state)

                for sz_bra in range(bra_size): # Loop over Sz values for the bra state ( < B | )
                    for sz_ket in range(ket_size): # Loop over Sz values for the ket state ( | A > )
                        
                        # Compute matrix indices for row and column, adjusting offsets
                        matrix_row = initial_row + sz_bra   
                        matrix_col = initial_column + sz_ket 

                        # Assign the value to the SOC matrix and ensure Hermitian symmetry
                        socmatrix[matrix_row][matrix_col] = complex(soc_dict[interstate][sz_bra][sz_ket])
                        socmatrix[matrix_col][matrix_row] = np.conj(socmatrix[matrix_row][matrix_col]) # Conjugate symmetry
                        # print('row:', matrix_row, ', col:', matrix_col, ', value: ', value)

                initial_column += ket_size
            else: 
                initial_column += ket_size 
        initial_row += bra_size

    cm_to_ev = constants.physical_constants['inverse meter-electron volt relationship'][0] * 100
    socmatrix_ev = socmatrix * cm_to_ev
    return socmatrix_ev


def get_spin_matrices(approx_spin_dict, nstates, states_length_sz):
    """
    Get spin matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the order of
    "selected_state". Functions "s2_to_s", "s_to_s2" and "spin_matrices" are taken from inside PyQChem.
    Matrix is formed as: 
    < B, S, Sz | HSOC | A, S, Sz' >   -->   < row, sz1 | HSOC | column, sz2 >
    from -S to +S from left to right. 
    :param: filee, selected_state.
    :return: spin_matr, standard_spin_mat.
    """
    def spin_matrices(spin):
        """
        Get spin matrices s_x, s_y, s_z between two spin states selected (s,m) and (s,m') such that
        sxx = < m' | s_x | m >, syy = < m' | s_y | m > and szz = < m' | s_z | m >
        :param: spin: total spin (s)
        :return: s_x, s_y, s_z
        """

        def are_equal(a, b, thresh=1e-4):
            return abs(a - b) < thresh

        def sz_values(spinn):
            return np.arange(-spinn, spinn + 1)

        # spin-multiplicities
        multiplicity = len(sz_values(spin))

        # initialize s_x, s_y, s_z
        s_x = np.zeros((multiplicity, multiplicity), dtype=complex)
        s_y = np.zeros((multiplicity, multiplicity), dtype=complex)
        s_z = np.zeros((multiplicity, multiplicity), dtype=complex)

        # build spin matrices
        for iii, sz_bra in enumerate(sz_values(spin)):
            for g, sz_ket in enumerate(sz_values(spin)):

                if are_equal(sz_bra, sz_ket):
                    s_z[iii, g] = sz_ket

                if are_equal(sz_bra, sz_ket + 1):
                    s_x[iii, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[iii, g] = -0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)

                if are_equal(sz_bra, sz_ket - 1):
                    s_x[iii, g] = 0.5 * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
                    s_y[iii, g] = 0.5j * np.sqrt(s_to_s2(spin) - sz_bra * sz_ket)
        return s_x, s_y, s_z

    ndim = sum(states_length_sz.values())
    final_spin_matrix = np.zeros((ndim, ndim, 3), dtype=complex)

    # Form the spin matrix
    initial_position = 0
    for state in nstates:

        # Form the spin matrix of the state
        s_state = approx_spin_dict[state]
        sx, sy, sz = spin_matrices(s_state)
        s_xyz_dict = {0: sx, 1: sy, 2: sz}

        # Locate the matrix in the final spin matrix
        state_size = states_length_sz[state]

        for dim in range(3): 
            for sz_row in range(state_size):
                for sz_col in range(state_size):

                    # Compute matrix indices for row and column, adjusting offsets
                    matrix_row = initial_position + sz_row 
                    matrix_col = initial_position + sz_col 

                    # Assign the value to the final matrix 
                    final_spin_matrix[matrix_row][matrix_col][dim] = s_xyz_dict[dim][sz_row][sz_col]
        initial_position += state_size
    
    # Form the standard spin matrix
    ndim_standard = states_length_sz[nstates[0]]
    standard_spin_matrices = np.zeros((ndim_standard, ndim_standard, 3), dtype=complex)

    s_state = approx_spin_dict[nstates[0]]
    sx, sy, sz = spin_matrices(s_state)
    s_xyz_dict = {0: sx, 1: sy, 2: sz}

    for dim in range(3): 
        for sz_row in range(ndim_standard):
            for sz_col in range(ndim_standard):
                standard_spin_matrices[sz_row][sz_col][dim] = s_xyz_dict[dim][sz_row][sz_col]
    return final_spin_matrix, standard_spin_matrices


def get_orbital_matrices(angmoment_dict, nstates, states_length_sz):
    """
    Construct the L matrix from the L list of json filee. Matrix is formed as: 
    < B | L | A >   -->   < row | HSOC | col >
    :param: angmoment_dict, len_sz, nstates
    :return: all_multip_lk
    """
    # Form the orbital matrices, with dimensions nstates * nstates
    ndim = sum(states_length_sz.values())
    orbital_matrix = np.zeros((ndim, ndim, 3), dtype=complex)

    initial_column = 0 
    for col_state in nstates:  # Iterate over ket | A > states (columns) 
        ket_size = states_length_sz[col_state] # Number of microstates in | A > 

        initial_row = 0 
        for row_state in nstates:  # Iterate over bra < B | states (rows) 
            bra_size = states_length_sz[row_state] # Number of microstates in < B | 

            # < B | larger than | A >, i.e. lower part triangle is calculated
            if row_state > col_state:
                
                # Create the interstate key as "A_B"
                interstate = str(col_state) + "_" + str(row_state)

                # Create offsets from sz difference:
                shift_row = ((bra_size - ket_size) // 2) if bra_size > ket_size else 0
                shift_col = ((ket_size - bra_size) // 2) if ket_size > bra_size else 0
                len_sz = ket_size if bra_size > ket_size else bra_size

                # Loop over Sz values for the bra state (< B |)
                for sz_bra in range(len_sz):
                        
                    # Compute matrix indices for row and column, adjusting offsets
                    matrix_row = initial_row + shift_row + sz_bra 
                    matrix_col = initial_column + shift_col + sz_bra 

                    for k in range(0, 3):
                        # Convert the L value to a complex number
                        value = complex((angmoment_dict[interstate][k]))

                        # Assign the value to the L matrix and ensure Hermitian symmetry
                        try:
                            orbital_matrix[matrix_row][matrix_col][k] = value
                            orbital_matrix[matrix_col][matrix_row][k] = np.conj(value)  # Conjugate symmetry

                        except IndexError: 
                            print("Index error in orbital matrix")
                            print("- initial_row, shift_row, sz_ket:", initial_row, shift_row, sz_bra)
                            print("- initial_column, shift_col, sz_ket", initial_column, shift_col, sz_bra)

                        # print('row:', matrix_row, ', col:', matrix_col, ', value: ', value)
                initial_row += bra_size
            
            else: 
                initial_row += bra_size 
        initial_column += ket_size
    return orbital_matrix


def from_json_to_matrices(outpuut_dict_selected):
    """
    Transform all the dictionaries to lists or matrices 
    """
    def get_approx_spin_dict(spin_dict):
        """
        Form the spin list with approximated spins, since they can come from non-pure spin eigenfunctions.
        Form the sz list from -s to +s in +1 intervals for ground state and largest multiplicity state. 
        :return: sz__ground, sz__maxspin, approx__spin_list
        """ 
        # Round the minimum spin to the nearest multiple of 0.5
        minimum_s = s2_to_s(min(list(spin_dict.values())))
        approx_minimum_s = round(minimum_s * 2) / 2 
        approx_minimum_s_list = np.arange(approx_minimum_s, 15* approx_minimum_s + 1).tolist()

        # Form a dictionary with the approximated spins 
        approx_spin_dict = {}
        for k, state_s2 in spin_dict.items():
            state_approx_s = min(approx_minimum_s_list, key=lambda x: abs(x - s2_to_s(state_s2)))
            approx_spin_dict.update({k: state_approx_s})

        # Warning for singlet ground states
        if approx_spin_dict[next(iter(approx_spin_dict))] == 0:
            raise ValueError("WARNING! It is not allowed the calculation of the g-tensor in a singlet ground state.")
        return approx_spin_dict

    def count_intervals(data):
        """
        This function calculates how many values exist in the range [-value, value] 
        for each key in the dictionary, using a step size of 1.

        Parameters:
        data (dict): A dictionary where keys are identifiers and values are numerical limits.

        Returns:
        dict: A dictionary with the same keys but values representing the count of intervals.
        """
        result = {}  # Dictionary to store the results

        for key, value in data.items():
            # Use a fixed step size of 1
            step = 1  

            # Calculate the number of values in the range [-value, value] with step size 1
            count = int((value - (-value)) / step) + 1  

            # Store the result in the dictionary
            result[key] = count  

        return result  # Return the final dictionary

    def print_all_matrices(matrices_dict):
        print('-------------------------')
        print('INPUT MATRICES')
        print('-------------------------')

        cm_to_ev = constants.physical_constants['inverse meter-electron volt relationship'][0] * 100
        # matrices_dict["soc"] = matrices_dict["soc"] * cm_to_ev
        print('SOC:')
        print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
                        for row in np.round((matrices_dict["soc"] [:,:] / cm_to_ev),3)]))
        print()

        print('Spin matrices: ')
        for ndim in range(0,3):
            print('Ndim: ', ndim)
            print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
                            for row in ((matrices_dict["spin"][:,:,ndim]))]))
            print()
        
        print('STARDARD SPIN MATRICES: ')
        for ndim in range(0,3):
            print('Ndim: ', ndim)
            print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
                            for row in ((matrices_dict["spin"][:,:,ndim]))]))
            print()

        print('Orbital angular matrices: ')
        for ndim in range(0, 3):
            print('Dim: ', ndim)
            print('\n'.join([''.join(['{:^8}'.format(item) for item in row]) \
                            for row in np.round(matrices_dict["orbital"][:, :, ndim], 4)]))
            print()
            
    # Get a dictionary with the "approximated" spin, i.e. the real value
    approx_spins = get_approx_spin_dict(outpuut_dict_selected["spin_dict"])

    # Dimension coming from the multiplicity of each state 
    states_lengthsz = count_intervals(approx_spins)
    selected_states = list(states_lengthsz.keys())  # Take the list of states 
    
    # Transform the dictionaries to matrices
    excit_energies = {key: value[1] for key, value in outpuut_dict_selected["energy_dict"].items()}

    soc_matrix = get_soc_matrix(outpuut_dict_selected["soc_matrix_dict"], selected_states, states_lengthsz)

    spin_matrix, standard_spin_matrix = get_spin_matrices(approx_spins, selected_states, states_lengthsz)

    orbital_matrix = get_orbital_matrices(outpuut_dict_selected["angmoment_dict"], selected_states, states_lengthsz)

    matrices__dict = {
    "excit_energies": excit_energies,
    "soc": soc_matrix,
    "spin": spin_matrix, 
    "standard_spin": standard_spin_matrix,
    "orbital": orbital_matrix,
    }

    # print_all_matrices(matrices__dict)

    return states_lengthsz, approx_spins, matrices__dict


def select_soc_order(soc_order, selected__states, states__lengthsz, matrices__dict):
    """
    Take all the SOC matrix, the first-order or the higher-order terms
    """

    new_matrix = np.zeros((len(matrices__dict["soc"][:,:]), len(matrices__dict["soc"][:,:])), dtype=complex)

    if soc_order == 0:
        new_matrix = matrices__dict["soc"]

    elif soc_order == 1:
        for row in range(1, nstates):
            for column in range(0, row): 
                # Only for those socs involving ground state interaction
                if row == 0 or column == 0: 
                    for sz1 in range(0, len_sz):
                        for sz2 in range(0, len_sz):
                            matrix_row = row * len_sz + sz1 
                            matrix_col = column * len_sz + sz2 

                            new_matrix[matrix_row][matrix_col] = matrices__dict["soc"][matrix_row][matrix_col]
                            new_matrix[matrix_col][matrix_row] = matrices__dict["soc"][matrix_col][matrix_row]

    # 2ND ORDER DOES NOT MAKES SENSE IF WE ARE NOT INCLUDED GROUND STATE SOCS, 
    # SINCE G-TENSOR WILL BE ALWAYS ZERO
    elif soc_order == 2:
        for row in range(1, nstates):
            for column in range(0, row): 
                # Only for those socs not involving ground state interaction
                if row != 0 and column != 0: 
                    for sz1 in range(0, len_sz):
                        for sz2 in range(0, len_sz):
                            matrix_row = row * len_sz + sz1 
                            matrix_col = column * len_sz + sz2 

                            new_matrix[matrix_row][matrix_col] = matrices__dict["soc"][matrix_row][matrix_col]
                            new_matrix[matrix_col][matrix_row] = matrices__dict["soc"][matrix_col][matrix_row]
    
    matrices__dict["soc"] = new_matrix
    # print()
    # print('SOC:')
    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #                 for row in np.round((matrices__dict["soc"][:,:]/(constants.physical_constants['inverse meter-electron volt relationship'][0] * 100)))]))
    return matrices__dict


def hermitian_test(matrix, sz_list):
    """
    Check if matrix is Hermitian. If not, "ValueError".
    :param: matrix, sz_list
    """
    # Assuming matrix and sz_list are already defined
    matrix_rounded = np.round(matrix, 4)
    conjugate_matrix_rounded = np.round(np.conjugate(matrix.T), 4)

    # Create a mask for the upper triangle including the diagonal
    mask = np.triu(np.ones(matrix.shape, dtype=bool))

    # Compare elements
    different_elements = (matrix_rounded != conjugate_matrix_rounded) & mask

    # Get the indices where the elements are different
    i_indices, j_indices = np.where(different_elements)

    # Compute states
    state_1 = i_indices // len(sz_list)
    state_2 = j_indices // len(sz_list)

    # Combine states into a list of tuples if needed
    states = list(zip(state_1, state_2))

    if states: print(states)


def get_hamiltonian_construction(states__lengthsz, eigenenergies, soc_matrix):
    """
    Construct Hamiltonian matrix with dimensions 'bra' x 'ket', with spin order (-Ms , +Ms) in the order of
    "selected_states".
    Make hermitian test to the matrix.
    :param: selected_states, eigenenergies, spin_orbit_coupling, sz_values
    :return: hamiltonian
    """ 
    hamiltonian = np.zeros((len(soc_matrix), len(soc_matrix)), dtype=complex)
    nstates = list(states__lengthsz.keys())

    initial_row = 0
    for row in range(len(nstates)):  # Iterate over all states (rows)
        row_state = nstates[row]   # < B | (bra state)
        bra_size = states__lengthsz[row_state] # Number of elements in < B |

        # Loop over Sz values for the bra state ( < B | )
        for sz_bra in range(bra_size):
                
            # Compute matrix indices for row and column, adjusting offsets
            matrix_row = initial_row + sz_bra 

            # Assign the value to the hamiltonian matrix 
            hamiltonian[matrix_row][matrix_row] = eigenenergies[row_state]
        initial_row += bra_size

    # Fill the off-diagonal with spin_orbit_coupling values
    hamiltonian += soc_matrix
    return hamiltonian


def diagonalization(matrix):
    """
    Diagonalize Hamiltonian. Eigenvectors-eigenvalues are ordered by weight coefficients. 
    Construct the diagonal matrix.
    :param: matrix
    :return: eigenvalues, eigenvectors, diagonal_matrix
    """
    def reordering_eigenvectors(eigenval, eigenvect):
        """
        Reorder eigenvectors and eigenenergies by eigenvectors weight coefficients.
        :param: eigenvalues, eigenvectors
        :return: eigenvalues, eigenvectors
        """
        change_order = np.zeros(len(eigenvect), dtype=complex)
        for v_1 in range(0, len(eigenvect)):
            for v_2 in range(v_1, len(eigenvect)):
                if abs(eigenvect[v_1, v_2]) > abs(eigenvect[v_1, v_1]):
                    change_order[:] = eigenvect[:, v_1]
                    eigenvect[:, v_1] = eigenvect[:, v_2]
                    eigenvect[:, v_2] = change_order[:]
                    
                    change_order.real[0] = eigenval[v_1]
                    eigenval[v_1] = eigenval[v_2]
                    eigenval[v_2] = change_order.real[0]
        return eigenval, eigenvect

    try:
        eigenvalues, eigenvectors = linalg.eigh(matrix)
    except np.linalg.LinAlgError:
        print("Some values are NaN in the QChem output")
    
    eigenvalues, eigenvectors = reordering_eigenvectors(eigenvalues, eigenvectors)
    rotation_inverse = np.linalg.inv(eigenvectors)
    diagonal_matrix = np.matmul(np.matmul(rotation_inverse, matrix), eigenvectors)
    # for i in range(0, 3): #range(len(eigenvalues)):
    #     eigval = round(eigenvalues[i], 3)
    #     eigvec = [f"{x:.3f}" for x in eigenvectors[:, i]]
    #     print(f"Eigenvalue {i+1}: {eigval}")
    #     print(f"Eigenvector {i+1}: [{', '.join(eigvec)}]")
    #     print()
    # exit()
    return eigenvalues, eigenvectors, diagonal_matrix


def angular_matrices_obtaining(eigenvectors, ndim, input_angular_matrix):
    """
    Angular matrix from non-relativistic states (states from Q-Chem) is expanded in the relativistic states. In
    < B(S,Sz) | Sx | A(S',Sz') >, < B(S,Sz)| corresponds to rows and | A(S',Sz') > to columns.
    :param: eigenvectors, input_angular_matrix, sz_list
    :return: angular_matrix
    """
    angular_matrix = np.zeros((ndim, ndim, 3), dtype=complex)
    conj_eigenvectors = np.conjugate(eigenvectors)

    # eigenvectors[:][ket_final] first creates a copy of the entire eigenvectors 
    # array (due to eigenvectors[:]) and then attempts to access the ket_final-th 
    # element from this copy. This operation does not correctly retrieve the desired eigenvector.
    # print(eigenvectors[:, 0])
    # for bra_eigen2 in range(len(eigenvectors)):
    #     print(bra_eigen2, 0, eigenvectors[bra_eigen2, 0])
    # exit()

    # print('Kramer state eigenvectors:' )
    # for i in range(ndim):
    #     formatted_row = ["{:.2f}".format(value) for value in eigenvectors[:, i]]
    #     print(f"{i}: " + " ".join(formatted_row))
    # print()

    for k in range(0, 3): # For each dimension 
        for ket_final in range(ndim): # Iterate over ket | A > microstates (columns) 
            for bra_final in range(ndim): # Iterate over bra < B | microstates (rows) 
                # print('ROW:', bra_final, ', COL: ', ket_final)

                for ket_eigen in range(len(eigenvectors[:, ket_final])): # dimension of eigenvector
                    for bra_eigen in range(len(eigenvectors[:, bra_final])): # dimension of eigenvector
                        coeff_ket = (eigenvectors[ket_eigen, bra_final])
                        coeff_bra = (conj_eigenvectors[bra_eigen, ket_final])

                        angular_values = input_angular_matrix[bra_eigen, ket_eigen, k]

                        elements = coeff_bra * coeff_ket * angular_values
                        angular_matrix[bra_final, ket_final, k] += np.sum(elements)

                #         print('coeff row:', bra_eigen, np.round(coeff_bra, 2), 
                #               ', coeff col: ', ket_eigen, np.round(coeff_ket, 2), 
                #               ', angular value: ', angular_values, 
                #               ', final value: ', elements)
                #     print()
                # print()
                # exit()
    return angular_matrix


def projection_technique(standard_spin_matrix, s_matrix, l_matrix, ppms=0, lande_factor = 2.002319304363):
    """
    g-shift with orbitals and spin angular momentum matrices.
    :param: standard_spin_matrix, s_matrix, l_matrix, sz_list, ground_sz
    :return: g_shift
    """
    def j_diagonalization(matrix):
        """
        J-matrix diagonalization giving i) eigenvectors ii) diagonal matrix.
        :param: matrix
        :return: eigenvalues, eigenvectors, diagonal_matrix
        """
        try:
            eigenvalues, eigenvectors = linalg.eigh(matrix)
            rotation_inverse = np.linalg.inv(eigenvectors)
            diagonal_matrix = np.matmul(np.matmul(rotation_inverse, matrix), eigenvectors)
            return eigenvectors, diagonal_matrix
        except np.linalg.LinAlgError:
            print("Some values are NaN in the QChem output")

    def trace_g_values(total_j, spin_matr):
        """
        Obtaining g_value with projection of the traces.
        :param: total_j, spin_matr
        :return: g_value
        """
        a = np.matmul(total_j, spin_matr)
        b = np.matmul(spin_matr, spin_matr)
        g_value = np.trace(a) / np.trace(b)
        return g_value

    def units_gshift(gvalues, ppm=0):
        """
            Pass from ppt to ppm the gvalues.
            :param: ppm, gvalues
            :return: gvalues
            """    
        if ppm == 0:
            gvalues = [i * 1000 for i in gvalues]
        elif ppm == 1:
            gvalues = [i * 1000000 for i in gvalues]
        return gvalues

    j_matrix = lande_factor * s_matrix + l_matrix
    
    # PROJECTION TECHNIQUE TO OBTAIN THE TRIANGULAR G-MATRIX:
    # 1) g-value zz
    j_matrix_rotation, j_matrix_diagonal_z = j_diagonalization(j_matrix[:, :, 2])
    g_matrix_triangular = np.zeros((3, 3), dtype=complex)
    g_matrix_triangular[2, 2] = trace_g_values(j_matrix_diagonal_z, standard_spin_matrix[:, :, 2])

    # 2) pseudospin z
    pseudospin_matrix = np.zeros((len(j_matrix[0, :]), len(j_matrix[:, 0]), 3), dtype=complex)
    pseudospin_matrix[:, :, 2] = j_matrix[:, :, 2] / g_matrix_triangular[2, 2]

    # 3) g-value yz
    g_matrix_triangular[1, 2] = trace_g_values(j_matrix[:, :, 1], pseudospin_matrix[:, :, 2])
    residue = j_matrix[:, :, 1] - g_matrix_triangular[1, 2] * pseudospin_matrix[:, :, 2]

    # 4) g-value yy
    j_matrix_rotation, j_matrix_transformed_y = j_diagonalization(j_matrix[:, :, 1])
    g_matrix_triangular[1, 1] = trace_g_values(j_matrix_transformed_y, standard_spin_matrix[:, :, 2])

    # 5) pseudospin y
    pseudospin_matrix[:, :, 1] = residue[:, :] / g_matrix_triangular[1, 1]

    # 6) g-value xy and xz
    g_matrix_triangular[0, 1] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 1])
    g_matrix_triangular[0, 2] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 2])

    # 7) g-value xx
    j_matrix_rotation, j_matrix_transformed_x = j_diagonalization(j_matrix[:, :, 0])
    g_matrix_triangular[0, 0] = trace_g_values(j_matrix_transformed_x, standard_spin_matrix[:, :, 2])

    # 8) from g_matrix_triangular to g-shifts
    g_matrix = np.matmul(g_matrix_triangular, np.transpose(g_matrix_triangular))
    g_matrix_eigenvalues, rotation_g_matrix, g_matrix_diagonal = diagonalization(g_matrix)

    g_shifts = np.zeros(3, dtype=complex)
    for i in range(0, 3):
        g_shifts[i] = (sqrt(g_matrix_diagonal[i, i]) - lande_factor) 

    g_shifts = units_gshift(g_shifts, ppms)
    return g_matrix, g_shifts


def from_matrices_to_gshift(states_lengthsz, dict_matrices, ppms=0, lande_factor = 2.002319304363):
    """
    From matrices to the g-shift value. 
    """
    def print_all_matrices(hamiltonian, eigenvector, combination_spin_matrix, combination_orbital_matrix, gmatrix):
        print('-------------------------')
        print('PROJECTION TECHNIQUE')
        print('-------------------------')
        
        print('Hamiltonian: ')
        print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
                            for row in np.round((hamiltonian[:, :]), 5)])) 
        print()

        print('Kramer state eigenvectors:' )
        for i in range(len(combination_spin_matrix)):
            formatted_row = ["{:.2f}".format(value) for value in eigenvector[:, i]]
            print(f"{i}: " + " ".join(formatted_row))
        print()

        print('Linear combination of Spin angular matrix:')
        for k in range(0, 3):
            print('Dimension: ', k)
            print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
                            for row in np.round((combination_spin_matrix[:, :, k]), 5)]))
            print()
        
        print('Linear combination of Orbital angular matrix:')
        for k in range(0, 3):
            print('Dimension: ', k)
            print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
                            for row in np.round((combination_orbital_matrix[:, :, k]), 5)]))
            print()
        # exit()
        
        print('G-matrix: ')
        print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
                    for row in np.round((gmatrix[:,:]), 10)]))
        print()

    hamiltonian = get_hamiltonian_construction(states_lengthsz, dict_matrices["excit_energies"], dict_matrices["soc"])

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

    angular_matrix_dimension = next(iter(states_lengthsz.values()), None)

    combination_spin_matrix = angular_matrices_obtaining(eigenvector, angular_matrix_dimension, dict_matrices["spin"])

    combination_orbital_matrix = angular_matrices_obtaining(eigenvector, angular_matrix_dimension, dict_matrices["orbital"])

    gmatrix, gshift = projection_technique(dict_matrices["standard_spin"], combination_spin_matrix, combination_orbital_matrix, ppms, lande_factor)
    
    # print_all_matrices(hamiltonian, eigenvector, combination_spin_matrix, combination_orbital_matrix, gmatrix)

    return gmatrix, gshift 


def print_g_calculation(filee, approxspin_dict, g_shift, g_tensor, soc_order=0, soc_option=0, ppms=0):
    """
    Printing results. 
    """
    def count_spin_states(approx_spin_dict):
        """
        Obtain a dictionary with multiplicities and the states with those multiplicities. 
        """
        spin_mapping = {
            0: 'singlet',
            0.5: 'doublet',
            1: 'triplet',
            1.5: 'quartet',
            2: 'quintet',
            2.5: 'sextet',
            3: 'septet',
            3.5: 'octet',
            4: 'nonet',
            4.5: 'decet',
            5: 'undecet',
            5.5: 'duodecet',
            6: 'tredecet',
            6.5: 'quattuordecet',
            7: 'quindecet',
            7.5: 'sexdecet',
            8: 'septendecet',
            8.5: 'octodecet',
            9: 'novemdecet',
            9.5: 'vigintet',
            10: 'unvigintet'
        }

        data_dict = {key: [] for key in spin_mapping.values()}
        for state, spin_state in approx_spin_dict.items():
            data_dict[spin_mapping[spin_state]].append(state) 
        return data_dict

    print()
    print("---------------------")
    print(" G-SHIFT RESULTS")
    print("---------------------")
    print("File selected: ", filee.replace(".json", ""))

    soc_options_dict = {0: "All BP Hamiltonian", 1: "One electron term", 2: "Bielectornic term"}
    soc_order_dict = {0: "All orders", 1: "First-order", 2: "Second-order"}
    print("SOC: ", soc_order_dict[soc_order], ',', soc_options_dict[soc_option])
    print("States: ", list(approxspin_dict.keys()))

    multip_dict = count_spin_states(approxspin_dict)
    print("States multiplicity: ")
    for key, value in multip_dict.items():
        if value: 
            print(f"- {key} = {len(value)}: {value}")

    print(
        f"g-factor (x y z dimensions) in {'ppt' if ppms == 0 else 'ppm'}:", 
        np.round(g_shift[0].real, 3), 
        np.round(g_shift[1].real, 3), 
        np.round(g_shift[2].real, 3)
    )
    print('')

    # print("g_tensor: ")
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #             for row in np.round(np.sqrt(g_tensor[:,:]), 6)]))
    # print('')


def save_picture(save_options, filee, title_main):
    """
    Function that shows the plot (save_options=0) or save it (save_options=1).
    :param save_options, file, main_title.
    :return: plot (saved or shown)
    """
    if save_options == 1:
        plt.plot()
        figure_name = filee + '_' + title_main + '.png'
        plt.savefig(figure_name, 
                    bbox_inches="tight", 
                    dpi=600
                    )
        plt.close()
    else:
        plt.plot()
        plt.show()


def gshift_estimation_loop(output_dict, ppm):
        """
        Calculate the g-values in all the selected states.
        :param: output_dict, ppm
        :return: presentation_list
        """
        def estimation_phormula(orbitmoment, socc, energy):
            return abs(-4 * complex(orbitmoment) * socc / energy)
        
        ground_state = (list(output_dict["socc_dict"].keys())[0][0])
        from_cm_to_ev = 100 / constants.physical_constants['electron volt-inverse meter relationship'][0]

        g_xx = {}
        g_yy = {}
        g_zz = {}

        for interstate, socc in output_dict["socc_dict"].items():
            if interstate.startswith(ground_state+"_"):

                socc = socc * from_cm_to_ev
                excited_state = int(interstate.split("_")[1])
                excited_energy = output_dict["energy_dict"][excited_state][1]

                units = 10**3 if ppm == 0 else 10**6

                g_xx.update({excited_state: estimation_phormula(output_dict["angmoment_dict"][interstate][0], socc, excited_energy) * units})
                g_yy.update({excited_state: estimation_phormula(output_dict["angmoment_dict"][interstate][1], socc, excited_energy) * units})
                g_zz.update({excited_state: estimation_phormula(output_dict["angmoment_dict"][interstate][2], socc, excited_energy) * units}) 

        g_shift_dict = {"gxx": g_xx, "gyy": g_yy, "gzz": g_zz}
        return g_shift_dict


def gtensor_state_pairs_analysis(outputdict, ppms, cut_off_gvalue=0, cut_off_configs=0, plots=0, savepicture=0, cut_off_socc=0, cut_off_angmom=0):
    """
    Obtaining a matrix with several data for each excited state. 
    The cut-off determines the fraction of the amplitude of the 1st configuration that need to have the other configurations to be shown in each state.
    :param: file_ms_notnull, cutoff
    :return: excit_matrix
    """
    def get_new_active_space(active_space, orbitals, alpha_elec):
        """
        Return the orbitals to be added to the active space list
        :param orbitals:
        :param active_space:
        :return: active space
        """
        if orbitals == '-':  # It is a singlet, there is no SOMO so add HOMO
            active_space.append(alpha_elec)

        if type(orbitals) == int:  # It is a doublet, add SOMO
            active_space.append(int(orbitals))

        elif type(orbitals) == str and orbitals != '-':  # It is higher multiplicities, add all SOMOs
            orbitals = orbitals.split(',')
            for j in range(0, len(orbitals)):
                active_space.append(int(orbitals[j]))
        
        # Delete repeated orbitals
        active_space_set = set(active_space) 
        active_space = sorted(list(active_space_set))
        return active_space

    def get_new_active_space_electrons(active_space, alpha, beta):
        """
        Electrons in the new active space considering if they are occupied
        or unoccupied observing the HOMO position
        """
        electrons = 0
        for i in range(0, len(active_space)):
            if active_space[i] <= beta:
                electrons += 2
            elif (active_space[i] > beta) and (active_space[i] <= alpha):
                electrons += 1
            elif active_space[i] > alpha:
                pass
        return electrons

    def get_bar_chart(file, x_list, y_list, x_title, y_title, main_title, save_pict):
        """
        Print Bar plots: y_list vs x_list.
        :param:
        :return:
        """
        # "matplotlib" help: https://aprendeconalf.es/docencia/python/manual/matplotlib/
        # https://chartio.com/resources/tutorials/how-to-save-a-plot-to-a-file-using-matplotlib/
        # https://pythonspot.com/matplotlib-bar-chart/
        fuente = 'serif'  # "Sans"
        medium_size = 16
        big_size = 20

        y_pos = (list(x_list))
        plt.bar(x_list, y_list, align='center', width=0.5, color='r', edgecolor="black")
        plt.xticks(y_pos)

        plt.title(main_title, fontsize=big_size, fontname=fuente)
        plt.xlabel(x_title, fontsize=medium_size, fontfamily=fuente)
        plt.ylabel(y_title, fontsize=medium_size, fontfamily=fuente)
        # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.grid(True)

        save_picture(save_pict, file, main_title)

    def excited_states_plots(excit_list, save_pict):
            file_string = str(sys.argv[1])

            state_list = [excit_list[0][0]]
            energy_list = [excit_list[0][5]]
            socc_list = [0.0]
            momentum_list = [0.0]

            for line in excit_list[1:]:
                if line[0] not in state_list:
                    state_list.append(line[0])
                    energy_list.append(line[5])
                    socc_list.append(line[6])
                    momentum_list.append(max([float(num.replace('j', '')) for num in line[7].split(',')], key=abs))

            get_bar_chart(file_string[:-4], state_list, energy_list, 'Electronic State',
                        'Excitation energy (eV)', 'energ_analysis', save_pict)
            get_bar_chart(file_string[:-4], state_list, momentum_list, 'Electronic State',
                        'Orbital angular momentum', 'orbitmoment_analysis', save_pict)
            get_bar_chart(file_string[:-4], state_list, socc_list, 'Electronic State',
                        'SOCC (cm-1)', 'socc_analysis', save_pict)

    # Obtain all the g-shift estimations in the different direccions
    gshift_dict = gshift_estimation_loop(outputdict, ppms)

    # Cut-off for the g-values obtained
    cut_gxx = cut_off_gvalue * abs(max(gshift_dict["gxx"].values(), key=abs))
    cut_gyy = cut_off_gvalue * abs(max(gshift_dict["gyy"].values(), key=abs))
    cut_gzz = cut_off_gvalue * abs(max(gshift_dict["gzz"].values(), key=abs))

    # Cut-off for the soccs or orbital angular momentum
    cut_socc = cut_off_socc * max(list(outputdict["socc_dict"].values()), key=abs)
    grstate = outputdict["transitions_dict"][0][0]["state"]
    grstate_angmom = [value for key, value in outputdict["angmoment_dict"].items() if str(grstate)+"_" in key]
    cut_angmom = cut_off_angmom * abs(max([max([complex(x) for x in sublist], key=abs) for sublist in grstate_angmom], key=abs))

    # Obtain the final list with the estimated g-values
    # mapping_symmetry = {index+1: key for index, key in enumerate(outputdict["energy_dict"].keys())}
    gshift_presentation = [] # Final presentation list with g-shift estimation
    gshift_states = [] # Taken to be used in pandas
    final_active_orbitals = [] # New active space to be used

    excitstates_presentation = [] # Final presentation list with states analysis
    excitstates_states = [] # Taken to be used in pandas

    for k1, transit_state in enumerate(outputdict["transitions_dict"]):

        # Take the configurations with larger or equal to amplitude cut-off for each state
        max_amplitude_transition = max(outputdict["transitions_dict"][k1], key=lambda d: abs(d['amplitude']))
        cut_config = cut_off_configs * abs(max_amplitude_transition['amplitude'])

        for transit_configuration in outputdict["transitions_dict"][k1]:

            # For configurations with a certain amplitude
            if abs(transit_configuration["amplitude"]) >= cut_config:

                # Data required for g-tensor presentation list
                state = transit_configuration["state"]
                configuration = transit_configuration["configuration"]
                somos = transit_configuration["SOMO"]
                amplitude = transit_configuration["amplitude"]
                elec_alpha = outputdict["scf_alpha_beta_electrons"][1]

                # Data for excited states analysis
                energy = round(outputdict["energy_dict"][state][0], 2)
                excited_energy = round(outputdict["energy_dict"][state][1], 4)
                spin = round(outputdict["spin_dict"][state], 2)

                if k1 == 0: # For the ground state, where there are no g-value
                    gshift_presentation.append([state, configuration, amplitude, somos, "--", "--", "--"])
                    final_active_orbitals = get_new_active_space(final_active_orbitals, somos, elec_alpha)
                    gshift_states.append(state)

                    excitstates_presentation.append([state, spin, configuration, amplitude, somos, energy, excited_energy, "--", "--"])
                    excitstates_states.append(state)

                else: # For the excited states
                    # Including the g-values above a cut off
                    # interstate = f"{min(grstate, state)}_{max(grstate, state)}"
                    g_xx = abs(gshift_dict["gxx"][state])
                    g_yy = abs(gshift_dict["gyy"][state])
                    g_zz = abs(gshift_dict["gzz"][state])

                    threshold = 10**(-6)
                    if any(abs(g) >= cut and abs(cut) >= threshold for g, cut in zip((g_xx, g_yy, g_zz), (cut_gxx, cut_gyy, cut_gzz))):
                        gshift_presentation.append([str(state), str(configuration), amplitude, somos, 
                                                round(g_xx, 3), round(g_yy, 3), round(g_zz, 3)])
                        final_active_orbitals = get_new_active_space(final_active_orbitals, somos, elec_alpha)
                        gshift_states.append(state)
                    
                    socc = round(outputdict["socc_dict"][f"{grstate}_{state}"], 3)
                    angmom_max = max([complex(x) for x in outputdict["angmoment_dict"][f"{grstate}_{state}"]], key=abs)
                    
                    if socc >= cut_socc and abs(angmom_max) >= cut_angmom:
                        angmom = ','.join(str(complex(round(c.real, 2), round(c.imag, 2))) 
                                for c in map(complex, outputdict["angmoment_dict"][f"{grstate}_{state}"]))
                        excitstates_presentation.append([state, spin, configuration, amplitude, somos, energy, excited_energy, socc, angmom])
                        excitstates_states.append(state)

    print("----------------------")
    print(" G-TENSOR ANALYSIS")
    print("----------------------")
    print("cut-offs: ")
    print('- g-shift estimation (%):', cut_off_gvalue)
    print('- configurations (%): ', cut_off_configs)
    print("- somos to add to new active space: ", final_active_orbitals)
    print()
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame(np.array(gshift_presentation), index=gshift_states,
            columns=['state', 'config.', 'amplitude', 'somo', 'gxx', 'gyy', 'gzz'])    
    print(df.to_string(index=False))
    print()

    # electrons = get_new_active_space_electrons(final_active_orbitals, elec_alpha, elec_beta)
    # print('Final active space (HOMO =', elec_alpha, '):', '[', electrons, ',', len(final_active_orbitals), '] ;',
    #       final_active_orbitals)
    # print()

    print("------------------------")
    print(" EXCITED STATES ANALYSIS")
    print("------------------------")
    print("cut-offs: ")
    print('- configurations (%): ', cut_off_configs)
    print('- soccs (%): ', cut_off_socc)
    print('- angular momentum (%): ', cut_off_angmom)
    print()
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame(np.array(excitstates_presentation), index=excitstates_states,
            columns=['state', 'spin', 'config.', 'amplitude', 'somo', 'energy', 'excit. energy', 'socc', 'angmom'])    
    print(df.to_string(index=False))

    if plots == 1:
        excited_states_plots(excitstates_presentation, savepicture)


def plot_g_tensor_vs_states(file, subtitle, presentation_matrix, x_title, y_title, main_title, save_options):
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_type = 1 # 0: plot, 1: bars

    # MAIN FEATURES:
    fuente = 'sans-serif'  # 'serif'
    small_size = 17 # 25
    legend_size = small_size + 3
    bigger_size = small_size + 3
    weight_selected = 'normal'
    line_width = 2
    marker_size = 10

    x = presentation_matrix[:, 0]  # First column for x-axis
    y1 = presentation_matrix[:, 1]  # Second column for the first category
    y2 = presentation_matrix[:, 2]  # Third column for the second category
    y3 = presentation_matrix[:, 3]  # Fourth column for the third category

    #################################
    ###   PLOT TYPE
    #################################
    if plot_type == 0:
        # MAJOR AND MINOR TICKS:
        # x_tick = int((max(x))) / 4
        # x_tick = int((max(x))) / 4
        # y_tick = int((max(x))) / 4
        # x_tick = 20
        # y_tick = 1
        # ax.xaxis.set_major_locator(MultipleLocator(x_tick))
        # ax.yaxis.set_major_locator(MultipleLocator(y_tick))

        # x_tick_min = x_tick / 2
        # y_tick_min = y_tick / 2
        # ax.xaxis.set_minor_locator(MultipleLocator(x_tick_min))
        # ax.yaxis.set_minor_locator(MultipleLocator(y_tick_min))

        # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        # LINES:
        # ax.plot(x, y1, 'r',
        #         label=r'$\mathregular{\Delta g_{xx}}$', linewidth=line_width)
        # ax.plot(x, y2, 'b', 
        #         label=r'$\mathregular{\Delta g_{yy}}$', linewidth=line_width)
        # ax.plot(x, y3, 'k',
        #         label=r'$\mathregular{\Delta g_{zz}}$', linewidth=line_width)

        # MARKERS: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
        ax.plot(x, y1, 'ro', markersize=marker_size)
        ax.plot(x, y2, 'bv', markersize=marker_size,
                markerfacecolor='none', markeredgewidth=1.5)
        ax.plot(x, y3, 'ks', markersize=marker_size)

        # Enable grid lines for both x and y axes
        plt.grid(axis='both', linestyle='--', linewidth=0.7, alpha=0.7)

    #################################
    ###   BAR PLOTS
    #################################
    elif plot_type == 1:
        # Set width of the bars
        bar_width = 0.25

        # Set the positions of the bars on the x-axis
        r1 = x  # np.arange(len(x))
        r2 = [x + bar_width for x in r1]
        r3 = [x + bar_width * 2 for x in r1]

        # Create the bar plot 
        plt.bar(r1, y1, width=bar_width, color='red', edgecolor='red', label=r'$\mathregular{\Delta g_{xx}}$')
        plt.bar(r2, y2, width=bar_width, color='blue', edgecolor='blue', label=r'$\mathregular{\Delta g_{yy}}$')
        plt.bar(r3, y3, width=bar_width, color='black', edgecolor='black', label=r'$\mathregular{\Delta g_{zz}}$')

    # CHANGING THE FONTSIZE OF TICKS
    plt.xticks(fontsize=small_size, weight=weight_selected)
    plt.yticks(fontsize=small_size, weight=weight_selected)
    # axis.set_major_locator(MaxNLocator(integer=True))

    # LIMIT TO AXIS:
    ax.set_xlim(min(x)-1, max(x)+1)

    # To put only integer numbers in the label: 
    # from matplotlib.ticker import MaxNLocator
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    max_value = np.maximum.reduce([y1.max(), y2.max(), y3.max()])
    min_value = np.minimum.reduce([y1.min(), y2.min(), y3.min()])
    # Choose the larger value in absolute terms to calculate the interval
    interval = 0.1 * max(abs(max_value), abs(min_value))
    ax.set_ylim(min_value-interval, max_value+interval)

    # LABELS:
    # labelpad: change the space between axis umbers and labels
    plt.xlabel(x_title, fontsize=bigger_size, fontfamily=fuente, labelpad=15,
               weight=weight_selected)
    plt.ylabel(y_title, fontsize=bigger_size, fontfamily=fuente, style='italic',
               weight=weight_selected, labelpad=15)
    # x_min = 0
    # x_max =  11
    # y_min = -45
    # y_max =  5
    # plt.xlim([x_min, x_max])  # Limit axis values
    # plt.ylim([y_min, y_max])  # Limit axis values

    # TITLE:
    # y = 1.05 change the space between title and plot
    # plt.title(main_title, fontsize=bigger_size, fontfamily=fuente, y=1.05)
    
    # LEGEND
    legend = plt.legend(fontsize=legend_size, 
                        fancybox=True, 
                        framealpha=0.5,
                        labelcolor='linecolor', 
                        loc='best', 
                        frameon=False, 
                        )
    # plt.legend(frameon=False)
    # frame = legend.get_frame().set_edgecolor('black')
    # frame = legend.get_frame().set_linewidth(1)
    # frame = legend.get_frame().set_facecolor('white')
    # frame.set_edgecolor('black')

    # plt.locator_params(nbins=10)
    # plt.grid()

    # Add an horizontal line in y = 0
    # ax.hlines(y=0, xmin=x_min, xmax=x_max, linewidth=line_width, color='k',
    #           linestyle='dotted')
    # dotted, dashed, solid, dashdot

    # Frame of the plot: https://e2eml.school/matplotlib_framing.html#spinestyle
    line_width = line_width - 0.8
    ax.spines["top"].set_linewidth(line_width)
    ax.spines["bottom"].set_linewidth(line_width)
    ax.spines["left"].set_linewidth(line_width)
    ax.spines["right"].set_linewidth(line_width)
    
    save_picture(save_options, file, subtitle)


def sum_over_state_plot(outputdict, gestimation, ppm, cutoff, saveplot):
    """
    Generate the sum-over-states plot, i.e. calculation of the g-tensor by including states
    from 1 to nstates, above a cutoff of g-value. 
    This can be done by:
    - gestimation = 0: effective Hamiltonian created between each pair of states
    - gestimation = 1: use of predictive phormula 
    :param: 
    :return: shows SOS plot
    """
    def filter_dictionary(dictionary, state):
        """
        Form a dictionary with only a pair of states: ground state and "state"
        """
        new_dict = {}
        ground_state = list(outputdict["energy_dict"].keys())[0]
        
        for name_dict in ["energy_dict", "spin_dict"]:
            new_dict[name_dict] = {k: v for k, v in dictionary[name_dict].items() if k == ground_state or k == state}

        for name_dict in ["soc_matrix_dict", "angmoment_dict"]:
            for k, v in dictionary[name_dict].items():
                if k in [f"{ground_state}_{state}", f"{state}_{ground_state}"]:
                    new_dict[name_dict] = {k: v}
        return new_dict
    
    filtered_gshifts = []

    if gestimation == 0:
        all_gshifts = []
        
        for excit_state in list(outputdict["energy_dict"].keys())[1:]:
            
            # Form a dictionary only for a pair of states
            filtered_dict = filter_dictionary(outputdict, excit_state)
            
            states__lengthsz, approxspin_dict, matrices_dict = from_json_to_matrices(filtered_dict)

            gmatrix, gshift = from_matrices_to_gshift(states__lengthsz, matrices_dict, ppm)
            
            all_gshifts.append([excit_state, 
                                (np.round(gshift[0].real, 3)),
                                (np.round(gshift[1].real, 3)), 
                                (np.round(gshift[2].real, 3))])

        # Convert to a NumPy array for efficient processing
        all_gshifts_array = np.array(all_gshifts, dtype=np.float64)

        # Step 1: Find the three maximum values in each of the last three columns and multiply by cutoff
        cutoffs = np.max(np.abs(all_gshifts_array[:, 1:]), axis=0) * cutoff

        # Step 2: Filter rows where at least one column meets or exceeds the threshold
        data = [
            [int(row[0])] + list(row[1:])  # Convert the first value to an integer and keep the rest as float
            for row in all_gshifts_array
            if any(abs(row[i]) >= cutoffs[i-1] for i in range(1, 4))
        ]

        # Convert np.float64 to regular float
        filtered_gshifts = [
            [row[0]] + [float(value) for value in row[1:]]  # Convert all elements except the first to floats
            for row in data
        ]

    elif gestimation == 1:
        gshift_dict = gshift_estimation_loop(outputdict, ppm)

        # Avoid to include g-tensor with zero
        threshold = 10**(-6)
        new_cutoff = cutoff if cutoff != 0 else threshold

        # Create the new dictionary with the maximum absolute value multiplied by cutoff
        cut_gvalues = {key: max((abs(v), v) for v in values.values())[1] * new_cutoff for key, values in gshift_dict.items()}
        
        # Take states with estimated g-shift higher than a cutoff
        for k, v in gshift_dict["gxx"].items():
            if any(abs(gshift_dict[key][k]) >= cut_gvalues[key] and cut_gvalues[key] >= threshold
                for key in ["gxx", "gyy", "gzz"]):
                    filtered_gshifts.append([int(k),
                                             (gshift_dict["gxx"][k]),
                                             (gshift_dict["gyy"][k]),
                                             (gshift_dict["gzz"][k])])

    print("------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("------------------------------")
    technique = {0: "Effective Hamiltonian", 1: "Estimation phormula"}
    print("Technique used:", technique.get(gestimation, "Unknown"))
    print("cut-offs g-value (%): ", cutoff)

    # Set display options to show all rows and columns
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame([row[0:4] for row in filtered_gshifts], [row[0] for row in filtered_gshifts], columns=['state','gxx','gyy','gzz'])
    print(df.to_string(index=False)) 

    # Compute the sum of all the elements
    perturbative_sum = [round(float(sum(x)), 3) for x in zip(*filtered_gshifts)][1:]  # Sum columns, skipping the first
    print('Total: ', perturbative_sum[0], perturbative_sum[1], perturbative_sum[2])

    y_title = r'$\Delta g, ppt$' if ppm == 0 else r'$\Delta g, ppm$'
    file_string = str(sys.argv[1]).split('.')[0]
    plot_title = 'sos_analysis: ' + file_string

    plot_g_tensor_vs_states(file_string,
                                'sos', 
                                np.array(filtered_gshifts, dtype=object),
                                '# roots',
                                y_title,
                                plot_title, 
                                saveplot)
    return filtered_gshifts


def filter_list(my_list, ncolumn, cutoff):
    """
    Get a list filtering the rows in "my_list" by the value "cutoff" of the column "ncolumn"
    arg: my_list, ncolumn, cutoff
    """
    filtered_list = []
    all_states = list(set([sublist[0] for sublist in my_list]))
    for state in all_states:
        configurations = [row for row in my_list if row[0] == state]
        max_value = max(abs(sublist[ncolumn]) for sublist in configurations)
        for config in configurations:
            if abs(config[ncolumn]) >= cutoff * max_value:
                filtered_list.append(config)
    return filtered_list
    

def measure_function_time(func, *args, **kwargs): 
    funct_and_args = partial(func, *args, **kwargs)
    tiempo = timeit.timeit(funct_and_args, number=1)
    nombre_funcion = func.__name__
    print(f"{nombre_funcion} tardó {tiempo:.4f} segundos en ejecutarse.")


def excitedstates_analysis(nstate, excitenergies, orbitmoment, soc, plot, cutoff_amp):
    """
    Get a table in pandas with information for each of the excited states. 
    """
    # max_orbitmoment = [abs(max(sublista, key=abs)) for sublista in orbitmoment]
    socc = [((sum((abs(complex(elemento)))**2 for sublista in i for elemento in sublista))**(1/2)) for i in soc]
    
    # Make a list of orbital angular momentums rounded by two decimals 
    np_data = np.array(orbitmoment)
    rounded_data = np.around(np_data, decimals=2)
    orbitmoment = rounded_data.tolist()

    presentation_list_all = []
    list_states = []
    for state in range(0, nstate):
        for transition in list(transitions_json[state]):
            if state == 0:
                presentation_list_all.append([state+1, excitenergies[state], "---", "---", 
                                        transition['transition/SOMO'], 
                                        transition['amplitude']])
            else:
                presentation_list_all.append([state+1, np.round(excitenergies[state],3), 
                                        orbitmoment[state-1], 
                                        np.round(socc[state-1],3),  
                                        transition['transition/SOMO'], 
                                        transition['amplitude']])
    
    # Filter all configurations by those with largest amplitude
    presentation_list = []
    presentation_list = filter_list(presentation_list_all, 5, cutoff_amp)

    print("--------------------------")
    print(" EXCITED STATES ANALYSIS")
    print("--------------------------")
    # Set display options to show all rows and columns
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame(presentation_list, columns=['state', 'energy','orbit mom','socc','transition/SOMO','amplitude'])
    print(df.to_string(index=False))
    print()

    if plot == 1:
        get_bar_chart(file[:-4], [i for i in range(1,nstate+1)], excitenergies, 'Electronic State',
                        'Excitation energy (eV)', 'energ_analysis', save_pict=0)
        orbit_max_values = [max(map(abs, sublist)) for sublist in orbitmoment]
        get_bar_chart(file[:-4], [i for i in range(2,nstate+1)], orbit_max_values, 'Electronic State',
                        'Orbital angular momentum', 'orbit_analysis', save_pict=0)
        get_bar_chart(file[:-4], [i for i in range(2,nstate+1)], socc, 'Electronic State',
                        'Spin-orbit coupling constants (cm-1)', 'soc_analysis', save_pict=0)


def from_gvalue_to_shift(lista):
    """
    Obtain the g-shifts from the g-values
    :param: list
    :return: g_shift
    """
    lande_factor = 2.002319304363
    g_shift = []
    for i in range(0, len(lista)):
        value = (lista[i] - lande_factor) * 10**3
        g_shift.append(value)
    print(np.round(g_shift, 2))


def comparison_s2(json_file, outputdict, save_options, s2__per_gshift, gshift__list):
    """
    Program to obtain the plot comparing the calculated <S2> saved in dictionaries and 
    the expected value <S2>. 
    """
    def plot_spin_states(file_name, presentation_matrix, x_title, y_title, main_title, save_options):
        fig, ax = plt.subplots(figsize=(10, 5))

        # MAIN FEATURES:
        fuente = 'sans-serif'  # 'serif'
        small_size = 17
        legend_size = small_size + 3
        bigger_size = small_size + 3
        weight_selected = 'normal'
        line_width = 2
        marker_size = 10

        x = presentation_matrix[:, 0]  # First column for x-axis
        y1 = presentation_matrix[:, 1]  # Second column for the first category
        y2 = presentation_matrix[:, 2]  # Third column for the second category

        # Set width of the bars
        bar_width = 0.25

        # Set the positions of the bars on the x-axis
        r1 = x  # np.arange(len(x))
        r2 = [x + bar_width for x in r1]
        r3 = [x + bar_width * 2 for x in r1]

        # Create the bar plot
        plt.bar(r1, y1, width=bar_width, color='red', edgecolor='red', label=r'Calc. $\mathregular{<S^{2}>}$')
        plt.bar(r2, y2, width=bar_width, color='blue', edgecolor='blue', label=r'Real $\mathregular{<S^{2}>}$')

        # CHANGING THE FONTSIZE OF TICKS
        plt.xticks(fontsize=small_size, weight=weight_selected)
        plt.yticks(fontsize=small_size, weight=weight_selected)
        # axis.set_major_locator(MaxNLocator(integer=True))

        # LIMIT TO AXIS:
        ax.set_xlim(min(x)-0.5, max(x)+1)
        max_value = np.maximum.reduce([y1.max(), y2.max()])
        min_value = np.minimum.reduce([y1.min(), y2.min()])
        ax.set_ylim(min_value-0.05*max_value, max_value+0.1*max_value)

        # LABELS:
        # labelpad: change the space between axis umbers and labels
        plt.xlabel(x_title, fontsize=bigger_size, fontfamily=fuente, labelpad=15,
                weight=weight_selected)
        plt.ylabel(y_title, fontsize=bigger_size, fontfamily=fuente, style='italic',
                weight=weight_selected, labelpad=15)
        # x_min = 0
        # x_max =  11
        # y_min = -45
        # y_max =  5
        # plt.xlim([x_min, x_max])  # Limit axis values
        # plt.ylim([y_min, y_max])  # Limit axis values

        # TITLE:
        # y = 1.05 change the space between title and plot
        # plt.title(main_title, fontsize=bigger_size, fontfamily=fuente, y=1.05)
        
        # LEGEND
        legend = plt.legend(fontsize=legend_size, 
                            fancybox=True, 
                            framealpha=0.5,
                            labelcolor='linecolor', 
                            loc='best', 
                            frameon=False, 
                            )
        # plt.legend(frameon=False)
        # frame = legend.get_frame().set_edgecolor('black')
        # frame = legend.get_frame().set_linewidth(1)
        # frame = legend.get_frame().set_facecolor('white')
        # frame.set_edgecolor('black')

        # plt.locator_params(nbins=10)
        # plt.grid()

        # Add an horizontal line in y = 0
        # ax.hlines(y=0, xmin=x_min, xmax=x_max, linewidth=line_width, color='k',
        #           linestyle='dotted')
        # dotted, dashed, solid, dashdot

        # Frame of the plot: https://e2eml.school/matplotlib_framing.html#spinestyle
        line_width = line_width - 0.8
        ax.spines["top"].set_linewidth(line_width)
        ax.spines["bottom"].set_linewidth(line_width)
        ax.spines["left"].set_linewidth(line_width)
        ax.spines["right"].set_linewidth(line_width)
        
        save_picture(save_options, file_name, 's2')

    # Take the real spins
    real_s2 = outputdict["spin_dict"]

    # If the spins come from even or odd multiplicity and 
    # create and approximate spins list for this. This is done because in some cases, 
    # triplets have S2 value that is closer to doublet than to triplet 
    ground_mult = int(2*(outputdict["spin_dict"][1])+1)
    if ground_mult % 2 == 0: # doublets, quartets...
        approx_s = [float(i) for i in np.arange(0.5, 5.5, 1)]
    else: # singlets, triplets...
        approx_s = [float(i) for i in np.arange(0, 5, 1)]

    # Obtain the plots with S2 real and approximated
    s_diff = {}
    final_table = []

    for k, v in real_s2.items():
        # From s2 to s
        s = s2_to_s(v)

        # Find the closest s and convert it to s2 again
        closest_s = min(approx_s, key=lambda x: abs(x - s))
        closest_s2 = s_to_s2(closest_s)

        # Take the value to a table 
        final_table.append([int(k), 
                            round(v, 3), 
                            round(closest_s2, 3), 
                            round(v-closest_s2, 3)])
    
    print('')
    print("------------------------------")
    print(" States <S2> comparison")
    print("------------------------------")

    # Set display options to show all rows and columns
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame(final_table, columns=['State','Calc. <S2>','Real <S2>', 'error'])
    print(df.to_string(index=False))
    print()

    if s2__per_gshift == 0:
        # Plot the spin relative errors for each of the states 
        plot_spin_states(json_file.replace(".out.json", ""), 
                        np.array(final_table, dtype=object), 
                        '# roots', 
                        '$\mathregular{<S^{2}>}$', 
                        '<S2> comparison', 
                        save_options)
    
    if s2__per_gshift == 1:
        # Compute the s2 contamination * gvalue 
        s2diff_list = []
        for i in range(0, len(gshift__list)):
            state = gshift__list[i][0]  
            xx = (final_table[i+1][1]-final_table[i+1][2]) * abs(gshift__list[i][1])
            yy = (final_table[i+1][1]-final_table[i+1][2]) * abs(gshift__list[i][2])
            zz = (final_table[i+1][1]-final_table[i+1][2]) * abs(gshift__list[i][3])            
            s2diff_list.append([state, round(xx, 3), round(yy, 3), round(zz, 3)])
        
        # Set display options to show all rows and columns
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        df = pd.DataFrame(s2diff_list, columns=['State', 'S^2 * g_{xx}', 'S^2 * g_{yy}', 'S^2 * g_{zz}'])
        print(df.to_string(index=False))
        print()

        # y_title = r'$(S^2_{calc} - S^2_{real})/\Delta g, ppt^{-1}$' # if ppm == 0 else r'$\Delta g, ppm$'
        y_title = r'$\sigma_{k}, ppt$'
        file_string = str(sys.argv[1]).split('.')[0]
        plot_title = 'sos_analysis: ' + file_string

        plot_g_tensor_vs_states(file_string,
                                's2cont', 
                                np.array(s2diff_list, dtype=object),
                                '# roots',
                                y_title,
                                plot_title, 
                                save_options)        


def remove_duplicate_positions(data_dict):
    """
    Checks if the second element (value at key index 1) contains values that differ by 2 or less.
    If such a pair is found, removes only the second occurrence in each pair from all lists.

    Parameters:
        data_dict (dict): Dictionary where keys map to lists of values.

    Returns:
        dict: A new dictionary with selected close-value positions removed from all lists.
    """
    keys = list(data_dict.keys())

    # Ensure there is a second entry in the dictionary
    if len(keys) < 2:
        return data_dict

    second_key = keys[1]
    second_list = data_dict[second_key]

    remove_indices = set()
    used = set()

    for i in range(len(second_list)):
        for j in range(i + 1, len(second_list)):
            if abs(second_list[i] - second_list[j]) <= 2 and j not in remove_indices and j not in used:
                remove_indices.add(j)
                used.add(i)
                break  # remove only one element per conflict

    if not remove_indices:
        return data_dict

    # Remove selected positions from all lists
    new_dict = {}
    for key, values in data_dict.items():
        new_dict[key] = [v for idx, v in enumerate(values) if idx not in remove_indices]

    return new_dict



def fit_polynomial(scaling_dict, degree=3):
    """
    Fits a polynomial to each of the three elements in the dictionary and returns the coefficients and R² values.
    the first element will correspond to the coefficient of x², the second element will correspond to the coefficient of x,
    and the third element will correspond to the constant term.
    
    Parameters:
    scaling_dict (dict): Dictionary where keys are x-values and values are lists of three y-values.
    degree (int): Degree of the polynomial fit.
    
    Returns:
    dict: A dictionary with polynomial coefficients and R² values for each y-series.
    """
    x_values = np.array(list(scaling_dict.keys()))
    y_values = np.array(list(scaling_dict.values()))

    poly_results = {}

    for i, key in enumerate(["gxx", "gyy", "gzz"]):
        coeffs = np.polyfit(x_values, y_values[:, i], degree)  # Fit polynomial
        p = np.poly1d(coeffs)  # Create polynomial function
        
        # Calculate R²
        y_pred = p(x_values)  # Predicted values
        ss_res = np.sum((y_values[:, i] - y_pred) ** 2)  # Residual sum of squares
        ss_tot = np.sum((y_values[:, i] - np.mean(y_values[:, i])) ** 2)  # Total sum of squares
        r2 = 1 - (ss_res / ss_tot)  # R² calculation
        
        poly_results[key] = {"coefficients": coeffs, "R^2": r2}  # Store results
    
    return poly_results


def plot_and_fit_scaling_dict(scaling_dict, name, xname, save_pic, degree=2):
    """
    Plots each of the three elements in the list against the key in the dictionary and fits a polynomial.
    The polynomial equation is displayed on the graph.
    
    Parameters:
    scaling_dict (dict): Dictionary where keys are x-values and values are lists of three y-values.
    degree (int): Degree of the polynomial fit.
    """
    # Letter size
    smalllet = 10

    x_values = np.array(list(scaling_dict.keys()))
    y_values = np.array(list(scaling_dict.values()))

    fig, ax = plt.subplots(figsize=(8, 6))
    
    labels = ["$\mathregular{\Delta g_{\perp}}$", "$\mathregular{\Delta g_{\parallel}}$", "y3"]
    markers = ['o', 's', '^']
    colors = ['r', 'b', 'b']
    
    for i in range(2):
        # Scatter plot
        ax.scatter(x_values, y_values[:, i], label=labels[i], marker=markers[i], color=colors[i])

        # Polynomial fit
        coeffs = np.polyfit(x_values, y_values[:, i], degree)
        poly_eq = np.poly1d(coeffs)

        # Generate smooth curve for fitting
        x_fit = np.linspace(min(x_values), max(x_values), 100)
        y_fit = poly_eq(x_fit)
        ax.plot(x_fit, y_fit, color=colors[i], linestyle="--", label="_nolegend_") # label=f"{labels[i]} fit"

        # Display equation on the plot
        # eq_str = " + ".join([f"{c:.3f}x^{degree - j}" if j < degree else f"{c:.3f}" 
        #                     for j, c in enumerate(coeffs)])
        # ax.text(0.05, 0.9 - i * 0.1, f"{labels[i]}: {eq_str}", 
        #         transform=ax.transAxes, fontsize=10, color=colors[i])
    
    ax.tick_params(axis='both', which='major', labelsize=smalllet + 4)

    ax.set_xlabel(xname, fontsize=smalllet + 6)
    ax.set_ylabel("$\mathregular{\Delta g}$, ppt", fontsize=smalllet + 6)
    ax.legend(fontsize=smalllet + 8)
    # ax.grid(True) 
    
    # Showing the plot 
    save_picture(save_pic, name, '')


def get_scaling_analysis(scaling__analysis, degree_fit, outpuut_dict, file, savepic, ppms):
    """
    SOC scaling analysis. 
    """
    states__lengthsz, approxspin_dict, matrices_dict = from_json_to_matrices(outpuut_dict)
    soc = matrices_dict["soc"].copy()
    orbmoment = matrices_dict["orbital"].copy()
    # spin_mat = matrices_dict["spin"].copy()

    a = 0
    b = 1.6
    interval = 0.2

    if scaling__analysis == 1:
        print("------------------------------")
        print(" SOC SCALING ANALYSIS")
        print("------------------------------")

        scaling_dict = {}
        for soc_factor in np.arange(a, b, interval):
            matrices_dict["soc"] = soc * soc_factor

            gmatrix, gshift = from_matrices_to_gshift(states__lengthsz, matrices_dict, ppms)

            scaling_dict.update({round(float(soc_factor), 1): [round(float(gshift[0].real), 3), 
                            round(float(gshift[1].real), 3), 
                            round(float(gshift[2].real), 3)]})
            
            print('Scaling factor: ', soc_factor, ', g-shifts: ', [round(float(gshift[0].real), 3), 
                            round(float(gshift[1].real), 3), 
                            round(float(gshift[2].real), 3)])
        # print()
        
        # Removing the second element from each value list
        new_scaling_dict = remove_duplicate_positions(scaling_dict)
        namee = file.split('_')[0]+"_soc_scaling"
        plot_and_fit_scaling_dict(new_scaling_dict, namee, "SOC Factor Scaling", savepic, degree_fit)

        # Print formatted output
        polynomial_coeffs = fit_polynomial(scaling_dict, degree_fit)
        print("\nSOC scaling - Polynomial fit (coefficient A corresponds to largest x^n):")

        for key, values in polynomial_coeffs.items():
            coeffs = values["coefficients"]
            r2 = values["R^2"]
            
            # Generate labels A, B, C, D, E... based on coefficient count
            labels = [chr(65 + i) for i in range(len(coeffs))]  # 65 = 'A'
            
            coeff_str = ", ".join(f"{label}={coeff:.6f}" for label, coeff in zip(labels, coeffs))
            print(f"{key}: {coeff_str}, R^2={r2:.6f}")

        print('---------------------')
        print()

    elif scaling__analysis == 2:
        print("------------------------------")
        print(" L SCALING ANALYSIS ")
        print("------------------------------")

        scaling_dict = {}
        for orbit_factor in np.arange(a, b, interval):
            matrices_dict["orbital"] = orbmoment * orbit_factor

            gmatrix, gshift = from_matrices_to_gshift(states__lengthsz, matrices_dict, ppms)

            scaling_dict.update({round(float(orbit_factor), 1): [round(float(gshift[0].real), 3), 
                            round(float(gshift[1].real), 3), 
                            round(float(gshift[2].real), 3)]})
            
            print('Scaling factor: ', round(orbit_factor, 2), ', g-shifts: ', [round(float(gshift[0].real), 3), 
                            round(float(gshift[1].real), 3), 
                            round(float(gshift[2].real), 3)])
            
        # Removing the second element from each value list
        new_scaling_dict = remove_duplicate_positions(scaling_dict)
        namee = file.split('_')[0]+"_l_scaling"
        plot_and_fit_scaling_dict(new_scaling_dict, namee, "L Factor Scaling", savepic, degree_fit)

        # Print formatted output
        polynomial_coeffs = fit_polynomial(scaling_dict, degree_fit)
        print("SOC scaling - Polynomial fit:")
        for key, values in polynomial_coeffs.items():
            coeffs = values["coefficients"]
            r2 = values["R^2"]
            print(f"{key}: A={coeffs[0]:.6f}, B={coeffs[1]:.6f}, C={coeffs[2]:.6f}, D={coeffs[3]:.6f}, R^2={r2:.6f}")

        print('---------------------')
        print()
    
    elif scaling__analysis == 3:
        print("----------------------------------------------")
        print(" SPIN SCALING ANALYSIS (Landé factor scaling)")
        print("----------------------------------------------")
        lande_factor = 2.002319304363
        a = 0
        b = lande_factor * 2 
        interval = b / 10 

        scaling_dict = {}
        for factor_lande in np.arange(a, b, interval):
            # matrices_dict["spin"] = spin_mat * spin_factor 

            gmatrix, gshift = from_matrices_to_gshift(states__lengthsz, matrices_dict, ppms, factor_lande)

            scaling_dict.update({round(float(factor_lande), 1): [round(float(gshift[0].real), 3), 
                            round(float(gshift[1].real), 3), 
                            round(float(gshift[2].real), 3)]})
            
            print('Scaling Landé factor: ', factor_lande, ', g-shifts: ', [round(float(gshift[0].real), 3), 
                            round(float(gshift[1].real), 3), 
                            round(float(gshift[2].real), 3)])
        print()
        
        # Removing the second element from each value list
        new_scaling_dict = remove_duplicate_positions(scaling_dict)
        namee = file.split('_')[0]+"_spin_scaling"
        plot_and_fit_scaling_dict(new_scaling_dict, namee, "Landé Factor value", savepic, degree_fit)

        # Print formatted output
        polynomial_coeffs = fit_polynomial(scaling_dict, degree_fit)
        print("Spin scaling - Polynomial fit:")
        for key, values in polynomial_coeffs.items():
            coeffs = values["coefficients"]
            r2 = values["R^2"]
            print(f"{key}: A={coeffs[0]:.6f}, B={coeffs[1]:.6f}, C={coeffs[2]:.6f}, D={coeffs[3]:.6f}, R^2={r2:.6f}")

        print('---------------------')
        print()

    elif scaling__analysis == 4:
        print("------------------------------")
        print(" SOC SCALING ANALYSIS")
        print("------------------------------")

        for orbit_factor in np.arange(a, b, interval):
            matrices_dict["orbital"] = orbmoment * orbit_factor
            print('Orbital factor: ', orbit_factor)

            scaling_dict = {}
            for soc_factor in np.arange(a, b, interval): 
                matrices_dict["soc"] = soc * soc_factor

                gmatrix, gshift = from_matrices_to_gshift(states__lengthsz, matrices_dict, ppms)

                scaling_dict.update({round(float(soc_factor), 1): [round(float(gshift[0].real), 3), 
                                round(float(gshift[1].real), 3), 
                                round(float(gshift[2].real), 3)]})
                
                print('Scaling factor: ', soc_factor, ', g-shifts: ', [round(float(gshift[0].real), 3), 
                                round(float(gshift[1].real), 3), 
                                round(float(gshift[2].real), 3)])
            print()
            
            # Removing the second element from each value list
            new_scaling_dict = remove_duplicate_positions(scaling_dict)
            namee = file.split('_')[0]+"_soc_scaling"+"_l_"+str(round(orbit_factor, 2))
            plot_and_fit_scaling_dict(new_scaling_dict, namee, "Scaling Factor SOC", savepic, degree_fit)

            # Print formatted output
            polynomial_coeffs = fit_polynomial(scaling_dict, degree_fit)
            print("")
            print("SOC scaling - Polynomial fit:")
            for key, values in polynomial_coeffs.items():
                coeffs = values["coefficients"]
                r2 = values["R^2"]
                print(f"{key}: A={coeffs[0]:.6f}, B={coeffs[1]:.6f}, C={coeffs[2]:.6f}, D={coeffs[3]:.6f}, R^2={r2:.6f}")

            print('---------------------')
            print()
