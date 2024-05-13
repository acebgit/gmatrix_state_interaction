import json
import sys 
import numpy as np
import pandas as pd
from scipy import constants

from projectmethod.parsers.parser_gtensor import get_hamiltonian_construction, diagonalization, \
    angular_matrices_obtention, g_factor_calculation
from projectmethod.parsers.parser_excitstates import get_bar_chart
from projectmethod.parsers.parser_plots import plot_g_tensor_vs_states

# INPUT FILE
# file = str(sys.argv[1])
file = '../../projectmethod_allparsers/test/qchem_eomcc.out'

######## G-TENSOR CALCULATION ########
g_calculation = 1
ppm = 0 # 0: ppt; 1: ppm
# state_selection = 1 # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
# states_ras = [1,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20]
# symmetry_selection = 'B1u'  # Symmetry selected states_selected
# soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

######## G-TENSOR ANALYSIS ########
excitanalysis_gvalue_cut = 0.1 # =0: not calculate; ≠0: cut-off between ground-excited states (% of maximum g-value in each dim)
# gestimation = 0 # 0: g-tensor calculation (projection procedure); 1: g-tensor estimation (g = -4 L SOC / E)

######## EXCITED STATES ANALYSIS ########
excitanalysis = 1
excitanalysis_plot = 0
# excitanalysis_config_cut = 0.5 # cut-off for configurations amplitude (% of maximum amplitude)
# excitanalysis_soc_cut = 0 # cut-off for soccs (cm-1)
# excitanalysis_angmoment_cut = 0 # cut-off for orbital angular momentum (cm-1)

######## SOS PLOTS ########
sos_analysis = 1 # SOS g-tensor plot: g-tensor calculation with n states
# gestimation_comparison = 0 # 1: SOS comparison between g-shift calculated and estimated


def extract_data_from_json(filee):
    """
    Get lists with the information from "json" output filee.
    :param filee:
    :return: total_energy, excitation_energy_list, spin_list, soc_list, orbital_momentum_list
    """
    with open(filee, 'r') as f:
        object_text = f.read()
    input_dict = json.loads(object_text)

    total_energy_list = []
    excitation_energy_list = []
    for i in input_dict['selected_energy_dict']:
        total_energy_list.append(float(input_dict['selected_energy_dict'][i][0]))
        excit_energy = float(input_dict['selected_energy_dict'][i][1]) 
        excitation_energy_list.append(excit_energy)

    spin_list = [float(input_dict['spin_dict'][i]) for i in input_dict['spin_dict']]
    soc_list = [input_dict['soclist_dict'][i] for i in input_dict['soclist_dict']]
    orbital_momentum_list = [[complex(input_dict['orbitalmomentlist_dict'][i][0]),
                            complex(input_dict['orbitalmomentlist_dict'][i][1]),
                            complex(input_dict['orbitalmomentlist_dict'][i][2])] 
                            for i in input_dict['orbitalmomentlist_dict']]
    transitions_list = [i for i in input_dict['transitions_dict']]
    return total_energy_list, excitation_energy_list, spin_list, soc_list, orbital_momentum_list, transitions_list


def s2_to_s(s2):
    """
    get total spin (s) from s^2
    :param: s2_all
    :return: total spin (s)
    """
    return 0.5 * (-1 + np.sqrt(1 + 4 * s2))


def get_input_data(spin_state):
    """
    Get: i) number of states, ii) a list from -sz to sz, where sz = maximum multiplicity,
    iii) a list from -sz to sz, where sz = ground state spin multiplicity
    :param spin_state:
    :return:
    """
    max_multip = float(s2_to_s(max(spin_state)))
    max_multip_szlist = list(np.arange(-max_multip, max_multip + 1, 1))
    s_ground = float(s2_to_s(spin_state[0]))
    if s_ground == 0:
        raise ValueError("Warning! It is not allowed the calculation of the g-tensor in a singlet ground state. "
                         "Ground state corresponds to the first of the included states selected.")
    ground_state_szlist = list(np.arange(-s_ground, s_ground + 1))

    return max_multip_szlist, ground_state_szlist


def from_soclist_socmatrix(soc_list, maxsz_list):
    """
    Construct the SOC matrix from the SOC list of json filee.
    :param soc_list:
    :param maxsz_list:
    :return: socmatrix
    """
    len_sz = len(maxsz_list)  # dimension determined by the maximum multiplicity
    nstate = len(soc_list) + 1  # +1 is the ground state
    socmatrix = np.zeros((nstate * len_sz, nstate * len_sz), dtype=complex)

    for i in range(1, nstate):  # ground state does not have soc, and there is no state j
        soc_list_sz_state1 = len(soc_list[i-1][0])
        soc_list_sz_state2 = len(soc_list[i-1])

        for sz_1 in range(0, soc_list_sz_state1):
            for sz_2 in range(0, soc_list_sz_state2):
                matrix_row = sz_1 + ((len_sz-soc_list_sz_state1)//2)
                matrix_col = i * len_sz + sz_2 + ((len_sz-soc_list_sz_state2)//2)

                # print('State', i, 'Matrix:', matrix_row, matrix_col, '--> List:', i-1, sz_2, sz_1)
                value = complex(soc_list[i-1][sz_2][sz_1])
                socmatrix[matrix_row][matrix_col] = np.conj(value)
                socmatrix[matrix_col][matrix_row] = value

    # print('SOC:')
    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #                 for row in np.round((socmatrix[:,:]))]))
    # exit()
    cm_to_ev = constants.physical_constants['inverse meter-electron volt relationship'][0] * 100
    socmatrix = socmatrix * cm_to_ev
    return socmatrix


def get_spin_matrices(states_spin, maxsz_list):
    """
    Get spin matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the order of
    "selected_state". Functions "s2_to_s", "s_to_s2" and "spin_matrices" are taken from inside PyQChem.
    :param: filee, selected_state.
    :return: spin_matr, standard_spin_mat.
    """
    def s_to_s2(spin):
        """
        get s2_all from total spin (s)
        :param: spin: total spin (s)
        :return: s2_all
        """
        return spin * (spin + 1)

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

    def expand_spin_matrices(s_x, s_y, s_z, max_mult, state_mult):
        """
        Expand (sxx, syy, szz) matrices with dimension of the st multiplicity to (sxx, syy, szz) with dimension of the
        maximum multiplicity of states.
        :param: s_x, s_y, s_z, max_sz_listt, state_mult
        :return: long_sx, long_sy, long_sz
        """
        long_sx = np.zeros((max_mult, max_mult), dtype=complex)
        long_sy = np.zeros((max_mult, max_mult), dtype=complex)
        long_sz = np.zeros((max_mult, max_mult), dtype=complex)

        multipl_difference = int((max_mult - state_mult) // 2)

        if multipl_difference != 0:
            for n_row in range(0, len(sx)):
                for n_column in range(0, len(sx)):
                    iii = n_row + multipl_difference
                    jj = n_column + multipl_difference
                    long_sx[iii, jj] = s_x[n_row, n_column]
                    long_sy[iii, jj] = s_y[n_row, n_column]
                    long_sz[iii, jj] = s_z[n_row, n_column]
        else:
            long_sx = s_x
            long_sy = s_y
            long_sz = s_z
        return long_sx, long_sy, long_sz

    def form_big_spin_matrix(statee, max_mult, sxx, syy, szz, spin_matr):
        """
        Get big spin matrix with dimensions "(len(selected_state) * max_sz_listt, len(selected_state) *
        max_sz_listt, 3)" with s_x, s_y, s_z.
        :param: st, selected_state, max_sz_listt, sxx, syy, szz, spin_matr
        :return: spin_matr
        """
        s_dim = 0
        initial_pos = statee * max_mult
        for row in range(0, max_mult):
            for column in range(0, max_mult):
                spin_matr[initial_pos, initial_pos, 0] = sxx[s_dim, s_dim]
                spin_matr[initial_pos, initial_pos + column, 0] = sxx[s_dim, s_dim + column]
                spin_matr[initial_pos + row, initial_pos, 0] = sxx[s_dim + row, s_dim]
                spin_matr[initial_pos + row, initial_pos + column, 0] = sxx[s_dim + row, s_dim + column]

                spin_matr[initial_pos, initial_pos, 1] = syy[s_dim, s_dim]
                spin_matr[initial_pos, initial_pos + column, 1] = syy[s_dim, s_dim + column]
                spin_matr[initial_pos + row, initial_pos, 1] = syy[s_dim + row, s_dim]
                spin_matr[initial_pos + row, initial_pos + column, 1] = syy[s_dim + row, s_dim + column]

                spin_matr[initial_pos, initial_pos, 2] = szz[s_dim, s_dim]
                spin_matr[initial_pos, initial_pos + column, 2] = szz[s_dim, s_dim + column]
                spin_matr[initial_pos + row, initial_pos, 2] = szz[s_dim + row, s_dim]
                spin_matr[initial_pos + row, initial_pos + column, 2] = szz[s_dim + row, s_dim + column]
        return spin_matr

    def get_standard_spin_matrix(spin_states, max_sz_listt, spin_mat):
        """
        Construct Standard Spin matrix from the previous Spin Matrix.
        :param: spin_states, max_sz_listt, spin_mat
        :return: standard_spin_mat
        """
        ground_multiplicity = int(2 * s2_to_s(spin_states[0]) + 1)
        standard_spin_mat = np.zeros((ground_multiplicity, ground_multiplicity, 3), dtype=complex)

        multip_difference = (max_sz_listt - ground_multiplicity) // 2
        for k in range(0, 3):
            for ii in range(0, ground_multiplicity):
                for jj in range(0, ground_multiplicity):
                    standard_spin_mat[ii, jj, k] = spin_mat[ii + multip_difference, jj + multip_difference, k]
        return standard_spin_mat

    len_sz = len(maxsz_list)  # dimension determined by the maximum multiplicity
    nstate = len(states_spin)
    spinmatrix = np.zeros((nstate * len_sz, nstate * len_sz, 3), dtype=complex)

    for state in range(0, nstate):
        # Form the spin matrix of the state
        s2_state = states_spin[state]
        s = s2_to_s(s2_state)
        sx, sy, sz = spin_matrices(s)

        # Expand the spin matrix of the st to the dimension of the maximum multiplicity (max_sz_listt)
        state_multip = int(2 * s + 1)
        sx, sy, sz = expand_spin_matrices(sx, sy, sz, len_sz, state_multip)

        # Mix (sxx,syy,szz) in one spin matrix:
        spinmatrix = form_big_spin_matrix(state, len_sz, sx, sy, sz, spinmatrix)

    # Take Standard Spin Matrix from Spin Matrix
    standardspin_matrix = get_standard_spin_matrix(states_spin, len_sz, spinmatrix)

    # print('\n'.join([''.join(['{:^8}'.format(item) for item in row])\
    #                 for row in np.round((standardspin_matrix[:,:,0]))]))
    # exit()
    return spinmatrix, standardspin_matrix


def get_orbital_matrices(states_momentum, maxsz_list):
    """
    Get orbitals angular momentum matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the
    order of "selected_state".
    :param: filee, totalstates, selected_states, sz_list
    :return: all_multip_lk
    """
    def get_selected_states_momentum(momentum_states):
        """
        Get Lk between the selected states selected in x,y,z dimensions.
        :param: selected_states, all_momentum
        :return: selected_momentum
        """
        n_states = len(momentum_states) + 1  # +1 because of the ground state
        selected_momentum = np.zeros((n_states, n_states, 3), dtype=complex)

        for i in range(1, n_states):
            for ndim in range(0, 3):
                selected_momentum[i][0][ndim] = momentum_states[i-1][ndim]
                selected_momentum[0][i][ndim] = np.conj(momentum_states[i - 1][ndim])
        return selected_momentum

    def get_all_multip_momentum(selected_momentums, max_szlist):
        """
        Get Lk between the selected states selected in x,y,z dimensions for doublets,
        i.e. in each row (column) there are state A,-1/2 and A,+1/2.
        :param: selected_momentums, max_szlist
        :return: big_orbit_matrix
        """
        n_states = len(selected_momentums)
        multip = len(max_szlist)
        big_orbit_matrix = np.zeros((n_states * multip, n_states * multip, 3), dtype=complex)

        for k in range(0, 3):
            for i in range(0, len(selected_momentums) * len(max_szlist), len(max_szlist)):
                for j in range(0, len(selected_momentums) * len(max_szlist), len(max_szlist)):

                    for multip in range(0, len(max_szlist)):
                        big_orbit_matrix[i + multip, j + multip][k] = \
                            selected_momentums[i // len(max_szlist)][j // len(max_szlist)][k]
        return big_orbit_matrix

    selected_lk = get_selected_states_momentum(states_momentum)
    all_multip_lk = get_all_multip_momentum(selected_lk, maxsz_list)

    # for ndim in range(0, 3):
    #     print()
    #     print('\n'.join([''.join(['{:^8}'.format(item) for item in row]) \
    #                      for row in np.round(all_multip_lk[:, :, ndim], 3)]))
    # exit()
    return all_multip_lk


def print_g_calculation(filee, totalstates, upper_g_tensor_results_ras):
    print("---------------------")
    print(" G-SHIFT RESULTS")
    print("---------------------")
    print("File selected: ", filee)
    print("Number of states: ", totalstates)
    print('g-factor (x y z dimensions):')
    print(np.round(upper_g_tensor_results_ras[0].real, 3), np.round(upper_g_tensor_results_ras[1].real, 3),
          np.round(upper_g_tensor_results_ras[2].real, 3))
    print('')


def from_json_to_matrices(json_energies, json_excitenergies, json_spin, json_soc, json_orbitmoment):
    max_sz_list, sz_ground = get_input_data(json_spin)

    soc_matrix = from_soclist_socmatrix(json_soc, max_sz_list)

    spin_matrix, standard_spin_matrix = get_spin_matrices(json_spin, max_sz_list)

    orbital_matrix = get_orbital_matrices(json_orbitmoment, max_sz_list)
    return max_sz_list, sz_ground, soc_matrix, spin_matrix, standard_spin_matrix, orbital_matrix


def from_matrices_to_gshift(excitenergies_json, soc, max_sz, spin, orbital, standard_spin, szground, ppms):
    hamiltonian = get_hamiltonian_construction(excitenergies_json, soc, max_sz)

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian)

    combination_spin_matrix = angular_matrices_obtention(eigenvector, spin, max_sz)

    combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital, max_sz)

    g_shift = g_factor_calculation(standard_spin, combination_spin_matrix, combination_orbital_matrix,
                                max_sz, szground, ppms)
    return g_shift


def excitedstates_analysis(nstate, excitenergies, orbitmoment, soc, plot):
    """
    Get a table in pandas with information for each of the excited states. 
    """
    max_orbitmoment = [abs(max(sublista, key=abs)) for sublista in orbitmoment]
    socc = [((sum((abs(complex(elemento)))**2 for sublista in i for elemento in sublista))**(1/2)) for i in soc]

    presentation_list = []
    list_states = []
    for state in range(0, nstate):
        for transition in list(transitions_json[state]):
            if state == 0:
                presentation_list.append([excitenergies[state], "---", "---", 
                                        transition['transition/SOMO'], 
                                        transition['amplitude']])
            else:
                presentation_list.append([np.round(excitenergies[state],3), 
                                        np.round(max_orbitmoment[state-1],3), 
                                        np.round(socc[state-1],3),  
                                        transition['transition/SOMO'], 
                                        transition['amplitude']])
            list_states.append(state+1)
    
    # Set display options to show all rows and columns
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame(presentation_list, index=list_states, columns=['energy', 'orbit mom', 'socc', 
                                                                    'transition/SOMO', 'amplitude'])
    
    print("--------------------------")
    print(" EXCITED STATES ANALYSIS")
    print("--------------------------")
    print(df)
    print()

    if plot == 1: 
        get_bar_chart(file[:-4], [i for i in range(1,nstate+1)], excitenergies, 'Electronic State',
                        'Excitation energy (eV)', 'energ_analysis', save_pict=0)
        get_bar_chart(file[:-4], [i for i in range(2,nstate+1)], max_orbitmoment, 'Electronic State',
                        'Orbital angular momentum', 'orbit_analysis', save_pict=0)
        get_bar_chart(file[:-4], [i for i in range(2,nstate+1)], socc, 'Electronic State',
                        'Spin-orbit coupling constants (cm-1)', 'soc_analysis', save_pict=0)


def sum_over_state_plot(energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json, ppm):
        """
        Generate the sum-over-states plot, i.e. calculation of the g-tensor by including states
        from 1 to nstates. 
        :param: 
        :return: shows SOS plot
        """
        presentation_list = []
        for state in range(1, len(energies_json)+1):
            # State properties
            energies = energies_json[0:state+1]
            excitenergies = excitenergies_json[0:state+1]
            spin = spin_json[0:state+1]

            # Interstate properties
            soc = soc_json[0:state]
            orbitmoment = orbitmoment_json[0:state]

            max_sz_list, sz_ground, soc_matrix, spin_matrix, standard_spin_matrix, orbital_matrix \
                = from_json_to_matrices(energies, excitenergies, spin, soc, orbitmoment)

            g_shift = from_matrices_to_gshift(excitenergies, soc_matrix, 
                                            max_sz_list, spin_matrix, orbital_matrix, 
                                            standard_spin_matrix, sz_ground, ppm)
            
            presentation_list.append([state, np.round(g_shift[0].real, 3),
                                            np.round(g_shift[1].real, 3), np.round(g_shift[2].real, 3)])
        presentation_matrix = np.array(presentation_list, dtype=object)

        print("------------------------------")
        print(" SUM-OVER-STATE ANALYSIS")
        print("------------------------------")
        df = pd.DataFrame([row[1:4] for row in presentation_list], list(range(1, nstates+1)), columns=['gxx','gyy','gzz'])
        print(df)
        print()

        y_title = r'$\Delta g, ppt$'
        if ppm == 1:
            y_title = r'$\Delta g, ppm$'
        plot_g_tensor_vs_states(file[:-4],presentation_matrix,'Electronic State',y_title,'sos_analysis', save_options=0)


def get_gtensor_analysis(energies_json, excitenergies_json, spin_json, soc_json, 
                         orbitmoment_json, transitions_json, ppm, cut_gvalue):
    """
    Obtaining a matrix with several data for each excited state. The cut-off determines the fraction of the amplitude
    of the 1st configuration that need to have the other configurations to be shown in each state.
    :param: file_ms_notnull, cutoff
    :return: excit_matrix
    """
    # Calculate the g-tensor for each ground-excited state pair and save it in g_shift_list
    g_shift_list = []
    for state in range(1, len(energies_json)):
        energies = [energies_json[0], energies_json[state]]
        excitenergies = [excitenergies_json[0], excitenergies_json[state]]
        spin = [spin_json[0], spin_json[state]]
        soc = [soc_json[state-1]]
        orbitmoment = [orbitmoment_json[state-1]]

        max_sz_list, sz_ground, soc_matrix, spin_matrix, standard_spin_matrix, orbital_matrix \
            = from_json_to_matrices(energies, excitenergies, spin, soc, orbitmoment)

        g_shift = from_matrices_to_gshift(excitenergies, soc_matrix, 
                                        max_sz_list, spin_matrix, orbital_matrix, 
                                        standard_spin_matrix, sz_ground, ppm)
        
        g_shift_list.append([np.round(g_shift[0].real, 3),np.round(g_shift[1].real, 3), np.round(g_shift[2].real, 3)])

    # Forming different list depending on if g-shift is calculated or not
    cut_gxx = cut_gvalue * max([abs(row[0]) for row in g_shift_list])
    cut_gyy = cut_gvalue * max([abs(row[1]) for row in g_shift_list])
    cut_gzz = cut_gvalue * max([abs(row[2]) for row in g_shift_list])

    # For the list with the data for all the configurations
    presentation_list = []
    list_states = []
    for state in range(0, len(energies_json)):
        for transition in list(transitions_json[state]):
            if state == 0:
                presentation_list.append([str(transition['transition/SOMO']),
                                          (transition['amplitude']),
                                          '---','---','---'])
                list_states.append(state+1)

            else:
                gxx = np.round(g_shift_list[state-1][0],3)
                gyy = np.round(g_shift_list[state-1][1],3)
                gzz = np.round(g_shift_list[state-1][2],3)
                if abs(gxx) >= cut_gxx or abs(gyy) >= cut_gyy or abs(gzz) >= cut_gzz:
                    presentation_list.append([str(transition['transition/SOMO']),
                                              (transition['amplitude']),
                                              gxx,gyy,gzz])
                    list_states.append(state+1)
    
    df = pd.DataFrame(presentation_list, index=list_states, columns=['transition/SOMO', 'amplitude', 
                                                                     'gxx','gyy','gzz'])
    
    print("--------------------------")
    print(" G-SHIFT IN STATES PAIRS")
    print("--------------------------")
    print(df)
    print()

file = file + ".json"
energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json, transitions_json = extract_data_from_json(file)
nstates = len(energies_json)

if g_calculation == 1:
    max_sz_list, sz_ground, soc_matrix, spin_matrix, standard_spin_matrix, orbital_matrix \
    = from_json_to_matrices(energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json)
    
    g_shift = from_matrices_to_gshift(excitenergies_json, soc_matrix, max_sz_list, spin_matrix, 
                                      orbital_matrix, standard_spin_matrix, sz_ground, ppm)
    
    print_g_calculation(file, nstates, g_shift)

if excitanalysis_gvalue_cut != 0:
    get_gtensor_analysis(energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json, 
                     transitions_json, ppm, excitanalysis_gvalue_cut)

if excitanalysis == 1: 
    excitedstates_analysis(nstates, excitenergies_json, orbitmoment_json, soc_json, excitanalysis_plot)

if sos_analysis == 1:
    sum_over_state_plot(energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json, ppm)
