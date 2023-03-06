#####################################
#          MODULES SELECTION
#####################################
# module load Python/3.7.4-Anaconda3-2019.10
import numpy as np
import sys
from tabulate import tabulate

#####################################
#            INPUT VALUES
#####################################
# G-TENSOR CALCULATION
g_calculation = 1
theory_level = 0 # 0: RAS ; 1: CASSCF ; 2: SF-DFT
input = '../o2_def2-TZVP_11_9/o2_def2-TZVP_11_9.out'
#input = str(sys.argv[1])

selected_states = 1 # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [1,2] # States to be included when "selected_states = 0"
symmetry_selection = 'Ag' # Symmetry selected states
selected_SOC = 0 # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

write_file = 0 # 0: write results directly; 1: write in output file
plots = 0

#####################################
#      G-VALUE CALCULATION
#####################################
if (g_calculation == 1):

    from g_read_triplets import get_number_of_states, get_eigenenergies, get_selected_states, \
        get_spin_orbit_couplings, get_spin_matrices, get_orbital_matrices

    from g_operations_triplets import get_Hamiltonian_construction, Hamiltonian_diagonalization, \
        angular_matrixes_obtention, g_factor_calculation, plot_obtention

    totalstates = get_number_of_states(input)

    states_ras = get_selected_states(input, totalstates, symmetry_selection, states_ras, selected_states)

    eigenenergies, excitation_energies_eV = get_eigenenergies(states_ras, input, theory_level, totalstates)
    excitation_energies_au = excitation_energies_eV / 27.211399

    spin_orbit_coupling = get_spin_orbit_couplings(selected_SOC, totalstates, states_ras, input)

    Ham_ras = get_Hamiltonian_construction(states_ras, excitation_energies_au, spin_orbit_coupling)

    eigenvalues, eigenvector, kramer_state = Hamiltonian_diagonalization(Ham_ras)

    spin_matrix = get_spin_matrices(input, states_ras)

    l_matrix = get_orbital_matrices(input, totalstates, states_ras)

    sigma_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramer_state, spin_matrix)

    lambda_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramer_state, l_matrix)

    G_tensor, G_tensor_results = g_factor_calculation( lambda_matrix , sigma_matrix )

    ##############################
    #      RAS OUTPUT PRINT
    ##############################
if (g_calculation == 1):
    print(np.round(G_tensor_results.real[0], 3), np.round(G_tensor_results.real[1], 3), \
          np.round(G_tensor_results.real[2], 3))

    output_file = input + '-gvalues.txt' # output name
    if (write_file == 1):
       sys.stdout = open(output_file, "w")

    print(" ")
    print("------------------------")
    print("     INPUT SECTION")
    print("------------------------")
    print("File selected: ", input)
    print("Number of states: ", totalstates)
    if (selected_states == 2):
        print("Symmetry: ", symmetry_selection)
        print("Selected states: ", states_ras)
    else:
        print("Selected states: ", states_ras)

    print(" ")
    print("------------------------")
    print("    OUTPUT SECTION")
    print("------------------------")
    print('g-factor (x y z dimensions):')
    print(np.round(G_tensor_results.real[0], 3), np.round(G_tensor_results.real[1], 3), \
          np.round(G_tensor_results.real[2], 3))
    print('')

    print('G_tensor matrix:')
    print(''.join(['{:^15}'.format(item) for item in ['x','y','z']]) )
    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
                     for row in np.round((G_tensor[:,:]),5)]))
    print('')
    print(np.round(G_tensor_results.real[0], 3), np.round(G_tensor_results.real[1], 3), \
          np.round(G_tensor_results.real[2], 3))
    print('')

    ##############################
    #        PLOTS PRINT
    ##############################

    if (plots == 1):
        name = input + '_energies.png'
        plot_obtention(name, states_ras, excitation_energies_eV)

        soc_plotting_data = np.zeros( len(states_ras) )

        for i in range(0, len(spin_orbit_coupling), 2):
            soc_spin_1 = np.absolute( spin_orbit_coupling[i, 0] )
            soc_spin_2 = np.absolute( spin_orbit_coupling[i+1, 0] )

            soc_plotting_data[i // 2] = max(soc_spin_1, soc_spin_2) * 27.211399

        name = input + '_SOC.png'
        plot_obtention(name, states_ras, soc_plotting_data)

