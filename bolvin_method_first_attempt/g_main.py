#####################################
#        MODULES SELECTION
#####################################

# module load Python/3.7.4-Anaconda3-2019.10

import numpy as np
import sys

from g_read import get_number_of_states, get_eigenenergies, get_spin_orbit_couplings
from g_read import get_spin_of_each_state, get_spin_matrices, get_orbital_matrices, get_states_same_symmetry

from g_operations import get_Hamiltonian_construction, hamiltonian_diagonalization
from g_operations import angular_matrices_obtention, g_factor_calculation

#####################################
#        INPUT VALUES
#####################################

input_file = 'ras_h2o.out'

write_file = 0  # 0: write results directly; 1: write them in output file
output_file = 'results.txt' # output name

select_state_by_symmetry = 0 # 0: use "state_ras" ; 1: use states by selected symmetry
symmetry_selection = 'B2' # Symmetry selected states
symmetry_ground_state = 1 # Ground state to be included in symmetry calculation

states_ras = [1,2,3] # States to be included when "select_state_by_symmetry = 0"

#####################################
#    HAMILTONIAN DIAGONALIZATION
#####################################

total_states = get_number_of_states(input_file)

# Take states by symmetry or, if not, verify that selected states are lower than maximum state number
if select_state_by_symmetry == 1:
    states_ras = get_states_same_symmetry(input_file, symmetry_ground_state, total_states, symmetry_selection)
else:
    for i in states_ras:
        if i <= 0 or i > total_states:
            print("The number of states selected cannot be higher than the total number of states calculated in QChem")
            sys.exit()

eigen_energies = get_eigenenergies(states_ras, input_file)

spin_orbit_coupling = get_spin_orbit_couplings(total_states, states_ras, input_file)
spin_orbit_coupling = spin_orbit_coupling / 219474.63  # From cm-1 to a.u.

Ham_ras = get_Hamiltonian_construction(states_ras, eigen_energies, spin_orbit_coupling)

la_Ham_ras, v_Ham_ras = hamiltonian_diagonalization(Ham_ras)

#############################################
#      SPIN AND LAMBDA MATRIXES FORMATION
#############################################

state_spin = get_spin_of_each_state(input_file)

spin_matrix = get_spin_matrices(input_file, states_ras, state_spin)
l_matrix = get_orbital_matrices(input_file, total_states, states_ras)

lambda_matrix = angular_matrices_obtention(la_Ham_ras, v_Ham_ras, l_matrix)
sigma_matrix = angular_matrices_obtention(la_Ham_ras, v_Ham_ras, spin_matrix)

##################################
#     G-FACTOR CALCULATION
##################################

G_tensor, G_tensor_results = g_factor_calculation( lambda_matrix , sigma_matrix )

#########################################
#         PRINTING
#########################################

if write_file == 1:
   sys.stdout = open(output_file, "w")

print("------------------------")
print("     INPUT SECTION")
print("------------------------")
print(" ")

print("File selected: ", input_file)
print("Number of states: ", total_states)
print("Selected states: ", states_ras)

print(" ")
print('Angular momentums (x,y,z):')
for k in range(0,3):
    print('Dimension: ', k)
    print('\n'.join([''.join(['{:^25}'.format(item) for item in row])
                     for row in np.round((l_matrix[:,:,k]),5)]))
    print(" ")

print("Spin angular momentums")
for k in range(0,3):
    print('Dimension: ', k)
    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])
                    for row in np.round((spin_matrix[:,:,k]), 5)]))
    print(" ")

print("------------------------")
print("   DATA OF INTEREST")
print("------------------------")
print(" ")

print("- Energy of the first selected state (a.u.) in RAS (state", states_ras[0], "): ")
print( - abs( la_Ham_ras[0] ) )
print(" ")

for i in range(1, len(states_ras) ):
    print("- Energy difference between state", states_ras[i], " and state", states_ras[0], "(cm-1):")
    print(np.round((eigen_energies[i * 2] - eigen_energies[0]) * 219474.63, 0))
print(" ")

print('SIGMA matrix with all spin angular momentums:')
print('\n'.join([''.join(['{:^25}'.format(item) for item in row])
                 for row in np.round((sigma_matrix[:,:]),8)]))
print(" ")

print('LAMBDA matrix with all spin angular momentums:')
print('\n'.join([''.join(['{:^25}'.format(item) for item in row])
                 for row in np.round((lambda_matrix[:,:]),5)]))
print("  ")

print("------------------------")
print("   g-factor results   ")
print("------------------------")
print(" ")

print('G_tensor matrix:')
print(''.join(['{:^25}'.format(item) for item in ['x','y','z']]))
print('\n'.join([''.join(['{:^25}'.format(item) for item in row])
                 for row in np.round((G_tensor[:,:]),5)]))
print(" ")

print("G-matrix eigenvalues: ")
for i in range(0, 3):
    print('dimension', i, ':', np.round(G_tensor_results.real[i], 3))
