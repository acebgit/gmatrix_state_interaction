"""
PROGRAM TO CALCULATE THE CORRELATION
BETWEEN ENERGIES-SOC-ORBITAL ANGULAR MOMENTUM
AND THE G-TENSOR VALUES
"""
import numpy as np
from scipy import stats

from parser_gtensor import get_number_of_states, get_selected_states, get_eigenenergies, get_spin_orbit_couplings, \
    from_energies_soc_to_g_values, get_hamiltonian_construction, diagonalization, get_spin_matrices, \
    get_orbital_matrices, angular_matrices_obtention, g_factor_calculation
from parser_plots import plot_g_tensor_vs_states

ras_input = 'triplets_molecules/o2_11_9_triplets.out'  # str(sys.argv[1])
state_selection = 0  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1, 5]  # States to be included when "states_option = 0"
symmetry_selection = 'A2'  # Symmetry selected states_selected
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

soc_analysis = 1
ener_analysis = 0
orbit_analysis = 0

def soc_and_gvalues_correlation(file, n_states, totalstate, excit_ener, socs, list_sz, ground_sz):
    """"
    Correlation analysis between SOC and gvalues. Variation is done only
    for imaginary part of SOC between ground and first excited state.
    :param: qchem_file, states_msnull, states_option, excit_energ, doublet_socs,socc_values
    :return: soc_gvalues_matrix, r_square of fit
    """
    presentation_tuple = []
    if len(n_states) != 2:
        raise ValueError("SOC correlation only works with two states")

    min_socc = 0
    max_socc = 0
    for i in range(0, len(socs)):
        if socs[0, i] != 0:
            min_socc = socs[0, i] / 10
            max_socc = socs[0, i] * 3
            break
    step_socc = 20

    for soc_value in np.linspace(min_socc, max_socc, step_socc):
        socs[3, 1] = soc_value.real + soc_value.imag * 1j
        socs[4, 2] = soc_value.real + soc_value.imag * 1j
        socs[4, 0] = -soc_value.real + soc_value.imag * 1j
        socs[5, 1] = -soc_value.real + soc_value.imag * 1j
        for i in range(0, len(socs)):
            for j in range(i, len(socs)):
                socs[i, j] = np.conjugate(socs[j, i])

        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) \
        #                  for row in np.round((socs[:, :]), 5)]))
        # print("----")
        # exit()

        g_shift = from_energies_soc_to_g_values(file, n_states, totalstate,
                                                excit_ener, socs, list_sz, ground_sz)
        presentation_tuple.append([abs(soc_value), np.round(g_shift.real[0], 3),
                                   np.round(g_shift.real[1], 3), np.round(g_shift.real[2], 3)])

    presentation_matrix = np.array(presentation_tuple, dtype=object)
    presentation_matrix[:, 0] = presentation_matrix[:, 0] * 27.211399  # from au to eV
    plot_g_tensor_vs_states(file, presentation_matrix, x_title='SOC (eV)', y_title='g-values (ppt)',
                            main_title='g-factor evolution with SOC', save_options=0)

    # R-squares obtention
    r_square = []
    for i in range(1, len(presentation_matrix[0])):
        res = stats.linregress(presentation_matrix[:, 0].astype(float), presentation_matrix[:, i].astype(float))
        r_square.append(res.rvalue ** 2)

    return presentation_matrix, r_square


def energy_and_gvalues_correlation(file, n_states, excit_energies, socs, list_sz, ground_sz):
    """"
    Correlation analysis between energies and gvalues. Change is done only for energy
    between ground and first excited state.
    :param: qchem_file, states_msnull, states_option, excit_energ, doublet_socs,SOCC_values
    :return: soc_gvalues_matrix, r_square of fit
    """
    presentation_tuple = []

    max_ener = (max(excit_energies)) * 3
    min_ener = (max(excit_energies)) / 10
    steps_ener = 25

    for ener_value in np.linspace(min_ener, max_ener, steps_ener):
        excit_energies[1] = ener_value

        g_shift = from_energies_soc_to_g_values(file, n_states, totalstates,
                                                excitation_energies_ras, socs, list_sz, ground_sz)

        presentation_tuple.append([ener_value * 27.211399, np.round(g_shift.real[0], 3),
                                   np.round(g_shift.real[1], 3), np.round(g_shift.real[2], 3)])
    presentation_matrix = np.array(presentation_tuple, dtype=object)
    plot_g_tensor_vs_states(file, presentation_matrix, x_title='Energy (eV)', y_title='g-values (ppt)',
                            main_title='g-factor evolution with Excitation Energy', save_options=0)

    return presentation_matrix


def orbitmomentum_and_gvalues_correlation(file, nstates, excit_energ, socs, list_sz, ground_sz):
    """"
    Correlation analysis between energies and gvalues. Change is done only for energy
    between ground and first excited state.
    :param: qchem_file, states_msnull, states_option, excit_energ, doublet_socs,SOCC_values
    :return: soc_gvalues_matrix, r_square of fit
    """
    hamiltonian_ras = get_hamiltonian_construction(nstates, excit_energ, socs, list_sz)

    eigenvalue, eigenvector, diagonal_mat = diagonalization(hamiltonian_ras)

    spin_matrix, standard_spin_matrix = get_spin_matrices(file, nstates)

    combination_spin_matrix = angular_matrices_obtention(eigenvector, spin_matrix, list_sz)

    orbital_matrix = get_orbital_matrices(file, totalstates, nstates, list_sz)
    orbital_matrix[:, :, :] = 0

    # Change L values
    presentation_tuple = []
    max_momentum = 0.99
    min_momentum = 0.01

    for l_value in np.linspace(min_momentum, max_momentum, 30):
        l_value = 0 + l_value * 1j

        orbital_matrix[3, 0, 0] = l_value
        orbital_matrix[4, 1, 0] = l_value
        orbital_matrix[5, 2, 0] = l_value
        for i in range(0, len(orbital_matrix)):
            for j in range(i, len(orbital_matrix)):
                orbital_matrix[i, j, 0] = np.conjugate(orbital_matrix[j, i, 0])
        # print('\n'.join([''.join(['{:^15}'.format(item) for item in row]) for row in
        #                  np.round((orbital_matrix[:, :, 0]), 5)]))
        # print("----")
        # exit()

        combination_orbital_matrix = angular_matrices_obtention(eigenvector, orbital_matrix, list_sz)
        g_shift = g_factor_calculation(standard_spin_matrix, combination_spin_matrix, combination_orbital_matrix,
                                       list_sz, ground_sz)
        presentation_tuple.append([abs(l_value), np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3),
                                   np.round(g_shift.real[2], 3)])

    presentation_matrix = np.array(presentation_tuple, dtype=object)
    plot_g_tensor_vs_states(file, presentation_matrix, x_title='Orbital Angular Momentum', y_title='g-values (ppt)',
                            main_title='g-factor evolution with orbital angular momentum', save_options=0)

    # R-squares obtention
    r_square = []
    for i in range(1, len(presentation_matrix[0])):
        res = stats.linregress(presentation_matrix[:, 0].astype(float), presentation_matrix[:, i].astype(float))
        r_square.append(res.rvalue ** 2)
    return presentation_matrix, r_square


# 1st part of the calculation
totalstates = get_number_of_states(ras_input)
states_ras = get_selected_states(ras_input, totalstates, states_ras, state_selection, symmetry_selection)
eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)
selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options)

# SOC CORRELATION
if soc_analysis == 1:
    soc_gvalues_matrix, r2 = soc_and_gvalues_correlation(ras_input, states_ras, totalstates, excitation_energies_ras,
                                                         selected_socs, sz_list, sz_ground)
    print('R-square in three dimmensions: ', np.round(r2[0], 2), np.round(r2[1], 3), np.round(r2[2], 3))

# ENERGY CORRELATION
if ener_analysis == 1:
    ener_gvalues_matrix = energy_and_gvalues_correlation(ras_input, states_ras, excitation_energies_ras, selected_socs,
                                                     sz_list, sz_ground)

# ORBITAL ANGULAR MOMENTUM CORRELATION
if orbit_analysis == 1:
    momentum_gvalues_matrix, r2_momentum = orbitmomentum_and_gvalues_correlation(ras_input, states_ras,
                                                                                 excitation_energies_ras, selected_socs,
                                                                                 sz_list, sz_ground)
    print('R-square in three dimmensions: ', np.round(r2_momentum[0], 2),
      np.round(r2_momentum[1], 3), np.round(r2_momentum[2], 3))
