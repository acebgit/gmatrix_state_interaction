#####################################
#          MODULES SELECTION
#####################################
import numpy as np
import sys

from projection_method.parsers_OPO.parser_gtensor import get_number_of_states, get_eigenenergies, get_selected_states, \
    get_socc_values, get_ground_state_orbital_momentum, get_symmetry_states, \
    get_spin_orbit_couplings
from projection_method.parsers_OPO.parser_gtensor import bolvin_from_energies_soc_to_g_values, print_g_calculation
from projection_method.parsers_OPO.parser_excitstates import get_excited_states_analysis, improved_active_space
from projection_method.parsers_OPO.parser_plots import get_bar_chart, bolvin_sos_analysis_and_plot

#####################################
#            INPUT VALUES
#####################################
# G-TENSOR CALCULATION
g_calculation = 1
ras_input = '\
doublets_molecules/h2o/h2o_def2tzvp_5_5.out'  # str(sys.argv[1])

selected_states = 1  # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry
states_ras = [1, 4, 5]  # States to be included when "states_option = 0"
symmetry_selection = 'A2'  # Symmetry selected states_selected
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

# EXCITED STATES ANALYSIS IN ras
excited_states_analysis = 0
new_active_space = 0
sos_analysis = 0
bar_plots = 0

# eom ANALYSIS AND ras-eom ENERGIES EXCHANGE
eom_information = 0
eom_input = '../eom_outputs/example_ras.out'

eom_change_energies = 0
ras_states_to_change = [2, 3, 4]
eom_states_to_change = [4, 7, 9]

# OUTPUT
write_file = 0  # 0: write results directly; 1: write in output file_ms_notnull
output_file = ras_input + '-gvalues.txt'
if write_file == 1:
    sys.stdout = open(output_file, "w")

#####################################
#      G-VALUE CALCULATION
#####################################
if g_calculation == 1:

    totalstates = get_number_of_states(ras_input)

    states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

    doublet_socs, sz_values, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options, bolvin=1)

    ras_upper_g_matrix, g_shift = bolvin_from_energies_soc_to_g_values(ras_input, states_ras,
                                                                       totalstates, excitation_energies_ras,
                                                                       doublet_socs, sz_values)

    print_g_calculation(ras_input, totalstates, selected_states, symmetry_selection, states_ras, g_shift)

#####################################
#        EXCITED STATE ANALYSIS
#####################################
if excited_states_analysis == 1:
    get_excited_states_analysis(ras_input)

if new_active_space == 1:
    improved_active_space(ras_input)

#####################################
#        PLOT ANALYSIS
#####################################

if sos_analysis == 1:
    bolvin_sos_analysis_and_plot(ras_input)

if bar_plots == 1:

    totalstates = get_number_of_states(ras_input)

    selected_states = 1
    states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection='None')

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

    doublet_socs, sz_values, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options, bolvin=1)

    ras_upper_g_matrix, g_shift = bolvin_from_energies_soc_to_g_values(ras_input, states_ras,
                                                                       totalstates,
                                                                       excitation_energies_ras,
                                                                       doublet_socs,
                                                                       sz_values)

    # Printing excitation energies versus orbitals symmetries:
    state_symmetries, ordered_state_symmetries = get_symmetry_states(ras_input, totalstates)

    get_bar_chart(ras_input, ordered_state_symmetries, excitation_energies_ras * 27.211399, 'State', 'Energy (eV)',
                  'Excitation energies')

    # Printing socc versus orbitals symmetries:
    socc_values = get_socc_values(ras_input, totalstates)
    get_bar_chart(ras_input, ordered_state_symmetries, socc_values / 8065.540107, 'State', 'Energy (eV)',
                  'Mean-field spin-orbit coupling constants')

    # Printing orbitals angular momentum versus orbitals symmetries:
    orbital_momentum = get_ground_state_orbital_momentum(ras_input, totalstates)
    get_bar_chart(ras_input, ordered_state_symmetries, orbital_momentum, 'State', 'Orbital angular momentum', '')
