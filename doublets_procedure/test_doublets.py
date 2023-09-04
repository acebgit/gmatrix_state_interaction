#####################################
#          MODULES SELECTION
#####################################
import numpy as np
import sys

from doublets_procedure.parser_init import *
from doublets_procedure.parser_gtensor import *
from doublets_procedure.parser_excitstates import *

#####################################
#            INPUT VALUES
#####################################
ras_input = '\
0_doublets_molecules/co2/co2_def2-TZVP_11_9.out'  # str(sys.argv[1])

selected_states = 1  # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [1, 4, 5]  # States to be included when "selected_states = 0"
symmetry_selection = 'A2'  # Symmetry selected states
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

totalstates = get_number_of_states(ras_input)

states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

doublet_socs, sz_values, sz_ground = get_spin_orbit_couplings(ras_input, totalstates, states_ras, soc_options, bolvin=1)

ras_upper_g_matrix, g_shift = bolvin_from_energies_soc_to_g_values(ras_input, states_ras,
                                                                   totalstates, excitation_energies_ras,
                                                                   doublet_socs, sz_values)

print_g_calculation(ras_input, totalstates, selected_states, symmetry_selection, states_ras, g_shift)
