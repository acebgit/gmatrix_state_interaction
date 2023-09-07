__author__ = 'Antonio Cebreiro-Gallardo'
from gtensor_doublets import *

ras_input = 'example_doublets.out'  # str(sys.argv[1])

selected_states = 0  # 0: use "state_ras" ; 1: use all states ; # 2: use states by selected symmetry "symmetry_selection"
states_ras = [1, 2, 3]  # States to be included when "selected_states = 0"
symmetry_selection = 'A1'  # Symmetry selected states
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

gfactor_obtention(ras_input, states_ras, selected_states, symmetry_selection, soc_options)
