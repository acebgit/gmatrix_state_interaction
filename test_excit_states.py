
from parser_excitstates import get_excited_states_analysis

ras_input = '\
doublets_molecules/mncn5no_2-/mn_prueba.out'

get_excited_states_analysis(ras_input, cutoff=0.9)
