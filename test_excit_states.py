
from parser_excitstates import get_excited_states_analysis

ras_input = '\
triplets_molecules/C64H28B4_4_4_concat.out'

get_excited_states_analysis(ras_input, cutoff=0.9)
