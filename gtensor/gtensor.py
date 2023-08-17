__author__ = 'Sven Kähler, Antonio Cebreiro-Gallardo'

from parser_rasci import *

file = 'fe_pyms2_def2tzvp_ras.out' # str(sys.argv[1])
totalstates, states_ras = get_number_of_states(file)
eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
doublet_socs, sz_values, sz_ground = get_spin_orbit_couplings(file, totalstates, states_ras)
g_shift = from_energies_soc_to_g_values(file, states_ras, totalstates, excitation_energies_ras, doublet_socs, sz_values)
print('g-factor:', np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3), np.round(g_shift.real[2], 3))
