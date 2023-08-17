__author__ = 'Sven Kähler, Antonio Cebreiro-Gallardo'

from parser_rasci import *

file = 'fe_pyms2_def2tzvp_ras.out' # str(sys.argv[1])

class MyEncoder(json.JSONEncoder):
    #https://itsourcecode.com/typeerror/typeerror-object-of-type-is-not-json-serializable-solved/?expand_article=1
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, complex):
            # return {'real': obj.real, 'image': obj.imag}
            a = obj.real
            b = obj.imag
            return {a, b}
        else:
            return super().default(obj)


def get_socs(output, totalstates, nstates):
    data = output['interstate_properties']
    socs = {}
    for i in range(2, totalstates+1):
        lista = data[(1, i)]['total_soc_mat']
        socs.update({nstates[i - 1]: lista})
    return socs


def get_energies(output, totalstates, nstates):
    data = output['excited_states']
    energy = {}
    for i in range(0, totalstates):
        lista = [data[i]['total_energy'], data[i]['excitation_energy']]
        energy.update({nstates[i]: lista})
    return energy


def get_orbital_angmoment(output, totalstates, nstates):
    data = output['interstate_properties']
    momentums = {}
    for i in range(2, totalstates+1):
        lista = data[(1, i)]['angular_momentum']
        momentums.update({nstates[i - 1]: lista})
    return momentums


with open(file, encoding="utf8") as f:
    output = f.read()
output = parser_rasci(output)

totalstates, states_ras = get_number_of_states(file)
rasAveTransSOCListDict = get_socs(output, totalstates, states_ras)
rasStateEnergiesDict = get_energies(output, totalstates, states_ras)
rasTransAngMomListDict = get_orbital_angmoment(output, totalstates, states_ras)

output_dict = {
    "rasStateEnergiesDict": rasStateEnergiesDict,
    "rasStateTotalAngMomDict": rasStateTotalAngMomDict,
    "rasTransAngMomListDict": rasTransAngMomListDict,
    "rasAveTransSOCListDict": rasAveTransSOCListDict
}

outfile_name = "prueba.json"
with open(outfile_name, 'w') as f:
    json.dump(output_dict, f,
          separators=(',', ':'),
          sort_keys=True,
          indent=4,
          cls=MyEncoder)
# https://itsourcecode.com/typeerror/typeerror-object-of-type-is-not-json-serializable-solved/?expand_article=1

exit()

eigenenergies_ras, excitation_energies_ras = get_eigenenergies(file, totalstates, states_ras)
doublet_socs, sz_values, sz_ground = get_spin_orbit_couplings(file, totalstates, states_ras)
g_shift = from_energies_soc_to_g_values(file, states_ras, totalstates, excitation_energies_ras, doublet_socs, sz_values)
print('g-factor:', np.round(g_shift.real[0], 3), np.round(g_shift.real[1], 3), np.round(g_shift.real[2], 3))
