import numpy as np
import pandas as pd
from generalization_projection_technique.parsers.parser_eom import get_eom_type, get_scf_energy, get_energies, get_maximum_amplitude_orbitals

file = '../../\
generalization_projection_technique/test/qchem_eomcc.out'
transition_cutoff = 0.8
ref_state = '1/A1'

eom_version = get_eom_type(file)

scf_reference_energy = get_scf_energy(file)

states, total_energies = get_energies(file, eom_version)

# allstates_excit_energies, allstates_momentums, allstates_spins, allstates_socs = get_interstate_properties(file)

# eom_socc = get_eom_socc_values(eom_input)

orbitals_in_transitions = get_maximum_amplitude_orbitals(file, eom_version, transition_cutoff)

presentation_list = []
row_list = []
for transition in list(total_energies.keys()):
    for trans in orbitals_in_transitions:
        if trans[0]['transition'] == transition:
            ener_excit = np.round(total_energies[transition][1] - total_energies[ref_state][1], 4)
            for subtrans in trans:
                presentation_list.append([total_energies[transition][0], ener_excit, subtrans['orbitals']])
                row_list.append(transition)
df = pd.DataFrame(np.array(presentation_list), index=row_list, columns=['Total energ', 'Excit. ener.', 'transition'])

print("")
print("------------------------")
print("    eom-CC ANALYSIS ")
print("------------------------")
print("SCF   reference energy: ", scf_reference_energy)
print()
print(df.to_string())
