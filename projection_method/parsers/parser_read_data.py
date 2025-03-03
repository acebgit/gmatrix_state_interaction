import json
import sys 
from PyQchem.pyqchem.parsers.parser_rasci import parser_rasci
from PyQchem.pyqchem.parsers.parser_cis import basic_cis

def output_json(filename, output_dict):
    outfile_name = filename + ".json"
    with open(outfile_name, 'w') as archivo:
        json.dump(output_dict, archivo, separators=(',', ':'), sort_keys=False, indent=4)


def check_dimmensions(dictionary, lista, dimmension):
    for i in lista:
        if dictionary[i]:
            if len(dictionary[i]) != dimmension:
                raise ValueError("The dimensions of the dictionary {} are incorrect.".format(i))


def gtensor_parser_rasci(outpuut):
    """
    Print all the states information in RASCI: 
    "energy_dict": total energy (a.u.) and excitation energy (eV),
    "soc_matrix_dict": SOC matrix (cm-1),
    "spin_dict": states spin,
    "angmoment_dict": orbital angular momentum in the three directions,
    "transitions_dict": transitions information: state, configuration, SOMO, amplitude
    """
    def add_count_to_elements(lst):
        count_dict = {}  # To store the count of each element
        new_list = []    # New list to store the result

        for element in lst:
            # Increment the count of the element or set to 1 if it doesn't exist
            count_dict[element] = count_dict.get(element, 0) + 1
            
            # Add the count and the element to the new list
            new_list.append(f"{count_dict[element]}{element}")
        return new_list

    def get_ras_configurations(dataa, outpuut, totalstatess, cutoff_amp=0):
        """
        Obtain the orbitals with the unpaired electrons in each relevant configuration of each state.
        :return: 
        """
        def from_ras_to_scf_order(beta_homo, ras_elec_orbs, initial_orbitals):
            """
            Change from RAS order to SCF order.
            Example: if in (1,2,3,4) orbital_list the active space selected is [1,3],
            RAS_order=(1,2,3,4) while SCF_order=(2,1,3,4), since RAS is ordered by RAS1-RAS2-RAS3.
            :param: homo_orbit, initial_scf_space, ras_orbitals
            :return: final_SCF_space
            """
            ras_scf_map = {}

            for scf_orbital in range(1, beta_homo + 1): # RAS1: from orbital 1 to last beta HOMO
                if scf_orbital not in ras_elec_orbs:
                    ras_scf_map.update({len(ras_scf_map)+1: scf_orbital})

            for scf_orbital in ras_elec_orbs: # RAS2: from last beta HOMO to last alpha HOMO
                ras_scf_map.update({len(ras_scf_map)+1: scf_orbital})

            for scf_orbital in range(beta_homo + 1, beta_homo + 100): # RAS3: all empty orbitals 
                if scf_orbital not in ras_elec_orbs:
                    ras_scf_map.update({len(ras_scf_map)+1: scf_orbital})

            final_orbitals = [ras_scf_map[i] for i in initial_orbitals]
            return ','.join(map(str, final_orbitals))

        # Take RAS_ACT_ORB 
        try:
            enum = outpuut.find('RAS_ACT_ORB')
            ras_elec_orb = [int(x) for x in (outpuut[enum:enum+50].split('[')[1].split(']')[0]).split(",")]
        except: # RAS_ACT_ORB automatically selected
            print("Make program for RAS_ACT_ORB automatically selected")
            exit()

        # Form a dictionary with states, configurations, transitions and amplitudes
        transitions = []
        for state in range(totalstatess):
            state_transitions = []
            for configuration in dataa['excited_states'][state]['configurations']:

                # Take alpha and beta electrons
                # enum = outpuut.find('There are')
                # elec_alpha, elec_beta = int(outpuut[enum:enum+70].split()[2]), int(outpuut[enum:enum+70].split()[5])
                
                # Take SOMOs
                amplitude = float(configuration['amplitude'])
                if abs(amplitude) >= cutoff_amp:

                    # Alpha and beta occupations of the state
                    config_alpha = configuration['occupations']['alpha']
                    config_beta = configuration['occupations']['beta']
                    elec_beta = config_beta.count(1) # Take number of beta electrons             

                    # Take orbitals with unpaired electrons
                    somo_orbitals_ras = [i+1 for i in range(len(config_alpha)) if config_alpha[i] != config_beta[i]]

                    # From RAS order to SCF order
                    somo_orbitals_scf = from_ras_to_scf_order(elec_beta, ras_elec_orb, somo_orbitals_ras)
                    
                    # if state == 3:
                    #     print(somo_orbitals_ras, ':', somo_orbitals_scf)
                    #     exit()
                    config_index = excitstates_dict[state]['configurations'].index(configuration)

                    state_transitions.append({'state': state+1,
                                        'configuration': config_index,
                                        'SOMO': somo_orbitals_scf,
                                        'amplitude': round(amplitude, 3)})
            transitions.append(state_transitions)
        return transitions

    # Take all data with "parser_qchem"
    data = parser_rasci(outpuut)

    # STATE PROPERTIES

    # Number of states
    excitstates_dict = data['excited_states']
    totalstates = len(excitstates_dict)

    # Symmetries of the states with their order number 
    states_symmetry_nonumbered = [excitstates_dict[i]['symmetry'] for i in range(totalstates)]
    states_symmetry = add_count_to_elements(states_symmetry_nonumbered)

    # Energy (in a.u.) and excitation energy (in eV)
    energy_dict = {states_symmetry[i]: [excitstates_dict[i]['total_energy'], excitstates_dict[i]['excitation_energy']] 
                   for i in range(totalstates)}

    # Spin multiplicity 
    spin_dict = {states_symmetry[i]: excitstates_dict[i]['multiplicity'] for i in range(totalstates)}

    # ORBITALS INVOLVED IN TRANSITIONS
    transitions_dict = get_ras_configurations(data, outpuut, totalstates)

    # INTERSTATE PROPERTIES
    soc_text = 'total_soc_mat' # total_soc_mat, 1e_soc_mat, 2e_soc_mat
    interstates_dict = data['interstate_properties']
    
    soc_matrix_dict = {}
    socc_dict = {}
    angmoment_dict = {}

    for i in range(1, totalstates+1):
        for j in range(1, totalstates+1):
            bra_sym = states_symmetry[i-1]
            ket_sym = states_symmetry[j-1]
            braket = bra_sym + "_" + ket_sym
            ketbra = ket_sym + "_" + bra_sym
            
            # If braket is written in code, then ketbra is skipped 
            if i != j and braket not in soc_matrix_dict and ketbra not in soc_matrix_dict: 
                # SOC matrix: Pass numbers from complex to strings 
                pair_socs = [[str(element) for element in sublist] for sublist in interstates_dict[(i, j)][soc_text]]
                soc_matrix_dict.update({braket: pair_socs})

                # SOCC elements
                socc_dict.update({braket: interstates_dict[(i, j)]['mf_socc']})

                # Orbital angular momentum L
                pair_angmoment = [str(element) for element in interstates_dict[(i, j)]['angular_momentum']]
                angmoment_dict.update({braket: pair_angmoment})

    output_dict = {
    "scf_alpha_beta_electrons": data["scf_electrons"], 
    "energy_dict": energy_dict,
    "soc_matrix_dict": soc_matrix_dict,
    "socc_dict": socc_dict, 
    "spin_dict": spin_dict,
    "angmoment_dict": angmoment_dict,
    "transitions_dict": transitions_dict
    }

    check_dimmensions(output_dict, ["energy_dict", "spin_dict"], totalstates)
    check_dimmensions(output_dict, ["soc_matrix_dict", "socc_dict", "angmoment_dict"], (totalstates * (totalstates - 1) // 2))
    check_dimmensions(output_dict, ["transitions_dict"], totalstates)
    return output_dict


def gtensor_parser_tddft(outpuut):
    """
    Print all the states information in TDDFT. 
    Ground state is state 1 instead of state 0. 
    "energy_dict": total energy (a.u.) and excitation energy (eV),
    "soc_matrix_dict": SOC matrix (cm-1),
    "spin_dict": states spin,
    "angmoment_dict": orbital angular momentum in the three directions,
    "transitions_dict": transitions information: state, configuration, SOMO, amplitude
    """
    def s_to_s2(s):
        """
        get s from s^2
        :param: s
        :return: s^2
        """
        return ((2*s + 1)**2 - 1) / 4

    def get_tddft_configurations(dataa, totalstatess, cutoff_amp=0):
        """
        Form a dictionary with states, configurations, transitions and amplitudes
        :param: dataa, outpuut, totalstatess, cutoff_amp
        :return: transitions
        """

        scf_alpha, scf_beta = dataa['scf_electrons'][0], dataa['scf_electrons'][1]
        scf_somos = ",".join(map(str, range(min(scf_alpha, scf_beta) + 1, max(scf_alpha, scf_beta) + 1)))
        transitions = [[{'state': 1,
                        'configuration': 1, 
                        'SOMO': scf_somos, 
                        'amplitude': 1.0}]]
        
        for state in range(totalstatess):
            state_transitions = []
            for configuration in dataa['excited_states'][state]['configurations']:

                # Take SOMOs of relevant configurations
                amplitude = float(configuration['amplitude'])

                if abs(amplitude) >= cutoff_amp:
                    # Take alpha and beta occupations
                    config_alpha = configuration['occupations']['alpha']
                    config_beta = configuration['occupations']['beta']

                    # Take SOMOs
                    somo_orbitals_list = [i+1 for i in range(len(config_alpha)) if config_alpha[i] != config_beta[i]]                    
                    somo_orbitals = ','.join(map(str, somo_orbitals_list))
                    config_index = excitstates_dict[state]['configurations'].index(configuration)
                    
                    state_transitions.append({'state': state+2,
                                        'configuration': config_index+1,
                                        'SOMO': somo_orbitals,
                                        'amplitude': round(amplitude, 3)})
            transitions.append(state_transitions)
        return transitions

    # Take all data with "parser_qchem"
    data = basic_cis(outpuut)

    # STATE PROPERTIES

    # Number of states 
    excitstates_dict = data['excited_states']
    totalstates = len(excitstates_dict)

    # Energy (in a.u.) and excitation energy (in eV) 
    energy_dict = {1: [data['scf_energy'], 0.000], # Zero state energy (SCF)
                   **{i+2: [excitstates_dict[i]['total_energy'], excitstates_dict[i]['excitation_energy']]
                   for i in range(totalstates)}}

    # Spin multiplicity 

    try: # In case of doublets, spin is a float
        spin_dict = {1: data['scf_multiplicity'],  # Zero state energy (SCF)
                    **{i+2: float(excitstates_dict[i]['multiplicity']) for i in range(totalstates)}}

    except: # In case of singlets, spin is a string
        mapping_spin = {'Singlet': s_to_s2(0), 'Triplet': s_to_s2(1), 'Quintet': s_to_s2(2), 'Heptet': s_to_s2(3)}
        spin_dict = {1: data['scf_multiplicity'],  # Zero state energy (SCF)
                 **{i+2: mapping_spin[excitstates_dict[i]['multiplicity']] for i in range(totalstates)}}
    
    # ORBITALS INVOLVED IN TRANSITIONS
    transitions_dict = get_tddft_configurations(data, totalstates)

    # INTERSTATE PROPERTIES
    soc_matrix_dict = {}
    socc_dict = {}
    angmoment_dict = {}

    soc_text = 'mf_soc_mat' # 1e_socc, mf_soc_mat
    socc_text = 'mf_socc' # 1e_socc, mf_socc

    for i in range(0, totalstates+1): # Totalstates + the ground state 
        for j in range(i+1, totalstates+1):
            
            try: # If there is SOC 
                complex_list = [[str(complex(real, imag)) for real, imag in sublist] for sublist in data['interstate_socs'][(i, j)][soc_text]]
                soc_matrix_dict.update({str(i+1)+'_'+str(j+1): complex_list})
                socc_dict.update({str(i+1)+'_'+str(j+1): data['interstate_socs'][(i, j)][socc_text]})

            except (KeyError, ValueError, TypeError): # If the SOC is 0 or if it is not defined (Singlet-Triplets)
                soc_matrix_dict.update({str(i)+'_'+str(j): 0})
                socc_dict.update({str(i)+'_'+str(j): 0})
            
            try: # If there is angular momentum 
                imaginary_list = [(str(num)+"j") for num in data['interstate_angmom'][(i, j)]['angular_momentum']]

                # Since it is written as ket-bra, the conjugated complex must be obtained 
                imaginary_changed = [str(complex(x).conjugate()) for x in imaginary_list]
                angmoment_dict.update({str(i+1)+'_'+str(j+1): imaginary_changed})
            
            except (KeyError, TypeError): # If the orbital momentum is 0 or if it is not defined (Singlet-Triplets)
                angmoment_dict.update({str(i+1)+'_'+str(j+1): ['0j', '0j', '0j']}) 

    output_dict = {
    "scf_alpha_beta_electrons": data["scf_electrons"], 
    "energy_dict": energy_dict,
    "soc_matrix_dict": soc_matrix_dict,
    "socc_dict": socc_dict, 
    "spin_dict": spin_dict,
    "angmoment_dict": angmoment_dict,
    "transitions_dict": transitions_dict
    }
    
    check_dimmensions(output_dict, ["energy_dict", "spin_dict", "transitions_dict"], totalstates+1)
    check_dimmensions(output_dict, ["soc_matrix_dict", "socc_dict", "angmoment_dict"], (totalstates * (totalstates + 1) // 2))
    return output_dict

file = str(sys.argv[1])

with open(file, encoding="utf8") as f:
        outpuut = f.read()

if bool(outpuut.find('R A S M A N 2')+1): # RAS method output
    output_dict = gtensor_parser_rasci(outpuut)

elif bool(outpuut.find('TDDFT')+1): # TDDFT method output
    output_dict = gtensor_parser_tddft(outpuut)

# Main code block that should not run when imported
if __name__ == "__main__":
    output_json(file, output_dict)