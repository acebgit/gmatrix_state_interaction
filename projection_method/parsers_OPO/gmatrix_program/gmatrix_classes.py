__author__ = 'Antonio Cebreiro'

import sys
import json
from PyQchem.pyqchem.parsers.parser_rasci import parser_rasci
from PyQchem.pyqchem.parsers.parser_cis import basic_cis
from projection_method.parsers_OPO.gmatrix_functions import (
    extract_data_from_json, get_selected_dict,
    from_json_to_matrices, from_matrices_to_gshift, print_g_calculation,
    gtensor_state_pairs_analysis, sum_over_state_plot,
    comparison_s2, get_scaling_analysis
)

class OutToJsonConverter:
    """Convert a Q-Chem .out file to a structured .json file."""

    def __init__(self, filename, cutoff_amp=0):
        self.filename = filename
        self.cutoff_amp = cutoff_amp
        self.outpuut = None
        self.output_dict = None

    # -----------------------------
    # File Handling
    # -----------------------------
    def read_file(self):
        """Read the content of the .out file."""
        with open(self.filename, encoding="utf8") as f:
            self.outpuut = f.read()

    def save_json(self):
        """Save the parsed output to a JSON file."""
        outfile_name = self.filename + ".json"
        with open(outfile_name, 'w') as f:
            json.dump(self.output_dict, f, separators=(',', ':'), sort_keys=False, indent=4)

    # -----------------------------
    # Helper Methods
    # -----------------------------
    @staticmethod
    def check_dimensions(dictionary, keys, dim):
        for key in keys:
            if dictionary[key] and len(dictionary[key]) != dim:
                raise ValueError(f"The dimensions of the dictionary {key} are incorrect.")

    @staticmethod
    def add_count_to_elements(lst):
        count_dict = {}
        new_list = []
        for element in lst:
            count_dict[element] = count_dict.get(element, 0) + 1
            new_list.append(f"{count_dict[element]}{element}")
        return new_list

    @staticmethod
    def s_to_s2(s):
        """Convert spin multiplicity string to S² value."""
        mapping = {'Singlet': 0.0, 'Triplet': 2.0, 'Quintet': 6.0, 'Heptet': 12.0}
        return mapping.get(s, 0.0)

    # -----------------------------
    # Parsing
    # -----------------------------
    def parse(self):
        """Detect method type and parse accordingly."""
        if "R A S M A N 2" in self.outpuut:
            self.output_dict = self._parse_rasci()
        elif "TDDFT" in self.outpuut:
            self.output_dict = self._parse_tddft()
        else:
            raise ValueError("Unrecognized output type in file.")

    # -----------------------------
    # RASCI Parsing
    # -----------------------------
    def _parse_rasci(self):
        data = parser_rasci(self.outpuut)
        excitstates_dict = data['excited_states']
        totalstates = len(excitstates_dict)
        states_symmetry = self.add_count_to_elements([s['symmetry'] for s in excitstates_dict])

        # Energy and spin dictionaries
        energy_dict = {states_symmetry[i]: [s['total_energy'], s['excitation_energy']]
                       for i, s in enumerate(excitstates_dict)}
        spin_dict = {states_symmetry[i]: s['multiplicity'] for i, s in enumerate(excitstates_dict)}

        # Transitions dictionary
        transitions_dict = self._get_ras_transitions(data, totalstates)

        # SOC, SOCC, angular momentum dictionaries
        soc_matrix_dict, socc_dict, angmoment_dict = self._get_interstate_properties_rasci(data, states_symmetry, totalstates)

        output = {
            "scf_alpha_beta_electrons": data["scf_electrons"],
            "energy_dict": energy_dict,
            "soc_matrix_dict": soc_matrix_dict,
            "socc_dict": socc_dict,
            "spin_dict": spin_dict,
            "angmoment_dict": angmoment_dict,
            "transitions_dict": transitions_dict
        }

        # Dimension checks
        self.check_dimensions(output, ["energy_dict", "spin_dict"], totalstates)
        self.check_dimensions(output, ["soc_matrix_dict", "socc_dict", "angmoment_dict"], totalstates*(totalstates-1)//2)
        self.check_dimensions(output, ["transitions_dict"], totalstates)
        return output

    def _get_ras_transitions(self, data, totalstates):
        """Compute RASCI transitions."""
        # Helper: map RAS → SCF order
        def from_ras_to_scf_order(beta_homo, ras_elec_orbs, initial_orbitals):
            ras_scf_map = {}
            for scf_orbital in range(1, beta_homo + 1):
                if scf_orbital not in ras_elec_orbs:
                    ras_scf_map[len(ras_scf_map)+1] = scf_orbital
            for scf_orbital in ras_elec_orbs:
                ras_scf_map[len(ras_scf_map)+1] = scf_orbital
            for scf_orbital in range(beta_homo + 1, beta_homo + 100):
                if scf_orbital not in ras_elec_orbs:
                    ras_scf_map[len(ras_scf_map)+1] = scf_orbital
            return ','.join(str(ras_scf_map[i]) for i in initial_orbitals)

        # Identify active orbitals
        try:
            enum = self.outpuut.find('RAS_ACT_ORB')
            ras_elec_orb = [int(x) for x in (self.outpuut[enum:enum+200].split('[')[1].split(']')[0]).split(",")]
        except:
            raise RuntimeError("RAS_ACT_ORB selection not implemented automatically")

        transitions = []
        for state in range(totalstates):
            state_transitions = []
            for config in data['excited_states'][state]['configurations']:
                amp = float(config['amplitude'])
                if abs(amp) >= self.cutoff_amp:
                    alpha = config['occupations']['alpha']
                    beta = config['occupations']['beta']
                    elec_beta = beta.count(1)
                    somo_ras = [i+1 for i in range(len(alpha)) if alpha[i] != beta[i]]
                    somo_scf = from_ras_to_scf_order(elec_beta, ras_elec_orb, somo_ras)
                    config_index = data['excited_states'][state]['configurations'].index(config)
                    state_transitions.append({
                        'state': state+1,
                        'configuration': config_index,
                        'SOMO': somo_scf,
                        'amplitude': round(amp, 3)
                    })
            transitions.append(state_transitions)
        return transitions

    def _get_interstate_properties_rasci(self, data, states_symmetry, totalstates):
        """Compute SOC, SOCC, and angular momentum for RASCI."""
        soc_matrix_dict = {}
        socc_dict = {}
        angmoment_dict = {}
        soc_text = 'total_soc_mat'
        interstates_dict = data['interstate_properties']

        for i in range(1, totalstates+1):
            for j in range(1, totalstates+1):
                bra, ket = states_symmetry[i-1], states_symmetry[j-1]
                braket, ketbra = bra+"_"+ket, ket+"_"+bra
                if i != j and braket not in soc_matrix_dict and ketbra not in soc_matrix_dict:
                    soc_matrix_dict[braket] = [[str(e) for e in sub] for sub in interstates_dict[(i, j)][soc_text]]
                    socc_dict[braket] = interstates_dict[(i, j)]['mf_socc']
                    angmoment_dict[braket] = [str(e) for e in interstates_dict[(i, j)]['angular_momentum']]
        return soc_matrix_dict, socc_dict, angmoment_dict

    # -----------------------------
    # TDDFT Parsing
    # -----------------------------
    def _parse_tddft(self):
        data = basic_cis(self.outpuut)
        excitstates_dict = data['excited_states']
        totalstates = len(excitstates_dict)

        # Energy dictionary
        energy_dict = {1: [data['scf_energy'], 0.0]}
        for i, s in enumerate(excitstates_dict):
            energy_dict[i+2] = [s['total_energy'], s['excitation_energy']]

        # Spin dictionary
        spin_dict = {1: data['scf_multiplicity']}
        for i, s in enumerate(excitstates_dict):
            try:
                spin_dict[i+2] = float(s['multiplicity'])
            except:
                spin_dict[i+2] = self.s_to_s2(s['multiplicity'])

        transitions_dict = self._get_tddft_transitions(data, totalstates)
        soc_matrix_dict, socc_dict, angmoment_dict = self._get_interstate_properties_tddft(data, totalstates)

        output = {
            "scf_alpha_beta_electrons": data["scf_electrons"],
            "energy_dict": energy_dict,
            "soc_matrix_dict": soc_matrix_dict,
            "socc_dict": socc_dict,
            "spin_dict": spin_dict,
            "angmoment_dict": angmoment_dict,
            "transitions_dict": transitions_dict
        }

        # Dimension checks
        self.check_dimensions(output, ["energy_dict", "spin_dict", "transitions_dict"], totalstates+1)
        self.check_dimensions(output, ["soc_matrix_dict", "socc_dict", "angmoment_dict"], (totalstates*(totalstates+1))//2)
        return output

    def _get_tddft_transitions(self, data, totalstates):
        """Compute TDDFT transitions."""
        scf_alpha, scf_beta = data['scf_electrons']
        scf_somos = ",".join(map(str, range(min(scf_alpha, scf_beta)+1, max(scf_alpha, scf_beta)+1)))
        transitions = [[{'state': 1, 'configuration': 1, 'SOMO': scf_somos, 'amplitude': 1.0}]]

        for state in range(totalstates):
            state_transitions = []
            for config in data['excited_states'][state]['configurations']:
                amp = float(config['amplitude'])
                if abs(amp) >= self.cutoff_amp:
                    alpha = config['occupations']['alpha']
                    beta = config['occupations']['beta']
                    somo = ','.join(str(i+1) for i in range(len(alpha)) if alpha[i] != beta[i])
                    config_index = data['excited_states'][state]['configurations'].index(config)
                    state_transitions.append({
                        'state': state+2,
                        'configuration': config_index+1,
                        'SOMO': somo,
                        'amplitude': round(amp, 3)
                    })
            transitions.append(state_transitions)
        return transitions

    def _get_interstate_properties_tddft(self, data, totalstates):
        soc_matrix_dict = {}
        socc_dict = {}
        angmoment_dict = {}

        soc_text, socc_text = 'mf_soc_mat', 'mf_socc'

        for i in range(totalstates+1):
            for j in range(i+1, totalstates+1):
                try:
                    soc_matrix_dict[f"{i+1}_{j+1}"] = [[str(complex(r, im)) for r, im in sub]
                                                        for sub in data['interstate_socs'][(i, j)][soc_text]]
                    socc_dict[f"{i+1}_{j+1}"] = data['interstate_socs'][(i, j)][socc_text]
                except:
                    soc_matrix_dict[f"{i+1}_{j+1}"] = 0
                    socc_dict[f"{i+1}_{j+1}"] = 0
                try:
                    angmoment_dict[f"{i+1}_{j+1}"] = [str(complex(str(num)+"j").conjugate()) 
                                                       for num in data['interstate_angmom'][(i, j)]['angular_momentum']]
                except:
                    angmoment_dict[f"{i+1}_{j+1}"] = ['0j', '0j', '0j']
        return soc_matrix_dict, socc_dict, angmoment_dict


class GTensorConfig:
    """Holds configuration parameters for g-tensor analysis."""
    def __init__(self):
        # --- State selection ---
        self.state_selection = 1  # 0: use "state_ras"; 1: use all states_selected; 2: use states_selected by symmetry
        self.initial_states = [1, 2]  # Example: list(range(1, 12)) ; can remove some manually
        self.symmetry_selection = 'B2'  # Symmetry for state selection

        # --- g-tensor calculation ---
        self.calculate_gshift = 1
        self.ppm = 0  # 0: ppt; 1: ppm
        self.soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix
        self.soc_orders = 0  # 0: All orders; 1: First-order only

        # --- g-tensor calculation by pairs ---
        self.cutoff_gvalue = 0  # !=0: cut-off between ground/excited states (% of max g-value)
        self.cutoff_config = 0  # cut-off for configuration amplitude (% of max amplitude)
        self.excit_plot = 0  # 0: do not show plot; 1: show plot

        # --- Sum-over-states (SOS) analysis ---
        self.sum_over_states_analysis = 0  # SOS g-tensor plot with n states
        self.sos_cutoff = 0
        self.g_estimation = 0
        self.save_plot = 1

        # --- <S²> analysis ---
        self.s2_comparison = 0
        self.s2_save_plot = 1
        self.s2_sos_cutoff = 0
        self.s2_per_gshift = 1

        # --- Scaling analysis ---
        self.scaling_analysis = 1  # 1: SOC; 2: Orbital; 3: Spin; 4: Orbital + SOC
        self.fit_degree = 3  # Polynomial fitting degree
        self.scaling_save_plot = 1


class GTensorPipeline:
    """Main pipeline for g-tensor analysis."""
    def __init__(self, file, config=None):
        # Automatically append .json if not present
        self.file = file if file.endswith(".json") else file + ".json"
        self.config = config if config else GTensorConfig()
        self.output_dict = None
        self.output_dict_selected = None

    def load_data(self):
        """Load raw JSON data."""
        self.output_dict = extract_data_from_json(self.file)

    def select_states(self):
        """Filter states according to configuration (initial states, symmetry, etc.)."""
        self.output_dict_selected = get_selected_dict(
            self.output_dict,
            len(self.output_dict["energy_dict"]),
            self.config.initial_states,
            self.config.state_selection,
            self.config.symmetry_selection
        )

    def run_gshift_calculation(self):
        """Perform g-shift calculation and print results."""
        states_lengthsz, approxspin_dict, matrices_dict = from_json_to_matrices(self.output_dict_selected)
        gmatrix, gshift = from_matrices_to_gshift(states_lengthsz, matrices_dict, self.config.ppm)
        print_g_calculation(
            self.file,
            approxspin_dict,
            gshift,
            gmatrix,
            self.config.soc_options,
            self.config.soc_orders,
            self.config.ppm
        )

    def run_pairs_analysis(self):
        """Perform g-tensor analysis by pairs of states (optional)."""
        gtensor_state_pairs_analysis(
            self.output_dict_selected,
            self.config.ppm,
            self.config.cutoff_gvalue,
            self.config.cutoff_config,
            self.config.excit_plot,
            savepicture=0
        )

    def run_sos_analysis(self):
        """Perform sum-over-states (SOS) analysis."""
        sum_over_state_plot(
            self.output_dict_selected,
            self.config.g_estimation,
            self.config.ppm,
            self.config.sos_cutoff,
            self.config.save_plot,
        )

    def run_s2_comparison(self):
        """Perform <S²> comparison and spin contamination analysis."""
        gshift_list = sum_over_state_plot(
            self.output_dict_selected,
            self.config.g_estimation,
            self.config.ppm,
            self.config.s2_sos_cutoff,
            saveplot=1
        )
        states_gtensor = [1]  # Ground state always included
        states_gtensor.extend(sublist[0] for sublist in gshift_list if sublist)

        # Re-select dictionary with the new state list
        self.output_dict_selected = get_selected_dict(
            self.output_dict,
            len(self.output_dict["energy_dict"]),
            states_gtensor,
            states__selection=0,
            sym__selected=0
        )

        comparison_s2(
            self.file,
            self.output_dict_selected,
            self.config.s2_save_plot,
            self.config.s2_per_gshift,
            gshift_list
        )

    def run_scaling_analysis(self):
        """Perform SOC/orbital/spin scaling analysis."""
        get_scaling_analysis(
            self.config.scaling_analysis,
            self.config.fit_degree,
            self.output_dict_selected,
            self.file,
            self.config.scaling_save_plot,
            self.config.ppm
        )

    def run(self):
        """Run the full g-tensor workflow according to configuration."""
        self.load_data()
        self.select_states()

        if self.config.calculate_gshift == 1:
            self.run_gshift_calculation()

        if self.config.cutoff_gvalue != 0:
            self.run_pairs_analysis()

        if self.config.sum_over_states_analysis == 1 and self.config.s2_comparison == 0:
            self.run_sos_analysis()

        if self.config.s2_comparison == 1:
            self.run_s2_comparison()

        if self.config.scaling_analysis != 0:
            self.run_scaling_analysis()


if __name__ == "__main__":
    # Example: python program.py example
    file = sys.argv[1]
    converter = OutToJsonConverter(file)
    converter.read_file()
    converter.parse()
    converter.save_json()

    # Example: python program.py example
    file = str(sys.argv[1])
    pipeline = GTensorPipeline(file)
    pipeline.run()
