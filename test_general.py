#####################################
#          MODULES SELECTION
#####################################
from parser_gtensor import gfactor_presentation
from parser_excitstates import get_excited_states_analysis
from parser_plots import sos_analysis_and_plot, gfactor_all_states

#####################################
#            INPUT VALUES
#####################################
g_calculation = 1
excited_states_analysis = 0
sos_analysis = 0
gfactor_excited_states = 0

ras_input = '\
roberto_molecules/C74H32B4/C74H32B4_8_8_singletref_concat_triplets.out'
ppms = 1

selected_states = 1  # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [2,1,3,4,5,6,7,8,9,10]  # States to be included when "selected_states = 0"
# [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
symmetry_selection = 'A2'  # Symmetry selected states
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

# OUTPUT
write_ras_input = 0  # 0: write results directly; 1: write in output qchem_file
output_ras_input = ras_input + '-gvalues.txt'
if write_ras_input == 1:
    sys.stdout = open(output_ras_input, "w")

###########################################
#      ACTING IN SEVERAL MOLECULES
###########################################
# several_molecules = 0
# path = "roberto_molecules"
# Import Module
# https://www.geeksforgeeks.org/how-to-read-multiple-text-files-from-folder-in-python/
# if several_molecules == 1:
#     os.chdir(path)  # Change the directory
#
#     g_list = []
#     # iterate through all file_ms_notnull
#     for file_ms_notnull in os.listdir():
#         qchem_file = f"{file_ms_notnull}"
#
#         totalstates = get_number_of_states(qchem_file)
#         states_ras = get_selected_states(qchem_file, totalstates, states_ras, selected_states, symmetry_selection)
#         eigenenergies_ras, excitation_energies_ras = get_eigenenergies(qchem_file, totalstates, states_ras)
#         selected_socs, sz_list, sz_ground = get_spin_orbit_couplings(qchem_file, totalstates, states_ras, soc_options)
#         g_shift = from_energies_soc_to_g_values(qchem_file, states_ras,
#                                                 totalstates, excitation_energies_ras,
#                                                 selected_socs, sz_list, sz_ground)
#
#         qchem_file = qchem_file.replace('.out', '')
#         g_list.append([qchem_file, np.round(g_shift.real[0]*1000, 3),
#                        np.round(g_shift.real[1]*1000, 3), np.round(g_shift.real[2]*1000, 3)])
#
#     sorted_list = sorted(g_list)
#     # for i in range(0, len(sorted_list)):
#     #     print(sorted_list[i, :])
#     print(tabulate(sorted_list, headers=["molecule", "gxx", "gyy", "gzz"]))
#     exit()

#      G-VALUE CALCULATION
if g_calculation == 1:
    gfactor_presentation(ras_input, states_ras, selected_states, symmetry_selection, soc_options, ppms)

#        PLOT ANALYSIS
if excited_states_analysis == 1:
    get_excited_states_analysis(ras_input, cutoff=0.8)

if sos_analysis == 1:
    sos_analysis_and_plot(ras_input, states_ras, selected_states, order_symmetry=1)

if gfactor_excited_states == 1:
    gfactor_all_states(ras_input, states_ras)
