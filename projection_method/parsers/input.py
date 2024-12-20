import sys 
from projection_method.parsers.gtensor_calculation import extract_data_from_json, get_selected_dict, \
    from_json_to_matrices, select_soc_order, from_matrices_to_gshift, print_g_calculation, gtensor_state_pairs_analysis, \
    sum_over_state_plot, from_gvalue_to_shift, comparison_s2

# from_gvalue_to_shift([1.9978])
# exit()

######## G-TENSOR CALCULATION ########
calculate_gshift = 1
ppm = 0 # 0: ppt; 1: ppm
state_selection = 1 # 0: use "state_ras" ; 1: use all states_selected ; 2: use states_selected by selected symmetry

initial_states = list(range(1, 31))
symmetry_selection = 'B2'  # Symmetry selected states_selected
soc_options = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix
soc_orders = 0 # 0: All orders; 1: First-order; 2: Higher-order 

######## G-TENSOR CALCULATION BY PAIRS ########
cutoff_gvalue = 0 # ≠0: cut-off between ground-excited states (% of maximum g-value in each dim) 
cutoff_config = 0 # cut-off for configurations amplitude (% of maximum amplitude)
excit_plot = 0 # 0: not show plot, 1: show plot 

######## SOS PLOTS ########
sum_over_states_analysis = 0 # SOS g-tensor plot: g-tensor calculation with n states
sos_cutoff = 0
g_estimation = 0
save_plot = 0

#### <S2> comparison ####
s2_comparison = 0

file = str(sys.argv[1]) + ".json"

output_dict = extract_data_from_json(file)

output_dict_selected = get_selected_dict(output_dict, len(output_dict["energy_dict"]), initial_states, state_selection, symmetry_selection)

# Organize: g_shift class, with the three methods of the three options
# In the class, give as input the input parameters at the beginning 
if calculate_gshift == 1:
    sz_ground, max_sz_list, approxspin_list, matrices_dict = from_json_to_matrices(output_dict_selected)

    matrices_dict = select_soc_order(matrices_dict, soc_orders)

    gmatrix, gshift = from_matrices_to_gshift(max_sz_list, sz_ground, matrices_dict, ppm)

    print_g_calculation(file, output_dict_selected, approxspin_list, gshift, gmatrix, ppm, soc_options, soc_orders)

if cutoff_gvalue != 0:
    gtensor_state_pairs_analysis(output_dict_selected, ppm, cutoff_gvalue, cutoff_config, excit_plot, savepicture=0)

if sum_over_states_analysis == 1:
    sum_over_state_plot(output_dict_selected, g_estimation, ppm, sos_cutoff, save_plot, showplot=1)

if s2_comparison == 1:
    # Get the states with the largest influence in the g-tensor
    gshift_list = sum_over_state_plot(output_dict_selected, g_estimation, ppm, sos_cutoff, save_plot, showplot=0)
    states_gtensor = [1] # First state is always added
    states_gtensor.extend(sublist[0] for sublist in gshift_list if sublist)

    # Form the output_dict with the previous states
    output_dict_selected = get_selected_dict(output_dict, len(output_dict["energy_dict"]), states_gtensor, states__selection=0, sym__selected=0)

    # Obtain the <S2> analysis
    comparison_s2(file, output_dict_selected, save_plot)
