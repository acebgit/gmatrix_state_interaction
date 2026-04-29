def save_picture(save_options, filee, title_main):
    """
    Function that shows the plot (save_options=0) or save it (save_options=1).
    :param save_options, file, main_title.
    :return: plot (saved or shown)
    """
    if save_options == 1:
        plt.plot()
        figure_name = filee + '_' + title_main + '.png'
        plt.savefig(figure_name, 
                    bbox_inches="tight", 
                    dpi=600
                    )
        plt.close()
    else:
        plt.plot()
        plt.show()


def plot_g_tensor_vs_states(file, subtitle, presentation_matrix, x_title, y_title, main_title, save_options):
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_type = 1 # 0: plot, 1: bars

    # MAIN FEATURES:
    fuente = 'sans-serif'  # 'serif'
    small_size = 17 # 25
    legend_size = small_size + 3
    bigger_size = small_size + 3
    weight_selected = 'normal'
    line_width = 2
    marker_size = 10

    x = presentation_matrix[:, 0]  # First column for x-axis
    y1 = presentation_matrix[:, 1]  # Second column for the first category
    y2 = presentation_matrix[:, 2]  # Third column for the second category
    y3 = presentation_matrix[:, 3]  # Fourth column for the third category

    #################################
    ###   PLOT TYPE
    #################################
    if plot_type == 0:
        # MAJOR AND MINOR TICKS:
        # x_tick = int((max(x))) / 4
        # x_tick = int((max(x))) / 4
        # y_tick = int((max(x))) / 4
        # x_tick = 20
        # y_tick = 1
        # ax.xaxis.set_major_locator(MultipleLocator(x_tick))
        # ax.yaxis.set_major_locator(MultipleLocator(y_tick))

        # x_tick_min = x_tick / 2
        # y_tick_min = y_tick / 2
        # ax.xaxis.set_minor_locator(MultipleLocator(x_tick_min))
        # ax.yaxis.set_minor_locator(MultipleLocator(y_tick_min))

        # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        # LINES:
        # ax.plot(x, y1, 'r',
        #         label=r'$\mathregular{\Delta g_{xx}}$', linewidth=line_width)
        # ax.plot(x, y2, 'b', 
        #         label=r'$\mathregular{\Delta g_{yy}}$', linewidth=line_width)
        # ax.plot(x, y3, 'k',
        #         label=r'$\mathregular{\Delta g_{zz}}$', linewidth=line_width)

        # MARKERS
        ax.plot(x, y1, 'ro', markersize=marker_size)
        ax.plot(x, y2, 'bv', markersize=marker_size,
                markerfacecolor='none', markeredgewidth=1.5)
        ax.plot(x, y3, 'ks', markersize=marker_size)

        # Enable grid lines for both x and y axes
        plt.grid(axis='both', linestyle='--', linewidth=0.7, alpha=0.7)

    #################################
    ###   BAR PLOTS
    #################################
    elif plot_type == 1:
        # Set width of the bars
        bar_width = 0.25

        # Set the positions of the bars on the x-axis
        r1 = x  # np.arange(len(x))
        r2 = [x + bar_width for x in r1]
        r3 = [x + bar_width * 2 for x in r1]

        # Create the bar plot 
        plt.bar(r1, y1, width=bar_width, color='red', edgecolor='red', label=r'$\mathregular{\Delta g_{xx}}$')
        plt.bar(r2, y2, width=bar_width, color='blue', edgecolor='blue', label=r'$\mathregular{\Delta g_{yy}}$')
        plt.bar(r3, y3, width=bar_width, color='black', edgecolor='black', label=r'$\mathregular{\Delta g_{zz}}$')

    # CHANGING THE FONTSIZE OF TICKS
    plt.xticks(fontsize=small_size, weight=weight_selected)
    plt.yticks(fontsize=small_size, weight=weight_selected)
    # axis.set_major_locator(MaxNLocator(integer=True))

    # LIMIT TO AXIS:
    ax.set_xlim(min(x)-1, max(x)+1)

    # To put only integer numbers in the label: 
    # from matplotlib.ticker import MaxNLocator
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    max_value = np.maximum.reduce([y1.max(), y2.max(), y3.max()])
    min_value = np.minimum.reduce([y1.min(), y2.min(), y3.min()])
    # Choose the larger value in absolute terms to calculate the interval
    interval = 0.1 * max(abs(max_value), abs(min_value))
    ax.set_ylim(min_value-interval, max_value+interval)

    # LABELS:
    # labelpad: change the space between axis umbers and labels
    plt.xlabel(x_title, fontsize=bigger_size, fontfamily=fuente, labelpad=15,
               weight=weight_selected)
    plt.ylabel(y_title, fontsize=bigger_size, fontfamily=fuente, style='italic',
               weight=weight_selected, labelpad=15)
    # x_min = 0
    # x_max =  11
    # y_min = -45
    # y_max =  5
    # plt.xlim([x_min, x_max])  # Limit axis values
    # plt.ylim([y_min, y_max])  # Limit axis values

    # TITLE:
    # y = 1.05 change the space between title and plot
    # plt.title(main_title, fontsize=bigger_size, fontfamily=fuente, y=1.05)
    
    # LEGEND
    legend = plt.legend(fontsize=legend_size, 
                        fancybox=True, 
                        framealpha=0.5,
                        labelcolor='linecolor', 
                        loc='best', 
                        frameon=False, 
                        )
    # plt.legend(frameon=False)
    # frame = legend.get_frame().set_edgecolor('black')
    # frame = legend.get_frame().set_linewidth(1)
    # frame = legend.get_frame().set_facecolor('white')
    # frame.set_edgecolor('black')

    # plt.locator_params(nbins=10)
    # plt.grid()

    # Add an horizontal line in y = 0
    # ax.hlines(y=0, xmin=x_min, xmax=x_max, linewidth=line_width, color='k',
    #           linestyle='dotted')
    # dotted, dashed, solid, dashdot

    line_width = line_width - 0.8
    ax.spines["top"].set_linewidth(line_width)
    ax.spines["bottom"].set_linewidth(line_width)
    ax.spines["left"].set_linewidth(line_width)
    ax.spines["right"].set_linewidth(line_width)
    
    save_picture(save_options, file, subtitle)


def sum_over_state_plot(outputdict, 
                        gestimation, 
                        ppm, 
                        sos_cutoff, 
                        sos_save_plot):
    """
    Generate the sum-over-states plot, i.e. calculation of the g-tensor by including states
    from 1 to nstates, above a cutoff of g-value. 
    This can be done by:
    - gestimation = 0: effective Hamiltonian created between each pair of states
    - gestimation = 1: use of predictive phormula 
    :param: 
    :return: shows SOS plot
    """
    def filter_dictionary(dictionary, state):
        """
        Form a dictionary with only a pair of states: ground state and "state"
        """
        new_dict = {}
        ground_state = list(outputdict["energy_dict"].keys())[0]
        
        for name_dict in ["energy_dict", "spin_dict"]:
            new_dict[name_dict] = {k: v for k, v in dictionary[name_dict].items() if k == ground_state or k == state}

        for name_dict in ["soc_matrix_dict", "angmoment_dict"]:
            for k, v in dictionary[name_dict].items():
                if k in [f"{ground_state}_{state}", f"{state}_{ground_state}"]:
                    new_dict[name_dict] = {k: v}
        return new_dict
    
    filtered_gshifts = []

    if gestimation == 0:
        all_gshifts = []
        
        for excit_state in list(outputdict["energy_dict"].keys())[1:]:
            
            # Form a dictionary only for a pair of states
            filtered_dict = filter_dictionary(outputdict, excit_state)
            
            states__lengthsz, approxspin_dict, matrices_dict = from_json_to_matrices(filtered_dict)

            gmatrix, gshift = from_matrices_to_gshift(states__lengthsz, matrices_dict, ppm)
            
            all_gshifts.append([excit_state, 
                                (np.round(gshift[0].real, 3)),
                                (np.round(gshift[1].real, 3)), 
                                (np.round(gshift[2].real, 3))])

        # Convert to a NumPy array for efficient processing
        all_gshifts_array = np.array(all_gshifts, dtype=np.float64)

        # Step 1: Find the three maximum values in each of the last three columns and multiply by cutoff
        cutoffs = np.max(np.abs(all_gshifts_array[:, 1:]), axis=0) * sos_cutoff

        # Step 2: Filter rows where at least one column meets or exceeds the threshold
        data = [
            [int(row[0])] + list(row[1:])  # Convert the first value to an integer and keep the rest as float
            for row in all_gshifts_array
            if any(abs(row[i]) >= cutoffs[i-1] for i in range(1, 4))
        ]

        # Convert np.float64 to regular float
        filtered_gshifts = [
            [row[0]] + [float(value) for value in row[1:]]  # Convert all elements except the first to floats
            for row in data
        ]

    elif gestimation == 1:
        gshift_dict = gshift_estimation_loop(outputdict, ppm)

        # Avoid to include g-tensor with zero
        threshold = 10**(-6)
        new_cutoff = sos_cutoff if sos_cutoff != 0 else threshold

        # Create the new dictionary with the maximum absolute value multiplied by cutoff
        cut_gvalues = {key: max((abs(v), v) for v in values.values())[1] * new_cutoff for key, values in gshift_dict.items()}
        
        # Take states with estimated g-shift higher than a cutoff
        for k, v in gshift_dict["gxx"].items():
            if any(abs(gshift_dict[key][k]) >= cut_gvalues[key] and cut_gvalues[key] >= threshold
                for key in ["gxx", "gyy", "gzz"]):
                    filtered_gshifts.append([int(k),
                                             (gshift_dict["gxx"][k]),
                                             (gshift_dict["gyy"][k]),
                                             (gshift_dict["gzz"][k])])

    print("------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("------------------------------")
    technique = {0: "Effective Hamiltonian", 1: "Estimation phormula"}
    print("Technique used:", technique.get(gestimation, "Unknown"))
    print("cut-offs g-value (%): ", sos_cutoff)

    # Set display options to show all rows and columns
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.DataFrame([row[0:4] for row in filtered_gshifts], [row[0] for row in filtered_gshifts], columns=['state','gxx','gyy','gzz'])
    print(df.to_string(index=False))

    # Compute the sum of all the elements
    perturbative_sum = [round(float(sum(x)), 3) for x in zip(*filtered_gshifts)][1:]  # Sum columns, skipping the first
    print('Total: ', perturbative_sum[0], perturbative_sum[1], perturbative_sum[2])

    # Make the title
    y_title = r'$\Delta g, ppt$' if ppm == 0 else r'$\Delta g, ppm$'
    # file_string = "str(sys.argv[1]).split('.')[0]
    # plot_title = 'sos_analysis: ' + file_string
    plot_g_tensor_vs_states(matrix=np.array(filtered_gshifts, dtype=object),
                            y_title=y_title,
                            save_option=sos_save_plot)
    return filtered_gshifts






