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






