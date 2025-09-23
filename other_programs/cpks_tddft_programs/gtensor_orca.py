import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from projection_method_functions.parsers_OPO.gmatrix_functions import save_picture

ppm = 0
save_options = 0


def convert_time_to_unit(time_list, target_unit):
    """
    Convert a time given in days, hours, minutes, seconds, and milliseconds to a specified unit.

    Parameters:
    - time_list: list [days, hours, minutes, seconds, milliseconds]
    - target_unit: str, target unit for conversion ('days', 'hours', 'minutes', 'seconds', 'milliseconds')

    Returns:
    - float: converted time in the target unit
    """
    # Conversion factors to seconds
    conversion_factors = {
        'days': 86400,          # 1 day = 86400 seconds
        'hours': 3600,          # 1 hour = 3600 seconds
        'minutes': 60,          # 1 minute = 60 seconds
        'seconds': 1,           # 1 second = 1 second
        'milliseconds': 1e-3    # 1 millisecond = 0.001 seconds
    }

    if target_unit not in conversion_factors:
        raise ValueError(f"Invalid target unit '{target_unit}'. Choose from 'days', 'hours', 'minutes', 'seconds', 'milliseconds'.")

    # Convert time components to seconds
    days, hours, minutes, seconds, milliseconds = time_list
    total_seconds = (
        days * conversion_factors['days'] +
        hours * conversion_factors['hours'] +
        minutes * conversion_factors['minutes'] +
        seconds +
        milliseconds * conversion_factors['milliseconds']
    )

    # Convert total seconds to target unit
    return total_seconds / conversion_factors[target_unit]


def gshift_dictionary_plots(filepath, data, gmagnitud, bar_width, saveoptions):
    # Letter size
    smalllet = 10
    medlet = 18
    biglet = 20

    # Extract data
    keys = list(data.keys())
    values = list(data.values())
    x_labels = ['$\mathregular{\Delta g_{xx}}$', '$\mathregular{\Delta g_{yy}}$', '$\mathregular{\Delta g_{zz}}$']
    # x_labels = [f"Elem {i+1}" for i in range(len(values[0]))]  # X-axis labels for each bar group

    # Bar plot parameters
    x = np.arange(len(x_labels)) * 2  # X locations for the groups with separation
    width = bar_width  # Width of the bars

    # Set default font to sans-serif
    plt.rcParams.update({'font.family': 'sans-serif'})

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot each set of bars with a black edge
    for i, key in enumerate(keys):
        ax.bar(x + i * width, values[i], width, label=key, edgecolor='black')

    # Set x-axis labels and title
    ax.set_xticks(x + width * (len(keys) - 1) / 2)  # Align group labels to the middle
    ax.set_xticklabels(x_labels, fontsize=medlet)
    ax.set_ylabel(gmagnitud, fontsize=medlet)
    # ax.set_title("Different contributions to the g-value", fontsize=biglet)
    ax.legend(fontsize=smalllet)

    # Increase the size of tick labels
    ax.tick_params(axis='both', which='major', labelsize=biglet)

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    # plt.show()
    save_picture(saveoptions, filepath, '')


if __name__ == "__main__":

    file_path = str(sys.argv[1])

    with open(file_path, 'r') as file:
        lines = file.readlines()  

    for i in range(len(lines)):
        if 'The g-matrix:' in lines[i]: # Take g-matrix
            gmatrix = []
            for j in range(0, 3):
                gmatrix.append([float(lines[i+1+j].strip().split()[0]), 
                                float(lines[i+1+j].strip().split()[1]), 
                                float(lines[i+1+j].strip().split()[2])])
        
        if 'Breakdown of the contributions' in lines[i]: # Take each contribution to the g-shift
            g_contributions = {}
            for line in lines[i+1:]:
                if '-----' in line:
                    break
                else:
                    term = line.split()[0]
                    gvalues = [float(num) for num in line.split()[1:4]]
                    g_contributions.update({term: gvalues})

        if 'Delta-g' in lines[i]: # Take g-shift
                term = lines[i].split()[0] 
                gvalues = [float(num) for num in lines[i].split()[1:4]]
                g_contributions.update({term: gvalues})
                
        if 'Orientation:' in lines[i]: # Take the orientation matrix
            oriented_list = []
            for j in range(0, 3):
                oriented_list.append([float(lines[i+1+j].strip().split()[1]), 
                                    float(lines[i+1+j].strip().split()[2]),
                                    float(lines[i+1+j].strip().split()[3])])

        if 'TOTAL RUN TIME' in lines[i]: # Take total calculation time 
            tiempo_string = lines[i].replace("TOTAL RUN TIME:","").split()[0::2]
            tiempo = [float(num) for num in tiempo_string]
            coverted_time = convert_time_to_unit(tiempo, 'seconds')

    # Reorder the g-values taking into account the orientation
    max_indices = [sublist.index(max(sublist, key=abs)) for sublist in oriented_list]
    g_contributions_oriented = {}
    for k, v in g_contributions.items():
        gshift = [v[max_indices[i]] for i in range(0,3)]
        g_contributions_oriented.update({k: gshift})

    # Put them in ppt or ppm
    if ppm == 0:
        g_contributions_oriented = {key: [value * 10**3 for value in values] for key, values in g_contributions_oriented.items()}
        g_magnitud = 'ppt'
    elif ppm == 1:
        g_contributions_oriented = {key: [value * 10**6 for value in values] for key, values in g_contributions_oriented.items()}
        g_magnitud = 'ppm'

    #####################################
    ## PRINTING RESULTS 
    #####################################
    print('Delta g contributions in ', g_magnitud)
    for k,v in g_contributions_oriented.items():
        # print(k, v)
        formatted_values = [f"{value:.2f}" for value in v]
        print(f"{k}: {', '.join(formatted_values)}")
    print()
    print('g-matrix:')
    df = pd.DataFrame(gmatrix)
    print(df.to_string(index=False, header=False))
    print()
    print('Orientation matrix:')
    df = pd.DataFrame(oriented_list)
    print(df.to_string(index=False, header=False))
    print()
    print("Time (seconds):", np.round(coverted_time, 3))

    # Plot only the RMC, DSO and PSO 1e and 2e
    del g_contributions_oriented['gel']
    del g_contributions_oriented['gDSO(tot)']
    del g_contributions_oriented['gPSO(tot)']
    gshift_dictionary_plots(file_path.replace('.out',''), 
                        g_contributions_oriented, 
                        '$\mathregular{\Delta g}$, '+g_magnitud, 
                        0.2, 
                        save_options)
