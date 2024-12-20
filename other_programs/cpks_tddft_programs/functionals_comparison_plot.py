import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from other_programs.cpks_tddft_programs.gtensor_orca import gshift_dictionary_plots

ppm = 0
save_options = 1

# Get the current directory
current_directory = os.getcwd()

# List all files in the current directory
files = [f for f in os.listdir(current_directory) if f.endswith('.txt')]

# Take g-values of each of the files and put them in list of dictionaries
gresults = []
all_molecules = []
all_methods = []

for file in files:
    molecule = file.split('_')[0]
    method = file.split('_')[1]
    functional = file.split('_')[2]

    all_molecules.append(molecule)
    all_methods.append(method)

    with open(file, 'r') as f:
        for line in f:
            if 'Delta-g' in line:
                gshift = [float(num) for num in re.findall(r'-?\d+\.\d+', line)]
            elif 'g-factor' in line:
                gshift = [float(num) for num in re.findall(r'-?\d+\.\d+', line)]
            
        gresults.append({'Molecule': molecule, 
                            'Method': method ,
                            'Functional': functional, 
                            'gshift': gshift})

all_molecules = list(set(all_molecules))
all_methods = list(set(all_methods))

g_magnitud = 'ppt' if ppm == 0 else 'ppm'

# Take all the values of each of the keys 
for molecule in all_molecules:
    for method in all_methods:
        
        # Filter dictionaries for a concrete molecule and method 
        dict_molecule = [d for d in gresults if d.get("Molecule") == molecule]
        dict_method = [d for d in dict_molecule if d.get("Method") == method]

        # Form a list of dictionaries to form the plot 
        plot_dict_1 = {}
        for element in dict_method: 
            plot_dict_1.update({element["Functional"]: element["gshift"]})

        # Sort the dictionary by keys alphabetically
        plot_dict = {key: plot_dict_1[key] for key in sorted(plot_dict_1, key=str.lower)}

        # Print results
        print('Molecule:', molecule, ', method: ', method)
        df = pd.DataFrame.from_dict(plot_dict, orient='index', columns=['gxx', 'gyy', 'gzz'])
        print(df.to_string(index=True, header=True))
        print()

        # Plot the values
        gshift_dictionary_plots(molecule + "_" + method + "_functionals", 
                                plot_dict, 
                                '$\mathregular{\Delta g}$, ' + g_magnitud, 
                                0.1, 
                                save_options)
