import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_path = 'phenal_results.out'  # Replace with your file path

hybrid_funct = ['B1LYP', 'B3LYP', 'O3LYP', 'X3LYP', 'B1P', 'B3P', 'B3PW', 'PW1PW', 'mPW1PW', 'mPW1LYP', 'PBE0', 'PW6B95', 'BHANDHLYP']
range_hybrid_funct = ['wB97', 'wB97X', 'wB97X-D3', 'wB97X-V', 'wB97X-D3BJ', 'wB97M-V', 'wB97M-D3BJ', 'CAM-B3LYP', 'LC-BLYP']

# Read the data from the file
data = pd.read_csv(file_path, delim_whitespace=True, header=None)


def plots_normalized_gtensors(x, y1, y2, y3, time, file_path, word_added):
    # Set positions of bar on X axis
    r1 = np.arange(len(x))
    r2 = [i + bar_width for i in r1]
    r3 = [i + bar_width for i in r2]

    # CREATE PLOT 1
    plt.figure(figsize=(10, 6))
    plt.bar(r1, y1, color='blue', width=bar_width, edgecolor='black', label=r'$\mathregular{\Delta g_{xx}}$')
    plt.bar(r2, y2, color='green', width=bar_width, edgecolor='black', label=r'$\mathregular{\Delta g_{yy}}$')
    plt.bar(r3, y3, color='red', width=bar_width, edgecolor='black', label=r'$\mathregular{\Delta g_{zz}}$')

    # Add xticks on the middle of the group bars
    plt.xlabel('Hybrid functional', fontweight='bold')
    plt.xticks([r + bar_width for r in range(len(x))], x)
    plt.ylabel('g-shifts normalized', fontweight='bold')
    plt.title(file_path+word_added)

    plt.legend()
    plt.grid(axis='y')

    figure_name = file_path.replace('.out', "_"+word_added+'_gshifts.png')
    plt.savefig(figure_name, dpi=300)
    plt.show()

    # CREATE PLOT 2
    plt.xlabel('Hybrid functional', fontweight='bold')
    plt.ylabel('Time (min)', fontweight='bold')
    plt.title(file_path+word_added)

    plt.plot(x, time, color='black', marker='o', label='Column 5')
    plt.grid(True)  # Add grid to both axes (x and y)

    figure_name = file_path.replace('.out', "_"+word_added+'_times.png')
    plt.savefig(figure_name, dpi=300)
    plt.show()


# Assign columns to variables
x = data[0]
y1_0 = data[1]
y2_0 = data[2]
y3_0 = data[3]
time = data[4]

# Normalized data
y1 = y1_0/np.abs(y1_0).max()
y2 = y2_0/np.abs(y2_0).max()
y3 = y3_0/np.abs(y3_0).max()

# Define the bar width
bar_width = 0.2

x_hyb = []
y1_hyb = []
y2_hyb = []
y3_hyb = []
time_hyb = []

x_rang = []
y1_rang = []
y2_rang = []
y3_rang = []
time_rang = []

for i in range(0, len(x)): 
    if x[i] in hybrid_funct:
        x_hyb.append(x[i])
        y1_hyb.append(y1[i])
        y2_hyb.append(y2[i])
        y3_hyb.append(y3[i])
        time_hyb.append(time[i])
 
    elif x[i] in range_hybrid_funct:
        x_rang.append(x[i])
        y1_rang.append(y1[i])
        y2_rang.append(y2[i])
        y3_rang.append(y3[i])
        time_rang.append(time[i])

print(data.to_string(index=False))

plots_normalized_gtensors(x_hyb, y1_hyb, y2_hyb, y3_hyb, time_hyb, file_path, 'hybrid')

plots_normalized_gtensors(x_rang, y1_rang, y2_rang, y3_rang, time_rang, file_path, 'range_hybrid')
