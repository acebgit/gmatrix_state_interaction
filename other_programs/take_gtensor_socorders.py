import os
import glob
import numpy as np
import pandas as pd

# Get the current working directory
current_directory = os.getcwd()

# Get all .txt files in the current directory
txt_files = glob.glob(os.path.join(current_directory, '*.txt'))

list_molecules = []
gmatrix_order0 = {}
gmatrix_order1 = {}

# Loop through all .txt files
for file in txt_files:
    # Ensure we're only reading regular files
    if os.path.isfile(file):
        
        # Open and read the file
        with open(file, 'r') as f:
            lines = f.readlines()
            # print(f"Content of {file}:\n{content}\n")
        
        for line in range(0,len(lines)):
            if "File selected" in lines[line]:
                molecule = lines[line].strip().split()[2].split("_")[0]
                list_molecules.append(molecule)
            elif "SOC:" in lines[line]:
                if "SOC:  All orders" in lines[line]:
                    soc = 0
                if "SOC:  First-order" in lines[line]:
                    soc = 1
            elif "g-factor" in lines[line]:
                if soc == 0:
                    gmatrix_order0[molecule] = lines[line+1].strip()
                elif soc == 1:
                    gmatrix_order1[molecule] = lines[line+1].strip()

molecule_lista = list(set(list_molecules))
presentation_list = []
for molecule in molecule_lista:
    presentation_list.append([molecule, gmatrix_order0[molecule], gmatrix_order1[molecule]])

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
presentation_matrix = np.array(presentation_list)
df = pd.DataFrame(presentation_matrix, index=None,
        columns=['molecule', 'All orders SOC', '1st order SOC'])

# Save DataFrame to a CSV file
df.to_csv('a_soc_orders.csv', index=False)
