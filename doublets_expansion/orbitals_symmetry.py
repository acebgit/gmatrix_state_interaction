#####################################
#          ORBITALS SYMMETRY
#####################################
# module load Python/3.7.4-Anaconda3-2019.10

import numpy as np

file = '../\
RASCI_results/cucl4_2-/cucl4_2-_def2tzvp_17_10_20_states.out'

word_search = ['-- Doubly Occupied --']
word_stop = '---------------------------------'
orbitals_symmetry = []
my_orbitals = [34, 35, 36, 37, 41, 43, 44, 45, 50, 51]

with open(file, encoding="utf8") as data:
    for line in data:
        if any(i in line for i in word_search):
            next_line = next(data)

            while word_stop not in next_line:
                if '--' in next_line:
                    next_line = next(data)

                next_line = next(data)  # Skip energies
                split_line = next_line.split()
                for i in range(0, len(split_line), 2):
                    element = ''.join(split_line[i:i+2])
                    orbitals_symmetry.append(element)
                next_line = next(data)

my_symmetry = []
for orbital in my_orbitals:
    my_symmetry.append(orbitals_symmetry[orbital-1])
print(my_symmetry)
