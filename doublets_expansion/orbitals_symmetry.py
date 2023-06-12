#####################################
#          ORBITALS SYMMETRY
#####################################
# module load Python/3.7.4-Anaconda3-2019.10

import numpy as np

file = '../\
../../Desktop/1_gfactor/calcs/transition_metal_complexes_scf/vo_h2o5_2+_scf_def2tzvp_1.out'
my_orbitals = [33, 37, 38, 39, 40, 42, 43, 44, 46, 49]

word_search = ['-- Doubly Occupied --', ' -- Occupied --']
word_stop = '---------------------------------'
orbitals_symmetry = []
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
