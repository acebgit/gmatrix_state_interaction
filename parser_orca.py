__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np

def get_gtensor(file, ppm):
    """
    Obtain the total number of states selected in RAS-CI in Q-Chem output.
    :param: file
    :return: totalstates
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    word_search = [' Delta-g ']
    g_shifts = []
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            elements = [float(line[1]), float(line[2]), float(line[3])]

            if ppm == 1:
                [g_shifts.append(np.round(i * 10**6, 2)) for i in elements]
            else:
                [g_shifts.append(np.round(i * 10**3, 2)) for i in elements]

            break
    return g_shifts

filee = 'orca_outs/quinoline_orca_dft.out'
gshifts = get_gtensor(filee, ppm = 0)
print(filee, ':', *gshifts)
