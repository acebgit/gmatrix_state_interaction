import sys
import numpy as np
import pandas as pd

options = 1 # 0) take g-shift and times; 1) take g-shift, g-matrix and orientation
ppm = 1
file_path = str(sys.argv[1])

with open(file_path, 'r') as file:
    lines = file.readlines()  

for i in range(len(lines)):
    if 'The g-matrix:' in lines[i]:
        gmatrix = []
        for j in range(0, 3):
            gmatrix.append([float(lines[i+1+j].strip().split()[0]), 
                            float(lines[i+1+j].strip().split()[1]), 
                            float(lines[i+1+j].strip().split()[2])])

    if 'Delta-g' in lines[i]:
            gshift = [float(lines[i].split()[1]),float(lines[i].split()[2]),float(lines[i].split()[3])]

    if 'Orientation:' in lines[i]:
        oriented_list = []
        for j in range(0, 3):
            oriented_list.append([float(lines[i+1+j].strip().split()[1]), 
                                  float(lines[i+1+j].strip().split()[2]),
                                  float(lines[i+1+j].strip().split()[3])])

    if 'TOTAL RUN TIME' in lines[i]:
        tiempo = lines[i].replace("TOTAL RUN TIME:","")
        minutes = float(tiempo.split()[0])*24*60+float(tiempo.split()[2])*60+float(tiempo.split()[4])+float(tiempo.split()[6])/60+float(tiempo.split()[8])/(60*1000)

# Reorder the g-values taking into account the orientation
max_indices = [sublist.index(max(sublist, key=abs)) for sublist in oriented_list]
gshift_oriented = [gshift[max_indices[i]] for i in range(0,3)]

if ppm == 0:
    gshift_oriented = [round((i * 10**3),3) for i in gshift_oriented]
elif ppm == 1:
    gshift_oriented = [round((i * 10**6),3) for i in gshift_oriented]

if options == 0:
    print("g-shift oriented:", gshift_oriented[0], gshift_oriented[1], gshift_oriented[2])
    # print("Time:", np.round(minutes, 3))

elif options == 1:
    print("g-shift oriented:", gshift_oriented[0], gshift_oriented[1], gshift_oriented[2])
    print()
    print('The g-matrix:')
    df = pd.DataFrame(gmatrix)
    print(df.to_string(index=False, header=False))
    print()
    print('Orientation:')
    df = pd.DataFrame(oriented_list)
    print(df.to_string(index=False, header=False))
