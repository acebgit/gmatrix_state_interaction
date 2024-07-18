import sys
import numpy as np

options = 1 # 0) take g-shift and times; 1) take g-matrix
ppm = 1
file_path = str(sys.argv[1])

with open(file_path, 'r') as file:
    lines = file.readlines()  

if options == 0:
    oriented_list = []
    for i in range(len(lines)):
        if 'Delta-g' in lines[i]:
                gshift = [float(lines[i].split()[1]),
                float(lines[i].split()[2]),
                float(lines[i].split()[3])]

        if 'Orientation:' in lines[i]:
            for j in range(i+1, i+4):
                oriented_list.append([abs(float(lines[j].split()[1])),
                abs(float(lines[j].split()[2])),
                abs(float(lines[j].split()[3]))])
            # .strip() function delete whitespaces at beginning and end
        
        if 'TOTAL RUN TIME' in lines[i]:
            tiempo = lines[i].replace("TOTAL RUN TIME:","")
            minutes = float(tiempo.split()[0])*24*60+float(tiempo.split()[2])*60+float(tiempo.split()[4])+float(tiempo.split()[6])/60+float(tiempo.split()[8])/(60*1000)

    max_indices = [sublist.index(max(sublist)) for sublist in oriented_list]
    gshift_oriented = [gshift[max_indices[i]] for i in range(0,3)]

    if ppm == 0:
        gshift_oriented = [np.round((i * 10**3),3) for i in gshift_oriented]
    elif ppm == 1:
        gshift_oriented = [np.round((i * 10**6),3) for i in gshift_oriented]

    print("g-shift oriented:", gshift_oriented)
    print("Time:", np.round(minutes, 3))

elif options == 1: 
    for i in range(len(lines)):
        if 'The g-matrix:' in lines[i]:
            gmatrix = []
            for j in range(0, 3):
                print(np.round(float(lines[i+1+j].strip().split()[0]), decimals=7), 
                np.round(float(lines[i+1+j].strip().split()[1]), decimals=7), 
                np.round(float(lines[i+1+j].strip().split()[2]), decimals=7))
            print()
