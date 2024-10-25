#################################################
#    PROGRAM TO PASS FROM MOLECULE FILE TO
#    XYZ FILE, INCLUDING CHANGES IN UNITS FROM
#    ANGSTROM TO BOHR
#################################################
# from scipy.constants import angstrom, physical_constants
import numpy as np
import sys

zmat = 0 
# 0) read and write coordinates 1) read z-matrix and write coordinates 2)read and write z-matrix

# # Get the Bohr radius (in meters)
# bohr_radius = physical_constants['Bohr radius'][0]
# # Conversion factor from angstroms to bohrs
# conversion_factor = angstrom / bohr_radius

file = str(sys.argv[1])
with open(file, 'r') as f:
    filelist = f.readlines()

output = []

molecule = file.replace('.molecule', '')
with open(molecule + ".xyz", 'a+') as filee:

    if zmat == 0:
        natoms = 0
        valid_lines = []

        for line in filelist:
            linesplit = line.split('!')[0].split()
            if len(linesplit) == 4:
                valid_lines.append(f"{' '.join(linesplit)}\n")
                natoms += 1
        
        filee.write(f"{natoms}\n\n")  # Write atom count and an empty line
        filee.write(''.join(valid_lines) + '\n')  # Write all valid lines

    if zmat == 1:
        for line in filelist:
            linesplit = line.split('!')[0].split()
            print(linesplit)

            if len(linesplit) == 1 and "$molecule" not in linesplit and "$end" not in linesplit: # Text lines
                filee.write(f"{linesplit[0]}\n")

            if len(linesplit) == 3: # Lines with distances
                newdist = float(linesplit[2]) 
                filee.write(f"{linesplit[0]:<5}{linesplit[1]:<5}{newdist:<13.6f}\n")
            
            elif len(linesplit) == 5: # Lines with distances
                newdist = float(linesplit[2]) 
                filee.write(f"{linesplit[0]:<5}{linesplit[1]:<5}{newdist:<13.6f}{linesplit[3]:<5}{linesplit[4]:<13}\n")

            elif len(linesplit) == 8: # Lines with distances
                newdist = float(linesplit[2]) 
                filee.write(f"{linesplit[0]:<5}{linesplit[1]:<5}{newdist:<13.6f}{linesplit[3]:<5}{linesplit[4]:<13}{linesplit[5]:<5}{linesplit[6]:<13}\n")

    elif zmat == 2:
        for line in filelist:
            linesplit = line.split('!')[0].split()

            if len(linesplit) == 1 and "$molecule" not in linesplit and "$end" not in linesplit: # Text lines
                filee.write(f"{linesplit[0]:<5}{'0':<5}{'0':<5}{'0':<5}{0:<13.6f}{0:<13.6f}{0:<13.6f}\n")

            elif len(linesplit) == 3: # Lines with the charge and multiplicity
                newdist = float(linesplit[2]) 
                filee.write(f"{linesplit[0]:<5}{linesplit[1]:<5}{'0':<5}{'0':<5}{newdist:<13.6f}{0:<13.6f}{0:<13.6f}\n")

            elif len(linesplit) == 5:
                newdist = float(linesplit[2]) 
                newang = float(linesplit[4])
                filee.write(f"{linesplit[0]:<5}{linesplit[1]:<5}{linesplit[3]:<5}{'0':<5}{newdist:<13.6f}{newang:<13.6f}{0:<13.6f}\n")
                
            elif len(linesplit) == 8:
                newdist = float(linesplit[2]) 
                newang = float(linesplit[4])
                newdied = float(linesplit[6])
                filee.write(f"{linesplit[0]:<5}{linesplit[1]:<5}{linesplit[3]:<5}{linesplit[5]:<5}{newdist:<13.6f}{newang:<13.6f}{newdied:<13.6f}\n")

        # f.write(f"{len(output)}\n")
        # f.write(f"{molecule}\n")
        # for item in output:
        #     f.write(f"{item}\n")