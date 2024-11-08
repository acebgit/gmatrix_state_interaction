#################################################
#    PROGRAM TO PASS FROM MOLECULE FILE TO
#    XYZ FILE, INCLUDING CHANGES IN UNITS FROM
#    ANGSTROM TO BOHR
#################################################
import sys

zmat = 2 
# 0) read and write coordinates 
# 1) from z-matrix (.molecule) to coordinates (.xyz)
# 2) from z-matrix (.molecule) to gzmtfile (.gzmt) 

# # Get the Bohr radius (in meters)
# bohr_radius = physical_constants['Bohr radius'][0]
# # Conversion factor from angstroms to bohrs
# conversion_factor = angstrom / bohr_radius

file = str(sys.argv[1])
with open(file, 'r') as f:
    filelist = f.readlines()

output = []
nmolecules = 0
molecules_map = {}
molecule = file.replace('.molecule', '')

if zmat == 0:
    with open(molecule + ".xyz", 'w') as filee:
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
    with open(molecule + ".xyz", 'w') as filee:
        for line in filelist:
            linesplit = line.split('!')[0].split()
            print(linesplit)

            if len(linesplit) == 1 and "$molecule" not in linesplit and "$end" not in linesplit: # Text lines
                filee.write(f"{linesplit[0]}\n")

            if len(linesplit) == 3: # Lines with distances
                newdist = float(linesplit[2]) 
                filee.write(f"{linesplit[0]:<3}{linesplit[1]:<3}{newdist:<10.6f}\n")
            
            elif len(linesplit) == 5: # Lines with distances
                newdist = float(linesplit[2]) 
                filee.write(f"{linesplit[0]:<3}{linesplit[1]:<3}{newdist:<10.6f}{linesplit[3]:<3}{linesplit[4]:<13}\n")

            elif len(linesplit) == 8: # Lines with distances
                newdist = float(linesplit[2]) 
                filee.write(f"{linesplit[0]:<3}{linesplit[1]:<3}{newdist:<10.6f}{linesplit[3]:<3}{linesplit[4]:<13}{linesplit[5]:<3}{linesplit[6]:<13}\n")

elif zmat == 2:
    with open(molecule + ".gzmt", 'w') as filee:
        for line in filelist:
            linesplit = line.split('!')[0].split()

            # Pass line or take the first atom symbol
            if len(linesplit) == 1 and "$molecule" not in linesplit and "$end" not in linesplit: # Text lines
                # filee.write(f"{linesplit[0][0]:<3}{'0':<3}{'0':<3}{'0':<3}{0:<10.6f}{0:<10.6f}{0:<10.6f}\n")
                filee.write(f"{linesplit[0][0]:<3}\n")
                nmolecules += 1
                molecules_map.update({nmolecules: linesplit[0]})   

            elif len(linesplit) == 3: # Lines with the charge and multiplicity
                newdist = float(linesplit[2])
                molecule1 = next((k for k, v in molecules_map.items() if v == linesplit[1]), None) 

                filee.write(f"{linesplit[0][0]:<3}{molecule1:<3}{newdist:<10.5f}\n")
                nmolecules += 1
                molecules_map.update({nmolecules: linesplit[0]}) 

            elif len(linesplit) == 5:
                newdist = float(linesplit[2]) 
                newang = float(linesplit[4])
                molecule1 = next((k for k, v in molecules_map.items() if v == linesplit[1]), None) 
                molecule2 = next((k for k, v in molecules_map.items() if v == linesplit[3]), None)
                
                filee.write(f"{linesplit[0][0]:<3}{molecule1:<3}{newdist:<10.5f}{molecule2:<3}{newang:<10.2f}\n")
                nmolecules += 1
                molecules_map.update({nmolecules: linesplit[0]})
                
            elif len(linesplit) == 7:
                newdist = float(linesplit[2]) 
                newang = float(linesplit[4])
                newdied = float(linesplit[6])
                molecule1 = next((k for k, v in molecules_map.items() if v == linesplit[1]), None) 
                molecule2 = next((k for k, v in molecules_map.items() if v == linesplit[3]), None)
                molecule3 = next((k for k, v in molecules_map.items() if v == linesplit[5]), None)

                filee.write(f"{linesplit[0][0]:<3}{molecule1:<3}{newdist:<10.5f}{molecule2:<3}{newang:<10.2f}{molecule3:<3}{newdied:<10.2f}\n")
                nmolecules += 1
                molecules_map.update({nmolecules: linesplit[0]})
