from scipy import constants
import sys 

file_path = str(sys.argv[1])
output_path = file_path #.replace(".out", "_2.out")

with open(file_path, 'r') as file:
    lines_before = file.readlines()

word_search = 'RAS-CI total energy for state'
total_energy = []
excit_energy = []
lines_after = []
line_to_jump = -1

for line in range(0, len(lines_before)):
    if "RAS-CI total energy for state" in lines_before[line]:
        lines_after.append(lines_before[line])

        total_energy.append(float(lines_before[line].split(":")[1]))
        if "*******" in lines_before[line+1].split("=")[1]:
            ener = (float(lines_before[line].split(":")[1]) - total_energy[0]) * constants.physical_constants['Hartree energy in eV'][0]
            lines_before[line+1] = lines_before[line+1].replace("***************", str(ener))
            lines_after.append(lines_before[line+1])
            line_to_jump = line+1
    
    elif line == line_to_jump:
        pass
    
    else:
        lines_after.append(lines_before[line])

lines_after_str = "".join(lines_after)
with open(output_path, 'w') as file:
    file.write(lines_after_str)
