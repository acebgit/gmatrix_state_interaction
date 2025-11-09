#########################################################
# PROGRAM FOR ORGANIC MOLECULES
# TO AUTOMATICALLY GENERATE THE RASCI
# INPUT WITH AN AUTOMATIC ACTIVE SPACE AROUND HOMO.
#########################################################
import sys

geometry_file = "C44H20B4.molecule"  # str(sys.argv[1])
ras_act = 8
ras_elec = 8
ras_multiplicity = 1

rasci_input = 'qchem_rasci.inp'
rasci_output = geometry_file.split('.')[0] + "_" + str(ras_elec) + "_" + str(ras_act) + ".inp"


def take_data(filee):
    electrons_dict = {
        'H': 1,
        'He': 2,
        'Li': 3,
        'Be': 4,
        'B': 5,
        'C': 6,
        'N': 7,
        'O': 8,
        'F': 9,
        'Ne': 10,
    }
    electrons = 0

    with open(filee, encoding="utf8") as ff:
        for linee in ff:
            element = linee.split()[0]
            if element in electrons_dict and linee != []:
                electrons += electrons_dict[element]
    orbitals = int(electrons / 2) + electrons % 2
    return electrons, orbitals


total_electrons, total_orbitals = take_data(geometry_file)

ras_elec_alpha = ras_elec//2 + ras_elec % 2 + (ras_multiplicity-1)//2
ras_elec_beta = ras_elec - ras_elec_alpha
if ras_elec_beta < 0:
    raise ValueError("ras_elec_beta cannot be negative. Choose another ras_multiplicity.")

total_spin = (ras_multiplicity-1)//2
homo = total_electrons//2 + total_electrons % 2 + total_spin
ras_act_orb = [i for i in range(homo - ras_elec//2, homo+ras_elec//2 + ras_elec % 2)]
ras_occ = ras_act_orb[0] - 1
print(homo, total_spin)
print(ras_act_orb)
exit()

with open(rasci_input, "r", encoding="utf8") as f:
    lines = f.readlines()

with open(rasci_output, "w", encoding="utf8") as f:
    for line in lines:
        if 'READ' in line:
            f.write('READ ' + str(geometry_file) + '\n')
        elif 'RAS_ELEC ' in line:
            f.write('RAS_ELEC ' + str(ras_elec) + '\n')
        elif 'RAS_ELEC_ALPHA ' in line:
            f.write('RAS_ELEC_ALPHA ' + str(ras_elec_alpha) + '\n')
        elif 'RAS_ELEC_BETA ' in line:
            f.write('RAS_ELEC_BETA ' + str(ras_elec_beta) + '\n')
        elif 'RAS_ACT ' in line:
            f.write('RAS_ACT ' + str(ras_act) + '\n')
        elif 'RAS_ACT_ORB ' in line:
            f.write('RAS_ACT_ORB ' + str(ras_act_orb) + '\n')
        elif 'RAS_OCC ' in line:
            f.write('RAS_OCC ' + str(ras_occ) + '\n')
        else:
            f.write(line)
