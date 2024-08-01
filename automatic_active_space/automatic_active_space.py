#########################################################
# PROGRAM FOR ORGANIC MOLECULES
# TO AUTOMATICALLY GENERATE THE RASCI
# INPUT WITH AN AUTOMATIC ACTIVE SPACE AROUND HOMOs
#########################################################
# BE CAREFUL: it does not work if pseudopotentials are used in basis set, since
# number of electrons in SCF calc are not the real number of electrons
import sys

# 2 options:
# 1) Obtain QChem RASCI input with an active space selected:
# - automatically: 1 arguments needed (molecule_file)
# - manually: 4 arguments needed (molecule_file, active orbitals, active electrons, multiplicity)

# 2) If a file with improved active space obtained with "parser_excited_state" program is provided, 
# the active space will be the one in this file. 
# Leave "file_new_active_spaces" empty if not improved active space is required

molecule_file = str(sys.argv[1])
file_new_active_spaces = "activespace_gtensor50.txt"


def take_data(filee):
    """
    Program to take the number of orbitals and electrons in the molecule file.
    :param: filee
    :return: electrons, multiplicity
    """
    def remove_numbers(input_string):
        return ''.join([char for char in input_string if not char.isdigit()])
    
    elements = {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Ga": 31,
        "Ge": 32,
        "As": 33,
        "Se": 34,
        "Br": 35,
        "Kr": 36,
        "Rb": 37,
        "Sr": 38,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "In": 49,
        "Sn": 50,
        "Sb": 51,
        "Te": 52,
        "I": 53
    }

    electrons = 0

    # Read "molecule" or "xyz" file 
    with open(filee, 'r') as file:
        lines = file.readlines()
    
    for line in range(0, len(lines)):

        if "molecule" in lines[line]: 
            charge = int(lines[line+1].split()[0])
            multiplicity = int(lines[line+1].split()[1])

        if not lines[line].isspace(): # Only continue if there is not a blanck line
            # By the time it is written in g-matrix, delete the numbers from the element name 
            element = remove_numbers(lines[line].split()[0])
            
            if element in elements:
                electrons += elements[element]
    
    electrons = electrons - charge
    return electrons, multiplicity


def take_ras_features(multiplicity):
    """
    Take active orbitals, active electrons and multiplicity manually or automatically. 
    :param: input_multiplicity
    :return: ras_act, ras_elec, selected_multiplicity
    """
    if len(sys.argv) == 2:
        select_multiplicity = multiplicity

        if multiplicity == 1 or multiplicity == 2:
            # With singlets or doublets, minimum active space is [2,2]
            ras__elec = 2 
        else: 
            # With the rest of multiplicities, minimum active space is [multip-1, multip-1]
            ras__elec = select_multiplicity - 1
        ras__act = ras__elec

    elif len(sys.argv) == 5:
        ras__elec = int(sys.argv[2])
        ras__act = int(sys.argv[3])
        select_multiplicity = int(sys.argv[4])
    
    else:
        raise ValueError("No manual or automatic selection has been done.")
    return ras__act, ras__elec, select_multiplicity


def take_alpha_beta_orbitals(select_multiplicity, ras__elec):
    """"
    Take the number of unpaired electrons, that with the electrons in RAS define the alpha and beta orbitals
    :param: selected_multiplicity, ras_elec
    :return: unpaired_electrons, ras_elec_alpha, ras_elec_beta
    """
    unpaired__electrons = select_multiplicity-1
    ras__elec_alpha = ras__elec//2 + ras__elec % 2 + unpaired__electrons//2
    ras__elec_beta = ras__elec - ras__elec_alpha

    if (ras__elec_alpha - ras__elec_beta) != unpaired__electrons:
        raise ValueError("Active space electrons and multiplicity are not compatible")
    elif ras__elec_beta < 0:
        raise ValueError("ras_elec_beta cannot be negative. Choose another ras_multiplicity.")
    return unpaired__electrons, ras__elec_alpha, ras__elec_beta


def take_domos_vomos_somos(total_electrns, unpaired_elects):
    ndomos = total_electrns // 2 - unpaired_elects // 2 
    domos = [i for i in range(1, ndomos+1)]
    somos = [i for i in range(len(domos)+1, len(domos)+unpaired_elects+1)]
    vomos = [i for i in range(len(domos)+unpaired_elects+1, len(domos)+unpaired_elects+20)]
    return domos, somos, vomos


def construct_active_space(total__electrons, unpaired__electrons, ras__elec, ras__act, domoss, somoss, vomoss):
    """
    Take the number of unpaired electrons, that with the electrons in RAS define the alpha and beta orbitals
    :param: selected_multiplicity, ras_elec
    :return: unpaired_electrons, ras_elec_alpha, ras_elec_beta
    """
    ras_actorb = []
    count_ras_elec = ras__elec
    count_ras_act = ras__act

    for orbital in somoss:
        if count_ras_elec > 0 and count_ras_act > 0:
            ras_actorb.append(orbital)
            count_ras_elec-=1
            count_ras_act-=1

    for i in range(len(domoss)-1, 1, -1):
        if count_ras_elec > 0 and count_ras_act > 0:
            ras_actorb.append(domoss[i])
            vaccant_index = abs(i - len(domoss)+1)
            ras_actorb.append(vomoss[vaccant_index])
            count_ras_elec-=2
            count_ras_act-=1

    ras_actorb = sorted(ras_actorb)
    return ras_actorb


def take_new_active_spaces(file__new_active_spaces):
    with open(file__new_active_spaces, 'r') as file:
        lines = file.readlines()

    active_dict = {}
    for line in range(0, len(lines)):
        if "out" in lines[line]: 
            output_file = lines[line].split("_")[0]+".molecule"
        
        if "Final active space" in lines[line]: 
            string_list = str(lines[line].split(";")[1].replace("[","").replace("]","").replace(",",""))
            actual_list = [int(i) for i in string_list.split()]
            active_dict[output_file] = actual_list
    return active_dict
    

# 1) Take the total number of electrons and orbitals
total_electrons, input_multiplicity = take_data(molecule_file)

ras_act, ras_elec, selected_multiplicity = take_ras_features(input_multiplicity)

# 2) Take the alpha and beta electrons
unpaired_electrons, ras_elec_alpha, ras_elec_beta = take_alpha_beta_orbitals(selected_multiplicity, ras_elec)

domos, somos, vomos = take_domos_vomos_somos(total_electrons, unpaired_electrons)

# 3) Form the active space 
if file_new_active_spaces: 
    # If we have file with new active spaces, redo the space
    ras_act_orb = take_new_active_spaces(file_new_active_spaces)[molecule_file]

    ras_act = len(ras_act_orb)
    
    ras_elec = sum(2 for orbit in ras_act_orb if orbit in domos)
    ras_elec += sum(1 for orbit in ras_act_orb if orbit in somos)

    unpaired_electrons, ras_elec_alpha, ras_elec_beta = take_alpha_beta_orbitals(selected_multiplicity, ras_elec)

else:
    ras_act_orb = construct_active_space(total_electrons, unpaired_electrons, ras_elec, ras_act, domos, somos, vomos)


qchem_rasci_input = '''$comment
 RAS-SF CALCULATION
$end

$molecule
 READ ARCHIVO
$end

$rem
!**CALCULATION TYPE***     
JOBTYPE            sp    
EXCHANGE           HF
CORRELATION        RASCI    !RASCI=Casanova; RASCI2=Zimmerman
BASIS              def2-TZVP
!BASIS2             STO-3G

!***SCF REFERENCE***
!SKIP_SCFMAN        FALSE
UNRESTRICTED       FALSE
!SCF_GUESS          core
!SCF_GUESS_MIX      FALSE
!SCF_GUESS_ALWAYS   FALSE
!MAX_SCF_CYCLES     500
!SCF_CONVERGENCE    10
!SCF_ALGORITHM      DIIS !DIIS,DM,DIIS_DM,GDM,DIIS_GDM,GD

!SYMMETRY           TRUE   !use of symmetry for calculating integrals
!SYM_IGNORE         FALSE  !turn off using symmetry
!SYM_TOL            5      !tolerance for determining point group
!USE_ABELIAN_SUBGROUP      !define to use abelian group

!THRESH             14
!SCF_FINAL_PRINT     1 !0,1,2,3

!****RAS-SF****
RAS_ROOTS          30       !number of RAS-CI roots to be computed
SET_ITER           100     !number of iterations in RASCI
RAS_ELEC           2       !number of electrons in RAS2
RAS_ELEC_ALPHA            !spin-α electrons in RAS2 
RAS_ELEC_BETA             !spin-β electrons in RAS2

RAS_ACT            2       !number of orbitals in RAS2 
RAS_ACT_ORB       1       !user-selected RAS2 orbitals
RAS_OCC            2       !number of orbitals in RAS1
N_FROZEN_CORE      FC     !Number of frozen core orbitals
!N_FROZEN_VIRTUAL   0      !Number of frozen virtual orbitals

!RAS_DO_HOLE        TRUE   !presence of hole excitations in RAS-CI wave function
!RAS_DO_PART        TRUE   !presence of particle excitations
!RAS_AMPL_PRINT     10     !threshold (×10 2 ) for the CI amplitudes

!RAS_GUESS_CS       0      !number of closed shell guess conﬁgurations
!RAS_SPIN_MULT      0      !spin multiplicity. Only for Ms = 0 
!RAS_SRDFT          FALSE  !short-range density functional RAS-CI
!RAS_NFRAG          0     !fragments of excitation analysis in RASCI

RAS_PRINT          2      !level of information
!RAS_NATORB         TRUE   !computation of the natural orbital occupancies
!RAS_NATORB_STATE  0      !Saves NO of i RAS computed state in .fchk ﬁle
!RAS_FOD            TRUE

!****RAS-SOC****
CALC_SOC           TRUE   !calculate SOC for EOM-CC, RAS-CI, ADC, TDDFT/TDA and TDDFT
!RAS_SOC_2E        TRUE   !two-electron mean-ﬁeld contribution, default=true
RAS_SOC_SYM_DENS   TRUE   !averaging of α and β densities

!****MEMORY****
MEM_TOTAL          50000 !Default 2000 MB
MEM_STATIC         2000  !Default 64 MB

!**Molecular Properties and Analysis **                                                                                          
GUI                0       !Checkpoints generation
STATE_ANALYSIS     TRUE
$end
'''
lines_list = qchem_rasci_input.splitlines()


rasci_output = molecule_file.split('.')[0] + "_" + str(ras_elec) + "_" + str(ras_act) + ".inp"
print("Electrons:", total_electrons, ", Multiplicity:", selected_multiplicity, ", Active space: ", ras_act_orb)

with open(rasci_output, "w", encoding="utf8") as f:
    for line in lines_list:
        if 'READ' in line:
            f.write('READ ' + str(molecule_file) + '\n')
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
            f.write('RAS_OCC ' + str(ras_act_orb[0] - 1) + '\n')
        else:
            f.write(line + '\n')
