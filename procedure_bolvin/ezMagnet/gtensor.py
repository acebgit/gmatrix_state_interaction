import numpy as np
import numpy.matlib
import json
import sys


# Constant, threshold and conversion factor definitions
G_E = -2.00231930436256
S_SQUARE_TOL = 0.05
INVCMTOAU= 1./219474.63068


# Load electronic state information from .json file_ms_notnull specified in first command
# line argument.
input_file_name = sys.argv[1] + ".json"
with open(input_file_name, 'r') as f:
    object_text = f.read()

input_dict = json.loads(object_text)

print_inputs = "print_inputs" in input_dict and input_dict["print_inputs"]
print_intermediates = ("print_intermediates" in input_dict) and input_dict[
    "print_intermediates"]

lowestEtargetState = "inputTargetState" not in input_dict
if not lowestEtargetState:
    inputTargetState = input_dict["inputTargetState"]

stateEnergiesDict = input_dict["stateEnergiesDict"]
stateTotalAngMomDict= input_dict["stateTotalAngMomDict"]
transAngMomListDict= input_dict["transAngMomListDict"]
aveTransSOCListDict = input_dict["aveTransSOCListDict"]


# Complex numbers in .json file_ms_notnull formatted as string; convert back to complex.
transAngMomVecDict = {key: np.array(
    [np.complex128(elem) for elem in transAngMomListDict[key]]) for key in transAngMomListDict}
aveTransSOCDict = {key: np.matrix(
    [[complex(elem) for elem in line] for line in aveTransSOCListDict[key]]) for key in aveTransSOCListDict}

if print_inputs:
    print(transAngMomListDict)
    print(aveTransSOCListDict)

if print_intermediates:
    for vector in transAngMomListDict:
        print("vector", vector)
    print(transAngMomVecDict)
    for matrix in aveTransSOCListDict.values():
        print("matrix", matrix)
    print("aveTransSOCDict", aveTransSOCDict)

# Determine label of target state either from .json file_ms_notnull or as label of lowest energy state
if lowestEtargetState:
    targetState = min(stateEnergiesDict, key=lambda x: stateEnergiesDict[x][0])
else:
    try:
        inputTargetState in stateEnergiesDict
    except:
        raise RuntimeError("Requested target state not found.")
    else:
       targetState = inputTargetState



# Collect all states_selected that differ in total angular momentum by less than
# S_SQUARE_TOL from targetState angular momentum.
totAng = stateTotalAngMomDict[targetState]
matchingStatesDict = {state: stateEnergiesDict[state] for state in filter(
    lambda x: abs(totAng - stateTotalAngMomDict[x]) < S_SQUARE_TOL, 
    stateTotalAngMomDict)}

# Sort states_selected matching in in angular momentum by energy.
stateList = sorted(matchingStatesDict, key=lambda x: stateEnergiesDict[x][0])
nStates = len(stateList)

if print_intermediates:
    #print("totalAngValues: " + str(totalAngValues))
    print("""Collecting states_selected of matching angular momentum and sorting them by 
          energy:""")
    print("totAng {ang} type(totAng) {angtype}".format(ang=totAng, angtype=type(totAng)))
    print("matchingStatesDict", matchingStatesDict)
    print("stateList", stateList)


def s_square(s):
    return s*(s+1.)

def spin_multiplity(angMom):
    '''Compute spin multiplicity that results in total spin within S_SQUARE_TOL of angMom'''
    maxssquare = float(np.ceil(2.*np.sqrt(angMom)))
    sspace = np.arange(0, maxssquare, 0.5)
    found_s_mult = False
    s_mult = -1.
    for s in sspace:
        if abs(s_square(s) - angMom) < S_SQUARE_TOL:
            s_mult = int(2*s+1)
            found_s_mult = True
            return s_mult
    if not found_s_mult:
        raise RuntimeError("Error: spin multiplicity cannot be determined; check for spin contamination")

# Find spin multiplicity corresponding to totAng and save in s_mult.
s_mult = spin_multiplity(totAng)

    
# Populate H0Mat, spinMat, angMomMat and SOCMat.
# inputs s_mult, nStates, stateEnergiesDict, transAngMomVecDict, aveTransSOCDict

H0Mat = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
spinMat_x = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
spinMat_y = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
spinMat_z = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
angMomMat_x = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
angMomMat_y = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
angMomMat_z = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
SOCMat = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)

# Calculate maximum spin projection for spin multiplicity of target state.
def m_max(s):
    return (s-1)/2
s = m_max(s_mult)

ms = np.linspace(-s, s, s_mult)

for j_col in range(nStates):
    for i_row in range(j_col+1):
        # spinMats are blockdiagonal with respect to symmetry
        if (i_row == j_col):
            # S_z operator contribution -> diagonal elements of spinMat_z
            for i in range(s_mult):
                H0Mat[i_row*s_mult+i,j_col*s_mult+i] = float(stateEnergiesDict[stateList[i_row]][0]) 
                spinMat_z[i_row*s_mult+i,j_col*s_mult+i] = -ms[i]
            # S_- operator contribution -> upper "subdiagonal" of spinMat_x/y
            for i in range(1, s_mult):
                spinMat_x[i_row*s_mult+i-1, j_col*s_mult+i] = -np.sqrt(s*(s+1) - ms[i]*(ms[i]-1))/2
                spinMat_y[i_row*s_mult+i-1, j_col*s_mult+i] = -np.sqrt(s*(s+1) - ms[i]*(ms[i]-1))*1.j/2
            # S_+ operator contribution -> lower "subdiagonal" of spinMat_x/y
            for i in range(s_mult-1):
                spinMat_x[i_row*s_mult+i+1, j_col*s_mult+i] = -np.sqrt(s*(s+1) - ms[i]*(ms[i]+1))/2
                spinMat_y[i_row*s_mult+i+1, j_col*s_mult+i] = -np.sqrt(s*(s+1) - ms[i]*(ms[i]+1))*-1.j/2
        else:
            # First symmetry label refers to Kets (cols and rows transposed).
            if (stateList[i_row]+"_"+stateList[j_col] in transAngMomVecDict):
                angMomVec = transAngMomVecDict[stateList[i_row]+"_"+stateList[j_col]]
                # Angular momenta are the diagonals of the off-diagonal blocks.
                for i in range(s_mult):
                    angMomMat_x[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[0]
                    angMomMat_y[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[1]
                    angMomMat_z[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[2]
                    angMomMat_x[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[0]*-1.
                    angMomMat_y[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[1]*-1.
                    angMomMat_z[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[2]*-1.
            if (stateList[j_col]+"_"+stateList[i_row] in transAngMomVecDict):
                angMomVec = transAngMomVecDict[stateList[j_col]+"_"+stateList[i_row]]
                for i in range(s_mult):
                    angMomMat_x[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[0]
                    angMomMat_y[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[1]
                    angMomMat_z[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[2]
                    angMomMat_x[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[0]*-1.
                    angMomMat_y[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[1]*-1.
                    angMomMat_z[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[2]*-1.
            if (stateList[i_row]+"_"+stateList[j_col] in aveTransSOCDict):
                socMatblock = aveTransSOCDict[stateList[i_row]+"_"+stateList[j_col]]
                SOCMat[j_col*s_mult:(j_col+1)*s_mult, i_row*s_mult:(i_row+1)*s_mult] = socMatblock
                SOCMat[i_row*s_mult:(i_row+1)*s_mult, j_col*s_mult:(j_col+1)*s_mult] = -1.*socMatblock
            if (stateList[j_col]+"_"+stateList[i_row] in aveTransSOCDict):
                socMatblock = aveTransSOCDict[stateList[j_col]+"_"+stateList[i_row]]
                SOCMat[i_row*s_mult:(i_row+1)*s_mult, j_col*s_mult:(j_col+1)*s_mult] = socMatblock
                SOCMat[j_col*s_mult:(j_col+1)*s_mult, i_row*s_mult:(i_row+1)*s_mult] = -1.*socMatblock

# Construct spin-orbit coupling matrix H.
H_soc_au = INVCMTOAU * SOCMat
H = H0Mat + H_soc_au

if print_intermediates:
    print("printing result of output processing:\n")
    print("H0Mat:\n"+str(H0Mat))
    print("spinMat_x:\n"+str(spinMat_x))
    print("spinMat_y:\n"+str(spinMat_y))
    print("spinMat_z:\n"+str(spinMat_z))
    print("angMomMat_x:\n"+str(angMomMat_x))
    print("angMomMat_y:\n"+str(angMomMat_y))
    print("angMomMat_z:\n"+str(angMomMat_z))
    print("SOCMat:\n"+str(SOCMat))
    print("H_soc_au:\n" + str(H_soc_au))

def sorted_diag(H):
    [e_vec_unord, c_mat_unord] = np.linalg.eigh(H)
    # Diagonalization may have changed order of states_selected. Reorder ensure that
    # state labels still correspond to LS-coupled state with largest coefficient
    # of that uncoupled state.
    col_max_coeff = np.where(abs(c_mat_unord) == numpy.amax(abs(c_mat_unord), axis=0)[0])
    reord = col_max_coeff[1].tolist()
    e_vec = e_vec_unord[reord]
    c_mat = c_mat_unord[:, reord]
    return e_vec, c_mat

# Compute spin orbit coupled energies and transformation matrix.
e_vec, c_mat = sorted_diag(H)
if print_intermediates:
    print("Spin coupled energy eigenvalues and coefficient of uncoupled states_selected:")
    print("e_vec", str(e_vec))
    print("c_mat", str(c_mat))

def trans_ls(c_mat, L_x_orig, L_y_orig, L_z_orig, S_x_orig, S_y_orig, S_z_orig):
    ''' Transform angular momentum and spin matrices into basis defined by c_mat.'''
    L_x_trans = c_mat.I @ L_x_orig @ c_mat
    L_y_trans = c_mat.I @ L_y_orig @ c_mat
    L_z_trans = c_mat.I @ L_z_orig @ c_mat
    S_x_trans = c_mat.I @ S_x_orig @ c_mat
    S_y_trans = c_mat.I @ S_y_orig @ c_mat
    S_z_trans = c_mat.I @ S_z_orig @ c_mat
    return L_x_trans, L_y_trans, L_z_trans, S_x_trans, S_y_trans, S_z_trans

# Transform spin and angular momentum matrices.
L_x_orig = np.matrix(angMomMat_x, dtype=np.cdouble)
L_y_orig = np.matrix(angMomMat_y, dtype=np.cdouble)
L_z_orig = np.matrix(angMomMat_z, dtype=np.cdouble)
S_x_orig = np.matrix(spinMat_x, dtype=np.cdouble)
S_y_orig = np.matrix(spinMat_y, dtype=np.cdouble)
S_z_orig = np.matrix(spinMat_z, dtype=np.cdouble)

L_x_trans, L_y_trans, L_z_trans, S_x_trans, S_y_trans, S_z_trans = trans_ls(
    c_mat, L_x_orig, L_y_orig, L_z_orig, S_x_orig, S_y_orig, S_z_orig)

if print_inputs:
    print("Untransformed angular momentum end spin matrices:")
    print("L_x_orig:\n" + str(L_x_orig))
    print("L_y_orig:\n" + str(L_y_orig))
    print("L_z_orig:\n" + str(L_z_orig))
    print("S_x_orig:\n" + str(S_x_orig))
    print("S_y_orig:\n" + str(S_y_orig))
    print("S_z_orig:\n" + str(S_z_orig))


if print_intermediates:
    print("Angular momentum end spin matrices in spin-orbit coupled basis:")
    print("L_x_trans:\n" + str(L_x_trans))
    print("L_y_trans:\n" + str(L_y_trans))
    print("L_z_trans:\n" + str(L_z_trans))
    print("S_x_trans:\n" + str(S_x_trans))
    print("S_y_trans:\n" + str(S_y_trans))
    print("S_z_trans:\n" + str(S_z_trans))


def calc_lambda(lx, ly, lz):
    Lambda = np.matlib.zeros((3, 3), dtype=np.cdouble)
    Lambda[0,0] = lx[1,0].real
    Lambda[1,0] = lx[1,0].imag
    Lambda[2,0] = lx[0,0]
    Lambda[0,1] = ly[1,0].real
    Lambda[1,1] = ly[1,0].imag
    Lambda[2,1] = ly[0,0]
    Lambda[0,2] = lz[1,0].real
    Lambda[1,2] = lz[1,0].imag
    Lambda[2,2] = lz[0,0]
    return Lambda

def calc_sigma(sx, sy, sz):
    Sigma = np.matlib.zeros((3, 3), dtype=np.cdouble)
    Sigma[0,0] = sx[1,0].real
    Sigma[1,0] = sx[1,0].imag
    Sigma[2,0] = sx[0,0]
    Sigma[0,1] = sy[1,0].real
    Sigma[1,1] = sy[1,0].imag
    Sigma[2,1] = sy[0,0]
    Sigma[0,2] = sz[1,0].real
    Sigma[1,2] = sz[1,0].imag
    Sigma[2,2] = sz[0,0]
    return Sigma

targetStateNo = stateList.index(targetState)
iup = targetStateNo*s_mult
idown = (targetStateNo+1)*s_mult
    
# Calculate Lambda, Sigma and g_combined from diagonal block of targetState
'''
        Compute G matrix for calculating g-tensor.
        G matrix is specified by block matrices of the sum of transformed
        angular momentum and spin matrices.
        The construction follows equations 12a-c and 13 in
        Bolvin, H. An alternative approach to the g-matrix: Theory and 
        applications ChemPhysChem 2006, 7, 1575–1589. DOI: 10.1002/cphc.200600051
'''
 
Lambda = calc_lambda(L_x_trans[iup:idown,iup:idown], 
                     L_y_trans[iup:idown,iup:idown], 
                     L_z_trans[iup:idown,iup:idown])
Sigma = calc_sigma(S_x_trans[iup:idown,iup:idown], 
                   S_y_trans[iup:idown,iup:idown], 
                   S_z_trans[iup:idown,iup:idown])

g_combined = 2. * Lambda + G_E * 2. * Sigma

if print_intermediates:
    print("Lambda\n" + str(Lambda))
    print("Sigma\n" + str(Sigma))
    print("g_combined\n" + str(g_combined))


# Diagonalize G_combined = g_combined.T @ g_combined
e_g, v_g = np.linalg.eigh(g_combined.T @ g_combined)
# Reorder in case of flip of cartesian coordinates during diagonalization.
v_g_max_coeff = np.where(abs(np.abs(v_g)) == numpy.amax(abs(v_g), axis=0))
g_reord = v_g_max_coeff[1].tolist()
e_g_reord = e_g[g_reord]

g_out_reord = np.sqrt(e_g_reord)
delta_g_reord = g_out_reord + G_E

print("g-tensor of {} state/ppt:".format(targetState), str(1000*delta_g_reord))


