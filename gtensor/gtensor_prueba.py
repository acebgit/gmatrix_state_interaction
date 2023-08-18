__author__ = 'Sven Kähler'

import numpy as np
import numpy.matlib
import json
import sys


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


def m_max(s):
    return (s - 1) / 2


def sorted_diag(H):
    [e_vec_unord, c_mat_unord] = np.linalg.eig(H)
    col_max_coeff = np.where(abs(c_mat_unord) == numpy.amax(abs(c_mat_unord), axis=0)[0])
    reord = col_max_coeff[1].tolist()
    e_vec = e_vec_unord[reord]
    c_mat = c_mat_unord[:, reord]
    return e_vec, c_mat


def trans_ls(c_mat, L_x_orig, L_y_orig, L_z_orig, S_x_orig, S_y_orig, S_z_orig):
    ''' Transform angular momentum and spin matrices into basis defined by c_mat.'''
    L_x_trans = c_mat.I @ L_x_orig @ c_mat
    L_y_trans = c_mat.I @ L_y_orig @ c_mat
    L_z_trans = c_mat.I @ L_z_orig @ c_mat
    S_x_trans = 2.0*c_mat.I @ S_x_orig @ c_mat
    S_y_trans = 2.0*c_mat.I @ S_y_orig @ c_mat
    S_z_trans = 2.0*c_mat.I @ S_z_orig @ c_mat
    return L_x_trans, L_y_trans, L_z_trans, S_x_trans, S_y_trans, S_z_trans


def compute_state_g_tensor(gsubblock_x, gsubblock_y, gsubblock_z):
    ''' Compute g-tensor from G matrix
        G matrix is specified by block matrices of the sum of transformed
        angular momentum and spin matrices.
        The construction follows equation 13 in
        Bolvin, H. An alternative approach to the g-matrix: Theory and
        applications ChemPhysChem 2006, 7, 1575–1589.
    '''
    g_mat = np.zeros((3,3), dtype=np.cdouble)
    g_mat[0,0] = 2. * np.sum(np.multiply(gsubblock_x, gsubblock_x.T))
    g_mat[0,1] = 2. * np.sum(np.multiply(gsubblock_x, gsubblock_y.T))
    g_mat[0,2] = 2. * np.sum(np.multiply(gsubblock_x, gsubblock_z.T))
    g_mat[1,0] = 2. * np.sum(np.multiply(gsubblock_y, gsubblock_x.T))
    g_mat[1,1] = 2. * np.sum(np.multiply(gsubblock_y, gsubblock_y.T))
    g_mat[1,2] = 2. * np.sum(np.multiply(gsubblock_y, gsubblock_z.T))
    g_mat[2,0] = 2. * np.sum(np.multiply(gsubblock_z, gsubblock_x.T))
    g_mat[2,1] = 2. * np.sum(np.multiply(gsubblock_z, gsubblock_y.T))
    g_mat[2,2] = 2. * np.sum(np.multiply(gsubblock_z, gsubblock_z.T))
    [g_eig_unord, g_vec_unord] = np.linalg.eigh(g_mat)
    g_vec_max_coeff     = np.where(abs(np.abs(g_vec_unord)) == numpy.amax(abs(g_vec_unord), axis=0))
    g_reord = g_vec_max_coeff[1].tolist()
    g_eig_reord = g_eig_unord[g_reord]
    return np.sqrt(g_eig_reord)

inputTargetState = ""
S_SQUARE_TOL = 0.05
INVCMTOAU= 1./219474.63068

# 1) TAKING THE INPUT FORMAT FILE
# input_file_name = sys.argv[1] + ".json"
input_file_name = "fe_pyms2_def2tzvp_ras.out.json"
# input_file_name = "ag_cucl4_2-_D2h_eomip_soc_def2-TZVP.in.out.json"
with open(input_file_name, 'r') as f:
    object_text = f.read()
input_dict = json.loads(object_text)

# 2) TAKE THE LOWEST ENERGY TRANSITION
lowestEtargetState = "inputTargetState" not in input_dict
if not lowestEtargetState:
    inputTargetState = input_dict["inputTargetState"]

StateEnergiesDict = input_dict["StateEnergiesDict"]
StateTotalAngMomDict= input_dict["StateTotalAngMomDict"]
TransAngMomListDict= input_dict["TransAngMomListDict"]
AveTransSOCListDict = input_dict["AveTransSOCListDict"]

eomTransAngMomVecDict = {key: np.array(
    [np.complex128(elem) for elem in TransAngMomListDict[key]]) for key in TransAngMomListDict}
eomAveTransSOCDict = {key: np.matrix(
    [[complex(elem) for elem in line] for line in AveTransSOCListDict[key]]) for key in AveTransSOCListDict}

if lowestEtargetState:
    targetState = min(StateEnergiesDict, key=lambda x: StateEnergiesDict[x][0])
else:
    try:
        inputTargetState in StateEnergiesDict
    except:
        raise RuntimeError("Requested target state not found.")
    else:
       targetState = inputTargetState

# 3) Collect all state that have the same total angular momentum as totAng in matchingStatesDict
totalAngValues = [s for s in list(StateTotalAngMomDict.values())]
totAng = StateTotalAngMomDict[targetState]
matchingStatesDict = {state: StateEnergiesDict[state] for state in filter(
    lambda x: abs(totAng - StateTotalAngMomDict[x]) < S_SQUARE_TOL, 
    StateTotalAngMomDict)}

stateList = sorted(matchingStatesDict, key=lambda x: StateEnergiesDict[x][0])
nStates = len(stateList)

# Find spin multiplicity corresponding to totAng and save in s_mult.
s_mult = spin_multiplity(totAng)

# Determine index of target state and save as targetStateNo.
targetStateNo = stateList.index(targetState)

s = m_max(s_mult)
ms = np.linspace(-s, s, s_mult)

H0Mat = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
spinMat_x = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
spinMat_y = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
spinMat_z = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
angMomMat_x = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
angMomMat_y = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
angMomMat_z = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)
SOCMat = np.matlib.zeros((nStates*s_mult, nStates*s_mult), dtype=np.cdouble)

for j_col in range(nStates):
    for i_row in range(j_col+1):
        # spinMats are blockdiagonal with respect to symmetry
        if (i_row == j_col):
            # S_z operator contribution -> diagonal elements of spinMat_z
            for i in range(s_mult):
                H0Mat[i_row*s_mult+i,j_col*s_mult+i] = float(StateEnergiesDict[stateList[i_row]][0]) 
                spinMat_z[i_row*s_mult+i,j_col*s_mult+i] = ms[i]
            # S_- operator contribution -> upper "subdiagonal" of spinMat_x/y
            for i in range(1, s_mult):
                spinMat_x[i_row*s_mult+i-1, j_col*s_mult+i] = np.sqrt(s*(s+1) - ms[i]*(ms[i]-1))/2
                spinMat_y[i_row*s_mult+i-1, j_col*s_mult+i] = np.sqrt(s*(s+1) - ms[i]*(ms[i]-1))*1.j/2
            # S_+ operator contribution -> lower "subdiagonal" of spinMat_x/y
            for i in range(s_mult-1):
                spinMat_x[i_row*s_mult+i+1, j_col*s_mult+i] = np.sqrt(s*(s+1) - ms[i]*(ms[i]+1))/2
                spinMat_y[i_row*s_mult+i+1, j_col*s_mult+i] = np.sqrt(s*(s+1) - ms[i]*(ms[i]+1))*-1.j/2
        else:
            # First symmetry label refers to Kets (cols and rows transposed).
            if (stateList[i_row]+"_"+stateList[j_col] in eomTransAngMomVecDict):
                angMomVec = eomTransAngMomVecDict[stateList[i_row]+"_"+stateList[j_col]]
                # Angular momenta are the diagonals of the off-diagonal blocks.
                for i in range(s_mult):
                    angMomMat_x[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[0]
                    angMomMat_y[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[1]
                    angMomMat_z[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[2]
                    angMomMat_x[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[0]*-1.
                    angMomMat_y[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[1]*-1.
                    angMomMat_z[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[2]*-1.
            if (stateList[j_col]+"_"+stateList[i_row] in eomTransAngMomVecDict):
                angMomVec = eomTransAngMomVecDict[stateList[j_col]+"_"+stateList[i_row]]
                for i in range(s_mult):
                    angMomMat_x[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[0]
                    angMomMat_y[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[1]
                    angMomMat_z[i_row*s_mult+i, j_col*s_mult+i] = angMomVec[2]
                    angMomMat_x[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[0]*-1.
                    angMomMat_y[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[1]*-1.
                    angMomMat_z[j_col*s_mult+i, i_row*s_mult+i] = angMomVec[2]*-1.
            if (stateList[i_row]+"_"+stateList[j_col] in eomAveTransSOCDict):
                socMatblock = eomAveTransSOCDict[stateList[i_row]+"_"+stateList[j_col]]
                SOCMat[j_col*s_mult:(j_col+1)*s_mult, i_row*s_mult:(i_row+1)*s_mult] = socMatblock
                SOCMat[i_row*s_mult:(i_row+1)*s_mult, j_col*s_mult:(j_col+1)*s_mult] = -1.*socMatblock
            if (stateList[j_col]+"_"+stateList[i_row] in eomAveTransSOCDict):
                socMatblock = eomAveTransSOCDict[stateList[j_col]+"_"+stateList[i_row]]
                SOCMat[i_row*s_mult:(i_row+1)*s_mult, j_col*s_mult:(j_col+1)*s_mult] = socMatblock
                SOCMat[j_col*s_mult:(j_col+1)*s_mult, i_row*s_mult:(i_row+1)*s_mult] = -1.*socMatblock

# Construct spin-orbit coupling matrix H.
H_soc_au = INVCMTOAU * SOCMat
H = H0Mat + H_soc_au

# Compute spin orbit coupled energies and transformation matrix.
e_vec, c_mat = sorted_diag(H)

# Transform spin and angular momentum matrices.
L_x_orig = np.matrix(angMomMat_x, dtype=np.cdouble)
L_y_orig = np.matrix(angMomMat_y, dtype=np.cdouble)
L_z_orig = np.matrix(angMomMat_z, dtype=np.cdouble)
S_x_orig = np.matrix(spinMat_x, dtype=np.cdouble)
S_y_orig = np.matrix(spinMat_y, dtype=np.cdouble)
S_z_orig = np.matrix(spinMat_z, dtype=np.cdouble)

L_x_trans, L_y_trans, L_z_trans, S_x_trans, S_y_trans, S_z_trans = trans_ls(
    c_mat, L_x_orig, L_y_orig, L_z_orig, S_x_orig, S_y_orig, S_z_orig)

G_x_combined = L_x_trans + S_x_trans
G_y_combined = L_y_trans + S_y_trans
G_z_combined = L_z_trans + S_z_trans

print("calculate g-tensor for target state", targetState)
iTarget = targetStateNo
iup = iTarget*s_mult
idown = (iTarget+1)*s_mult
g_out = compute_state_g_tensor(
    G_x_combined[iup:idown,iup:idown], G_y_combined[iup:idown,iup:idown], 
    G_z_combined[iup:idown,iup:idown])
print(str(g_out))
