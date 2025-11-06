#########################################################
# THIS PROGRAM GENERATES A "JSON" FILE                  #
# (WITH SOCS, ENERGIES, SPIN AND ORBITAL ANGULAR        #
# MOMENTUMS OF ALL THE STATES) FROM THE EOM Q-CHEM      #
# OUTPUT.                                               #
#########################################################

__author__ = 'Sven Kähler'

import sys, re, string
import numpy as np
import json
from decimal import Decimal

# flag for extra print
print_intermediates = False

# open EOM output file_ms_notnull
lowestEtargetState = False
if len(sys.argv) == 3:
    targetState = sys.argv[2]
else:
    lowestEtargetState = True
unsortStateList = []
stateList = []

with open(sys.argv[1], 'r') as f:
    contents = f.readlines()

# regular expression defitions
eomRe = re.compile(r'\s+Solving for EOM(?:IP|EA)-CCSD(?:\/MP2)? (\w[\d\'\w]*) transitions\.')
eomStateRe = re.compile(r'\s+EOM(?:IP|EA) transition')
symLabelRe = re.compile(r'(\d+/[A-B][\d\'\w]*)')
eomStateEnergyRe = re.compile(r'\s+Total energy = (-?\d+\.\d*) a\.u\.')
realValueRe = re.compile(r'(-?\d+\.\d*)')
statePropRe = re.compile(r'\s+Excited state properties for  EOM(?:IP|EA)-(?:CCSD|MP2) transition (\d+/\w[\d\'\w]*)')
stateAngMomRe = re.compile(r'\s+Angular momentum \(a\.u\.\) against gauge origin\:')
stateAngMomVecRe = re.compile(r'\s+\(X (-?\d+\.\d+)i, Y (-?\d+\.\d+)i, Z (-?\d+\.\d+)i\)')
totalAngMomRe = re.compile(r'\s+<S\^2> = (\d+\.\d+)')
transPropRe = re.compile(r'Start computing the transition properties')
stateARe = re.compile(r' State A: eom(?:ip|ea)_(?:ccsd|mp2)/(?:a|b): (\d+/\w[\d\'\w]*)')
stateBRe = re.compile(r' State B: eom(?:ip|ea)_(?:ccsd|mp2)/(?:a|b): (\d+/\w[\d\'\w]*)')
transAngMomRe = re.compile(r' Transition angular momentum against gauge origin \(a\.u\.\)\:')
transAngMomVecABRe = re.compile(r'\s*A->B:  \(X (-?\d\.\d+)i, Y (-?\d\.\d+)i, Z (-?\d\.\d+)i\)')
transAngMomVecBARe = re.compile(r'\s*B->A:  \(X (-?\d\.\d+)i, Y (-?\d\.\d+)i, Z (-?\d\.\d+)i\)')
socTransitionABRe = re.compile(r'   A(Ket)->B(Bra) transition SO matrices')
socTransitionBARe = re.compile(r'   B(Ket)->A(Bra) transition SO matrices')
aveTransSOCRe = re.compile(r'\s*Arithmetically averaged transition SO matrices')
mfSOCRe = re.compile(r'\s*Mean-field SO \(cm-1\)')
oneelSOCRe = re.compile(r'\s*One-electron SO \(cm-1\)')
matElemRe = re.compile(r'\s*Actual matrix elements:')
indexKetRowRe = re.compile(r'(?:\s*\|Sz=(-?\d\.\d+)>)+')
socMatElemRe = re.compile(r'\s+<Sz=(-?\d\.\d+)\|(?:\(-?\d+\.\d+,-?\d+\.\d+\))+')
complexPairRe = re.compile(r'\((-?\d+\.\d+),(-?\d+\.\d+)\)')
indexBraRe = re.compile(r'\s+<Sz=(-?\d\.\d+)|')
transPropDivRe = re.compile(r'\s+---+')

# initialize lists for captured properties
eomStates = []
eomStateLabels = []
eomStateEnergies = []
eomStateEnergiesDict = {}
eomStateTotalAngMomDict = {}
eomStateAngMomVecs = []
eomStateAngMomVecsDict = {} 
eomTransAngMomVecDict = {}
eomAveTransSOCDict = {}
nEOMStates = 0
foundEOM = False
foundState = False
foundEOMprop = False
foundStateEOMProp = False
foundStateAngMom = False
foundEOMtransProp = False
foundALabel = False
foundBLabel = False
foundStateA = False
foundStateB = False
foundTransAngMom = False
foundTransSOCAB = False
foundTransSOCBA = False
foundAveTransSOC = False
foundMatElem = False
foundMFSOC = False
found1elSOC = False

def printFoundStatePropMessage():
    print("found eom prop")
    print("eomStateLabels and eomStateEnergies should be complete")
    print("eomStateLabels:")
    print(str(eomStateLabels))
    print("eomStateEnergies:")
    print(str(eomStateEnergies))
    print("eomStateEnergiesDict:")
    print(str(eomStateEnergiesDict))

def printFoundTransPropMessage():
    print("found trans prop")
    print("stateAngMomVecsDict should be complete")
    print("eomStateAngMomVecsDict:")
    print(str(eomStateAngMomVecsDict))

# iteratre through file_ms_notnull search for state properties and start of transtion properties section
for line in contents:
    if eomRe.match(line):
        foundEOM = True
        foundState = False
        foundStateEOMProp = False
    # when EOM section of output found look excited state and symmetry label
    if foundEOM:
        if eomStateRe.match(line):
            match = eomStateRe.match(line)
            symLabel = re.findall(symLabelRe, line)[0]
            eomStateLabels.append(re.findall(symLabelRe, line))
            nEOMStates = nEOMStates + 1
            foundState = True
    # when excited state information found, search for excitation energy
    if foundState:
        if eomStateEnergyRe.match(line):
            energy_list = [float(energy) for energy in re.findall(realValueRe, line)]
            eomStateEnergies.append(energy_list)
            eomStateEnergiesDict[symLabel] = energy_list
    if statePropRe.match(line):
        symLabel = re.findall(symLabelRe, line)[0]
        if(not foundStateEOMProp and print_intermediates):
            printFoundStatePropMessage()
        foundStateEOMProp = True
    if (foundStateEOMProp and stateAngMomRe.match(line)):
        foundStateAngMom = True
    if (foundStateAngMom and stateAngMomVecRe.match(line)):
        foundStateAngMom = False
        eomStateAngMomVecs.append(re.findall(realValueRe, line))
        eomStateAngMomVecsDict[symLabel] = re.findall(realValueRe, line)
    # when EOM section found (not in SCF or MP2/CCSD) look for total spin angular momenta
    if (foundStateEOMProp and totalAngMomRe.match(line)):
        eomStateTotalAngMomDict[symLabel] = float(re.findall(realValueRe, line)[0])
    if (transPropRe.match(line)):
        if (not foundEOMtransProp and print_intermediates):
            printFoundTransPropMessage()
        foundEOMtransProp = True
        break;


foundState = False
for line in contents:
    if stateARe.match(line):
        symLabelA = re.findall(symLabelRe, line)[0]
        foundALabel = True
        foundStateA = True
    if stateBRe.match(line):
        symLabelB = re.findall(symLabelRe, line)[0]
        foundBLabel = True
        foundStateB = True
    if (foundStateA and foundStateB):
        if transAngMomRe.match(line):
            foundTransAngMom = True
        if (foundTransAngMom and transAngMomVecABRe.match(line)):
            transAngMomVecAB = [Decimal(i) for i in re.findall(realValueRe, line)]
        if (foundTransAngMom and transAngMomVecBARe.match(line)):
            transAngMomVecBA = [Decimal(i) for i in re.findall(realValueRe, line)]
        if (oneelSOCRe.match(line)):
            found1elSOC = True
            foundMFSOC = False
        if (mfSOCRe.match(line)):
            foundMFSOC = True
            found1elSOC = False
        if (aveTransSOCRe.match(line)):
            foundAveTransSOC = True;
        if (foundAveTransSOC and matElemRe.match(line)):
            foundMatElem = True;
        if (foundMatElem and indexKetRowRe.match(line)):
            ketIndices = re.findall(realValueRe, line)
            if found1elSOC:
                oneelSOCElems = []
                oneelSOCColIndices = []
            if foundMFSOC:
                mfSOCElems = []
                mfSOCColIndices = []
        if (foundMatElem and found1elSOC and socMatElemRe.match(line)):
            match = indexBraRe.match(line)
            oneelSOCColIndices.extend(re.findall(realValueRe, match.group()))
            oneelSOCElems.append(re.findall(complexPairRe, line))
        if (foundMatElem and foundMFSOC and socMatElemRe.match(line)):
            match = indexBraRe.match(line)
            mfSOCColIndices.extend(re.findall(realValueRe, match.group()))
            mfSOCElems.append([complex(float(p[0]),float(p[1])) for p in re.findall(complexPairRe, line)])
        # once end of transition properties section reached, process parsed properties:
        # average AB and BA angular momentum transition vector
        # save both, angular momemta and spin-orbit couplings to dictionaries
        if (transPropDivRe.match(line)):
            foundALabel = False
            foundBLabel = False
            foundStateA = False
            foundStateB = False
            foundMatElem = False
            foundAveTransSOC = False
            foundTransAngMom = False
            aveTransAngMomVec = np.array(
                [complex(0.,(ab-ba)/Decimal(2)) for ab,ba in zip(transAngMomVecAB,transAngMomVecBA)])
            eomTransAngMomVecDict[symLabelA + "_" + symLabelB] = aveTransAngMomVec
            eomAveTransSOCDict[symLabelA + "_" + symLabelB] = np.matrix(mfSOCElems)

def matrix_to_string_list(matrix):
    return [[str(elem) for elem in line] for line in matrix.tolist()]
    
# prepare data for json export - convert complex numbers to strings and matrices and vectors to lists
eomTransAngMomListDict = {x : [str(elem) for elem in list(y)] for x, y in eomTransAngMomVecDict.items()}
eomAveTransSOCListDict = {x : matrix_to_string_list(y) for x, y in eomAveTransSOCDict.items()}

output_dict = {
    "stateEnergiesDict" : eomStateEnergiesDict,
    "stateTotalAngMomDict" : eomStateTotalAngMomDict,
    "transAngMomListDict" : eomTransAngMomListDict,
    "aveTransSOCListDict" : eomAveTransSOCListDict
}

outfile_name = sys.argv[1] + ".json"
with open(outfile_name, 'w') as f:
    json.dump(output_dict, f, 
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)

