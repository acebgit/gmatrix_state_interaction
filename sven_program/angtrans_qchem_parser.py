import sys, re, string
import numpy as np
import json
from dicttoxml import dicttoxml
from decimal import Decimal

lowestEtargetState = False
print("len(sys.argv)= " + str(len(sys.argv)))
if len(sys.argv) == 3:
    targetState = sys.argv[2]
else:
    lowestEtargetState = True
unsortStateList = []
stateList = []

with open(sys.argv[1], 'r') as f:
    contents = f.readlines()

# find start of excited state results in file
#eomipRe = re.compile(r'\s+Solving for EOMIP-CCSD (\w[\d\'\w]*) transitions\.')
#eomRe = re.compile(r'\s+Solving for EOMIP-CCSD (\w[\d\'\w]*) transitions\.')
#eomRe = re.compile(r'\s+Solving for EOM(?:IP|EA)-(?:CCSD|MP2) (\w[\d\'\w]*) transitions\.')
eomRe = re.compile(r'\s+Solving for EOM(?:IP|EA)-CCSD(?:\/MP2)? (\w[\d\'\w]*) transitions\.')
#eomipStateRe = re.compile(r'\s+EOMIP transition (\d+/\w[\d\'\w]*)')
#eomipStateRe = re.compile(r'\s+EOMIP transition')
#eomStateRe = re.compile(r'\s+EOMIP transition')
eomStateRe = re.compile(r'\s+EOM(?:IP|EA) transition')
#symLabelRe = re.compile(r'(\d+/\w[\d\'\w]*)')
symLabelRe = re.compile(r'(\d+/[A-B][\d\'\w]*)')
#eomipStateEnergyRe = re.compile(r'\s+Total energy = (-?\d+\.\d*) a\.u\.')
eomStateEnergyRe = re.compile(r'\s+Total energy = (-?\d+\.\d*) a\.u\.')
realValueRe = re.compile(r'(-?\d+\.\d*)')
#statePropRe = re.compile(r'\s+Excited state properties for  EOMIP-CCSD transition (\d+/\w[\d\'\w]*)')
#statePropRe = re.compile(r'\s+Excited state properties for  EOMIP-CCSD transition (\d+/\w[\d\'\w]*)')
statePropRe = re.compile(r'\s+Excited state properties for  EOM(?:IP|EA)-(?:CCSD|MP2) transition (\d+/\w[\d\'\w]*)')
stateAngMomRe = re.compile(r'\s+Angular momentum \(a\.u\.\) against gauge origin\:')
stateAngMomVecRe = re.compile(r'\s+\(X (-?\d+\.\d+)i, Y (-?\d+\.\d+)i, Z (-?\d+\.\d+)i\)')
totalAngMomRe = re.compile(r'\s+<S\^2> = (\d+\.\d+)')
transPropRe = re.compile(r'Start computing the transition properties')
#stateARe = re.compile(r' State A: eomip_ccsd/a: (\d+/\w[\d\'\w]*)')
#stateARe = re.compile(r' State A: eomip_mp2/a: (\d+/\w[\d\'\w]*)')
#stateARe = re.compile(r' State A: eomip_(?:ccsd|mp2)/a: (\d+/\w[\d\'\w]*)')
#stateARe = re.compile(r' State A: eom(?:ip|ea)_(?:ccsd|mp2)/a: (\d+/\w[\d\'\w]*)')
stateARe = re.compile(r' State A: eom(?:ip|ea)_(?:ccsd|mp2)/(?:a|b): (\d+/\w[\d\'\w]*)')
#stateBRe = re.compile(r' State B: eomip_ccsd/a: (\d+/\w[\d\'\w]*)')
#stateBRe = re.compile(r' State B: eomip_mp2/a: (\d+/\w[\d\'\w]*)')
#stateBRe = re.compile(r' State B: eomip_(?:ccsd|mp2)/a: (\d+/\w[\d\'\w]*)')
#stateBRe = re.compile(r' State B: eom_(?:ccsd|mp2)/a: (\d+/\w[\d\'\w]*)')
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
#realPartComplexValueRe = re.compile(r'\((-?\d+\.\d+),-?\d+\.\d+\)')
#imaginaryPartComplexValueRe = re.compile(r'\(-?\d+\.\d+,(-?\d+\.\d+)\)')
indexBraRe = re.compile(r'\s+<Sz=(-?\d\.\d+)|')
transPropDivRe = re.compile(r'\s+---+')

#eomipStates = []
eomStates = []
#eomipStateLabels = []
eomStateLabels = []
#eomipStateEnergies = []
eomStateEnergies = []
#eomipStateEnergiesDict = {}
eomStateEnergiesDict = {}
#eomipStateTotalAngMomDict = {}
eomStateTotalAngMomDict = {}
#eomipStateAngMomVecs = []
eomStateAngMomVecs = []
#eomipStateAngMomVecsDict = {} 
eomStateAngMomVecsDict = {} 
#eomipTransAngMomVecDict = {}
eomTransAngMomVecDict = {}
#eomipAveTransSOCDict = {}
eomAveTransSOCDict = {}
#nEOMIPStates = 0
nEOMStates = 0
#foundEOMIP = False
foundEOM = False
foundState = False
#foundEOMIPprop = False
foundEOMprop = False
#foundStateEOMIPProp = False
foundStateEOMProp = False
foundStateAngMom = False
#foundEOMIPtransProp = False
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
    #print("stateList (= list of states with same S^2 as target) should be complete")
    #print("unsortStateList:")
    #print(unsortStateList)

def printFoundTransPropMessage():
    print("found trans prop")
    print("stateAngMomVecsDict should be complete")
    print("eomStateAngMomVecsDict:")
    print(str(eomStateAngMomVecsDict))

for line in contents:
    #if eomipRe.match(line):
    if eomRe.match(line):
        #foundEOMIP = True
        foundEOM = True
        foundState = False
        #foundStateEOMIPProp = False
        foundStateEOMProp = False
        #log.debug("foundEOMIP: ")
        #log.debug("foundEOM: ")
        #log.debug(line)
    #if foundEOMIP:
    if foundEOM:
        #print("looking for eomipStateRe")
        #if eomipStateRe.match(line):
        if eomStateRe.match(line):
            #match = eomipStateRe.match(line)
            match = eomStateRe.match(line)
            symLabel = re.findall(symLabelRe, line)[0]
            #log.debug("symLabel: " + symLabel + " in stateLabels")
            #eomipStateLabels.append(re.findall(symLabelRe, line))
            eomStateLabels.append(re.findall(symLabelRe, line))
            #nEOMIPStates = nEOMIPStates + 1
            nEOMStates = nEOMStates + 1
            #log.debug("found eom state: ")
            #log.debug(line)
            #log.debug("number of found eom states: " + str(nEOMIPStates))
            #log.debug("number of found eom states: " + str(nEOMStates))
            #log.debug("eom state labels: " + str(eomipStateLabels))
            #log.debug("eom state labels: " + str(eomStateLabels))
            foundState = True
    if foundState:
        #if eomipStateEnergyRe.match(line):
        if eomStateEnergyRe.match(line):
            #match = eomStateEnergyRe.match(line)
            energy_list = [float(energy) for energy in re.findall(realValueRe, line)]
            eomStateEnergies.append(energy_list)
            eomStateEnergiesDict[symLabel] = energy_list
            print("found eom state energy")
            print(line)
            #print("eom state energies: " + str(eomStateEnergies))
            print("eom state energies dict:", eomStateEnergiesDict[symLabel])
    if statePropRe.match(line):
        symLabel = re.findall(symLabelRe, line)[0]
        #if(not foundStateEOMIPProp):
        if(not foundStateEOMProp):
            printFoundStatePropMessage()
        #foundStateEOMIPProp = True
        foundStateEOMProp = True
        #log.debug("found eom prop")
        print("found eom prop")
        #log.debug(line)
    #if (foundStateEOMIPProp and stateAngMomRe.match(line)):
    if (foundStateEOMProp and stateAngMomRe.match(line)):
        foundStateAngMom = True
        #log.debug("found state ang mom")
        #log.debug(line)
    if (foundStateAngMom and stateAngMomVecRe.match(line)):
        foundStateAngMom = False
        #log.debug("found state ang mom vec")
        #log.debug(line)
        #eomipStateAngMomVecs.append(re.findall(realValueRe, line))
        eomStateAngMomVecs.append(re.findall(realValueRe, line))
        #eomipStateAngMomVecsDict[symLabel] = re.findall(realValueRe, line)
        eomStateAngMomVecsDict[symLabel] = re.findall(realValueRe, line)
        #log.debug("printing symLabel and AngMomVec: " + symLabel + ", " + str(re.findall(realValueRe, line)))
    #if (foundStateEOMIPProp and totalAngMomRe.match(line)):
    if (foundStateEOMProp and totalAngMomRe.match(line)):
        print("found total ang mom")
        print(line)
        #eomipStateTotalAngMomDict[symLabel] = re.findall(realValueRe, line)
        eomStateTotalAngMomDict[symLabel] = float(re.findall(realValueRe, line)[0])
        print("symLabel {} and totalAngMom {}".format(symLabel, eomStateTotalAngMomDict[symLabel]))
    if (transPropRe.match(line)):
        #if (not foundEOMIPtransProp):
        if (not foundEOMtransProp):
            printFoundTransPropMessage()
        #foundEOMIPtransProp = True
        foundEOMtransProp = True
        #log.debug("found trans prop")
        #log.debug(line)
        break;

foundState = False
for line in contents:
    if stateARe.match(line):
        #log.debug("found state A")
        #log.debug(line)
        symLabelA = re.findall(symLabelRe, line)[0]
        foundALabel = True
        #log.debug("symLabelA: " + symLabelA + " in stateLabels")
        foundStateA = True
        #log.debug("symLabelA: " + symLabelA)
    if stateBRe.match(line):
        #log.debug("found state B")
        #log.debug(line)
        symLabelB = re.findall(symLabelRe, line)[0]
        foundBLabel = True
        #log.debug("symLabelB: " + symLabelB + " in stateLabels")
        foundStateB = True
        #log.debug("symLabelB: " + symLabelB)
    if (foundStateA and foundStateB):
        if transAngMomRe.match(line):
            foundTransAngMom = True
            #log.debug("found trans ang mom")
        if (foundTransAngMom and transAngMomVecABRe.match(line)):
            transAngMomVecAB = [Decimal(i) for i in re.findall(realValueRe, line)]
            #transAngMomVecAB = [complex(0.,float(i)) for i in re.findall(realValueRe, line)]
            #log.debug("transAngMomVecAB: " + str(transAngMomVecAB))
            #log.debug(line)
        if (foundTransAngMom and transAngMomVecBARe.match(line)):
            transAngMomVecBA = [Decimal(i) for i in re.findall(realValueRe, line)]
            #transAngMomVecBA = [complex(0.,float(i)) for i in re.findall(realValueRe, line)]
            #log.debug("transAngMomVecBA: " + str(transAngMomVecBA))
            #log.debug(line)
        #if (socTransitionABRe.match(line)):
        #    foundTransSOCAB = True
        #    foundTransSOCBA = False
        #    print("found trans SOC AB")
        #if (socTransitionBARe.match(line)):
        #    foundTransSOCBA = True
        #    foundTransSOCAB = False
        #    print("found trans SOC BA")
        if (oneelSOCRe.match(line)):
            found1elSOC = True
            foundMFSOC = False
            #log.debug("found 1el SOC")
        if (mfSOCRe.match(line)):
            foundMFSOC = True
            found1elSOC = False
            #log.debug("found mf SOC")
        if (aveTransSOCRe.match(line)):
            foundAveTransSOC = True;
            #log.debug("found averaged trans SOC")
        if (foundAveTransSOC and matElemRe.match(line)):
            foundMatElem = True;
            #log.debug("found mat elem")
        if (foundMatElem and indexKetRowRe.match(line)):
            #log.debug("found index ket row")
            ketIndices = re.findall(realValueRe, line)
            #log.debug(line)
            #log.debug("ket indices: " + str(ketIndices))
            if found1elSOC:
                oneelSOCElems = []
                oneelSOCColIndices = []
            if foundMFSOC:
                mfSOCElems = []
                mfSOCColIndices = []
        if (foundMatElem and found1elSOC and socMatElemRe.match(line)):
            #log.debug("found 1el SOC mat elems")
            #log.debug(line)
            match = indexBraRe.match(line)
            oneelSOCColIndices.extend(re.findall(realValueRe, match.group()))
            oneelSOCElems.append(re.findall(complexPairRe, line))
            #log.debug("1el SOC col indices: " + str(oneelSOCColIndices))
            #log.debug("1el SOC mat elems: " + str(oneelSOCElems))
        if (foundMatElem and foundMFSOC and socMatElemRe.match(line)):
            #log.debug("found mf SOC mat elems")
            #log.debug(line)
            match = indexBraRe.match(line)
            mfSOCColIndices.extend(re.findall(realValueRe, match.group()))
            mfSOCElems.append([complex(float(p[0]),float(p[1])) for p in re.findall(complexPairRe, line)])
            #log.debug("mf SOC col indices: " + str(mfSOCColIndices))
            #log.debug("mf SOC mat elems: " + str(mfSOCElems))
        if (transPropDivRe.match(line)):
            #log.debug("found trans prop div line")
            foundALabel = False
            foundBLabel = False
            foundStateA = False
            foundStateB = False
            foundMatElem = False
            foundAveTransSOC = False
            foundTransAngMom = False
            # Taking stock (and printing) what has been captured
            #log.debug("Captured data for this state pair")
            #log.debug("State label A: " + str(symLabelA))
            #log.debug("State label B: " + str(symLabelB))
            #log.debug("Trans ang mom AB: " + str(transAngMomVecAB))
            #log.debug("Trans ang mom BA: " + str(transAngMomVecBA))
            #log.debug("ket spin indices: " + str(ketIndices))
            #log.debug("bra spin indices (mf): " + str(mfSOCColIndices))
            #log.debug("mf SOC mat elems: " + str(mfSOCElems))
            #log.debug("state pair to be included: " + symLabelA + " " + symLabelB)
            aveTransAngMomVec = np.array(
                [complex(0.,(ab-ba)/Decimal(2)) for ab,ba in zip(transAngMomVecAB,transAngMomVecBA)])
            print("aveTrangAngMomVec", aveTransAngMomVec)
            #aveTransAngMomVec = (np.array(transAngMomVecAB) - np.array(transAngMomVecBA))/2.0
            #eomipTransAngMomVecDict[symLabelA + "_" + symLabelB] = aveTransAngMomVec
            eomTransAngMomVecDict[symLabelA + "_" + symLabelB] = aveTransAngMomVec
            #eomipAveTransSOCDict[symLabelA + "_" + symLabelB] = np.matrix(mfSOCElems)
            eomAveTransSOCDict[symLabelA + "_" + symLabelB] = np.matrix(mfSOCElems)

def matrix_to_string_list(matrix):
    return [[str(elem) for elem in line] for line in matrix.tolist()]
    

print("finished parsing")
#print("eomStates:\n" + str(eomStates))
#print("eomStateLabels:\n" + str(eomStateLabels))
#print("eomStateEnergies:\n" + str(eomStateEnergies))
#print("eomStateEnergiesDict:\n" + str(eomStateEnergiesDict))
print("eomStateEnergiesDict:\n" + str(eomStateEnergiesDict))
print("eomStateTotalAngMomDict:\n" + str(eomStateTotalAngMomDict))
#print("eomStateAngMomVecsDict:\n" + str(eomStateAngMomVecsDict))
eomTransAngMomListDict = {x : [str(elem) for elem in list(y)] for x, y in eomTransAngMomVecDict.items()}
print("eomTransAngMomVecDict:\n" + str(eomTransAngMomListDict))
eomAveTransSOCListDict = {x : matrix_to_string_list(y) for x, y in eomAveTransSOCDict.items()}
print("eomAveTransSOCDict:\n" + str(eomAveTransSOCListDict))

output_dict = {
    "eomStateEnergiesDict" : eomStateEnergiesDict,
    "eomStateTotalAngMomDict" : eomStateTotalAngMomDict,
    "eomTransAngMomListDict" : eomTransAngMomListDict,
    "eomAveTransSOCListDict" : eomAveTransSOCListDict
}

# If python version >= 3.9 sys.argv[1].removesuffix(".out") would be nice.

outfile_name = sys.argv[1] + ".json"
with open(outfile_name, 'w') as f:
    json.dump(output_dict, f, 
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)

xmloutfile_name = sys.argv[1] + ".xml"
xml = dicttoxml(output_dict)
with open(xmloutfile_name, 'w') as f:
    f.write(str(xml))