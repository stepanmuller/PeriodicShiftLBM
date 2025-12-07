import numpy as np
import copy
import sympy as sp
import os

from readKnownDistributionsCpp import *
from periodicBCCpp import *
from runMBBCCpp import *
from uInletConsistencyCpp import *
from pOutletConsistencyCpp import *
from getLatex import *

f    =  [ f"f_{{{i}}}" for i in range(27) ]
cx   =  [ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 ]
cy   =  [ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 ]
cz   =  [ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 ]
w    =  [ "8/27", 
		  "2/27", "2/27", "2/27 ", "2/27", "2/27", "2/27", 
		  "1/54", "1/54", "1/54", "1/54", "1/54", "1/54", "1/54", "1/54", "1/54", "1/54", "1/54", "1/54", 
		  "1/216", "1/216", "1/216", "1/216", "1/216", "1/216", "1/216", "1/216" ]
abcs =	[ "000", 
		  "100", "010", "001", 
		  "200", "020", "002", "011", "101", "110", 
		  "111", "210", "201", "120", "021", "102", "012", 
		  "022", "202", "220", "211", "121", "112" ]

def getBCdata(normal):
	fk = []
	fu = []

	cxyz = np.vstack((np.array(cx), np.array(cy), np.array(cz))).T
	identifier = np.min(cxyz * normal, axis=1)
	fkIndexes = np.where(identifier >= 0)[0]
	fuIndexes = np.where(identifier < 0)[0]
	for i in fkIndexes:
		fk.append(f[i])
	for i in fuIndexes:
		fu.append(f[i])

	mLabels = []
	mfkRows = []
	mfuRows = []

	def getMoment(mLabels, mfkRows, mfuRows, abc):
		label = "m_{" + abc + "}"
		mLabels.append(label)
		a, b, c = int(abc[0]), int(abc[1]), int(abc[2])
		cResult = np.array(cx, dtype=int)**a * np.asarray(cy, dtype=int)**b * np.asarray(cz, dtype=int)**c
		mfkRow = (cResult[fkIndexes]).tolist()
		mfkRows.append(mfkRow)
		mfuRow = (cResult[fuIndexes]).tolist()
		mfuRows.append(mfuRow)

	for abc in abcs:
		getMoment(mLabels, mfkRows, mfuRows, abc)

	uniqueMfuRows = []
	mGroups = []

	for mfuRow, mLabel in zip(mfuRows, mLabels):
		if mfuRow not in uniqueMfuRows:
			uniqueMfuRows.append(mfuRow)
			mGroups.append([mLabel])
		else:
			index = uniqueMfuRows.index(mfuRow)
			mGroups[index].append(mLabel)

	def getChosenMoments():
		chosenMoments = []
		U = []
		K = []
		numberOfEquations = len(fu)
		i = 0
		while len(chosenMoments) < numberOfEquations:
			if i >= len(mfuRows):
				print("getChosenMoments Error: Ran out of independent groups")
				return None
			MfuRow = uniqueMfuRows[i]
			previousRank = np.linalg.matrix_rank(np.array(U))
			attemptedU = copy.deepcopy(U)
			attemptedU.append(MfuRow)
			attemptedRank = np.linalg.matrix_rank(np.array(attemptedU))
			if attemptedRank > previousRank:
				U.append(MfuRow)
				chosenMoment = mGroups[i][0]
				chosenMoments.append(chosenMoment)
				globalIndex = mLabels.index(chosenMoment)
				K.append(mfkRows[globalIndex])
			i += 1
		return chosenMoments, U, K

	chosenMoments, U, K = getChosenMoments()

	USympy = sp.Matrix(U)
	UInvSympy = USympy.inv()
	UInv = UInvSympy.tolist()

	def getMeqDict():
		meqDict = {}
		meqDict["m_{000}"] = "rho"
		meqDict["m_{100}"] = "rho * ux"
		meqDict["m_{010}"] = "rho * uy"
		meqDict["m_{001}"] = "rho * uz"
		meqDict["m_{200}"] = "(1.f/3.f) * rho + rho * ux * ux"
		meqDict["m_{020}"] = "(1.f/3.f) * rho + rho * uy * uy"
		meqDict["m_{002}"] = "(1.f/3.f) * rho + rho * uz * uz"
		meqDict["m_{011}"] = "rho * uy * uz"
		meqDict["m_{101}"] = "rho * ux * uz"
		meqDict["m_{110}"] = "rho * ux * uy"
		meqDict["m_{111}"] = "0.f"
		meqDict["m_{210}"] = "(1.f/3.f) * rho * uy + rho * ux * ux * uy"
		meqDict["m_{201}"] = "(1.f/3.f) * rho * uz + rho * ux * ux * uz"
		meqDict["m_{120}"] = "(1.f/3.f) * rho * ux + rho * ux * uy * uy"
		meqDict["m_{021}"] = "(1.f/3.f) * rho * uz + rho * uy * uy * uz"
		meqDict["m_{102}"] = "(1.f/3.f) * rho * ux + rho * ux * uz * uz"
		meqDict["m_{012}"] = "(1.f/3.f) * rho * uy + rho * uy * uz * uz"
		meqDict["m_{022}"] = "(1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz"
		meqDict["m_{202}"] = "(1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz"
		meqDict["m_{220}"] = "(1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy"
		meqDict["m_{211}"] = "(1.f/3.f) * rho * uy * uz"
		meqDict["m_{121}"] = "(1.f/3.f) * rho * ux * uz"
		meqDict["m_{112}"] = "(1.f/3.f) * rho * ux * uy"
		return meqDict

	meqDict = getMeqDict()
	
	return fk, fu, mLabels,	mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict

os.makedirs("results", exist_ok=True)

allNormals = [	[1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1],
				[0,1,1], [1,0,1], [1,1,0], [0,-1,-1], [-1,0,-1], [-1,-1,0],
				[0,1,-1], [1,0,-1], [1,-1,0], [0,-1,1], [-1,0,1], [-1,1,0],
				[1,1,1], [-1,-1,-1], 
				[-1,1,1], [1,-1,1], [1,1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1] ]

#### Read known distributions
allLines = []
for i, normal in enumerate(allNormals):
	fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict = getBCdata(normal)
	allLines += readKnownDistributionsCpp(i, normal, fk, fu, chosenMoments, K, UInv, meqDict)

with open("results/readFKnown.hpp", "w") as file:
    file.write("\n".join(allLines))

#### MBBC
allLines = []
for i, normal in enumerate(allNormals):
	fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict = getBCdata(normal)
	allLines += runMBBCCpp(i, normal, fk, fu, chosenMoments, K, UInv, meqDict)

with open("results/applyMBBC.hpp", "w") as file:
    file.write("\n".join(allLines))

allNormals = [	[1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1] ]

#### uInletConsistency
allLines = []
for i, normal in enumerate(allNormals):
	fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict = getBCdata(normal)
	allLines += uInletConsistencyCpp(i, normal, fk, mfkRows)

with open("results/getInletConsistency.hpp", "w") as file:
    file.write("\n".join(allLines))
    
#### pOutletConsistency
allLines = []
for i, normal in enumerate(allNormals):
	fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict = getBCdata(normal)
	allLines += pOutletConsistencyCpp(i, normal, fk, mLabels, mfkRows)

with open("results/getOutletConsistency.hpp", "w") as file:
    file.write("\n".join(allLines))
    
#### Periodic BC
allLines = []
for i, normal in enumerate(allNormals):
	fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict = getBCdata(normal)
	allLines += periodicBCCpp(i, normal, fu)

with open("results/applyPeriodicBC.hpp", "w") as file:
    file.write("\n".join(allLines))



#### Latex Tables
normal = [0, 0, 1]
fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict = getBCdata(normal)
latexCode = getLatex(f, cx, cy, cz, w, normal, fk, fu, mLabels,	mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict)

