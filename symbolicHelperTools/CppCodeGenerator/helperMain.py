import numpy as np
import copy
import sympy as sp
import os

from applyMBBC import *
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

mLabels = [	'm_{000}',
			'm_{100}', 'm_{010}', 'm_{001}',
			'm_{200}', 'm_{020}', 'm_{002}', 'm_{011}', 'm_{101}', 'm_{110}', 
			'm_{111}', 'm_{210}', 'm_{201}', 'm_{120}', 'm_{021}', 'm_{102}', 'm_{012}', 
			'm_{022}', 'm_{202}', 'm_{220}', 'm_{211}', 'm_{121}', 'm_{112}']

meq = 	["rho",
		"rho * ux",
		"rho * uy",
		"rho * uz",
		"(1.f/3.f) * rho + rho * ux * ux",
		"(1.f/3.f) * rho + rho * uy * uy",
		"(1.f/3.f) * rho + rho * uz * uz",
		"rho * uy * uz",
		"rho * ux * uz",
		"rho * ux * uy",
		"0.f",
		"(1.f/3.f) * rho * uy + rho * ux * ux * uy",
		"(1.f/3.f) * rho * uz + rho * ux * ux * uz",
		"(1.f/3.f) * rho * ux + rho * ux * uy * uy",
		"(1.f/3.f) * rho * uz + rho * uy * uy * uz",
		"(1.f/3.f) * rho * ux + rho * ux * uz * uz",
		"(1.f/3.f) * rho * uy + rho * uy * uz * uz",
		"(1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz",
		"(1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz",
		"(1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy",
		"(1.f/3.f) * rho * uy * uz",
		"(1.f/3.f) * rho * ux * uz",
		"(1.f/3.f) * rho * ux * uy"
		]	

meqLin = [
		["1", "0", "0", "0"],
		["0", "1", "0", "0"],
		["0", "0", "1", "0"],
		["0", "0", "0", "1"],
		["1/3", "0", "0", "0"],
		["1/3", "0", "0", "0"],
		["1/3", "0", "0", "0"],
		["0", "0", "0", "0"],
		["0", "0", "0", "0"],
		["0", "0", "0", "0"],
		["0", "0", "0", "0"],
		["0", "0", "1/3", "0"],
		["0", "0", "0", "1/3"],
		["0", "1/3", "0", "0"],
		["0", "0", "0", "1/3"],
		["0", "1/3", "0", "0"],
		["0", "0", "1/3", "0"],
		["1/9", "0", "0", "0"],
		["1/9", "0", "0", "0"],
		["1/9", "0", "0", "0"],
		["0", "0", "0", "0"],
		["0", "0", "0", "0"],
		["0", "0", "0", "0"]
		]

meqOrder = [0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4]

Q = sp.Matrix([
	[1, 0, 0, 0],
	[0, 1, 0, 0],
	[0, 0, 1, 0],
	[0, 0, 0, 1],
	[sp.Rational(1, 3), 0, 0, 0],
	[sp.Rational(1, 3), 0, 0, 0],
	[sp.Rational(1, 3), 0, 0, 0],
	[0, 0, 0, 0],
	[0, 0, 0, 0],
	[0, 0, 0, 0],
	[0, 0, 0, 0],
	[0, 0, sp.Rational(1, 3), 0],
	[0, 0, 0, sp.Rational(1, 3)],
	[0, sp.Rational(1, 3), 0, 0],
	[0, 0, 0, sp.Rational(1, 3)],
	[0, sp.Rational(1, 3), 0, 0],
	[0, 0, sp.Rational(1, 3), 0],
	[sp.Rational(1, 9), 0, 0, 0],
	[sp.Rational(1, 9), 0, 0, 0],
	[sp.Rational(1, 9), 0, 0, 0],
	[0, 0, 0, 0],
	[0, 0, 0, 0],
	[0, 0, 0, 0],
])

Qu = sp.Matrix([
	[0, 0, 0],
	[1, 0, 0],
	[0, 1, 0],
	[0, 0, 1],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, sp.Rational(1, 3), 0],
	[0, 0, sp.Rational(1, 3)],
	[sp.Rational(1, 3), 0, 0],
	[0, 0, sp.Rational(1, 3)],
	[sp.Rational(1, 3), 0, 0],
	[0, sp.Rational(1, 3), 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 0],
])

Qrho = sp.Matrix([
	[1],
	[0],
	[0],
	[0],
	[sp.Rational(1, 3)],
	[sp.Rational(1, 3)],
	[sp.Rational(1, 3)],
	[0],
	[0],
	[0],
	[0],
	[0],
	[0],
	[0],
	[0],
	[0],
	[0],
	[sp.Rational(1, 9)],
	[sp.Rational(1, 9)],
	[sp.Rational(1, 9)],
	[0],
	[0],
	[0],
])

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

	mfk = []
	mfu = []

	def getMoment(mfk, mfu, abc):
		a, b, c = int(abc[0]), int(abc[1]), int(abc[2])
		cResult = np.array(cx, dtype=int)**a * np.asarray(cy, dtype=int)**b * np.asarray(cz, dtype=int)**c
		mfkRow = (cResult[fkIndexes]).tolist()
		mfk.append(mfkRow)
		mfuRow = (cResult[fuIndexes]).tolist()
		mfu.append(mfuRow)

	for abc in abcs:
		getMoment(mfk, mfu, abc)
		
	mfk = sp.Matrix(mfk)
	mfu = sp.Matrix(mfu)

	### First part: MBBC
	rowsToKeep = [0]
	i = 1
	while len(rowsToKeep) < len(fu):
		if i >= len(mfu):
			print("getChosenMoments Error: Ran out of independent groups")
			return None
		previousRank = mfu.row(rowsToKeep).rank()
		attemptedRank = mfu.row(rowsToKeep + [i]).rank()
		if attemptedRank > previousRank:
			rowsToKeep.append(i)
		i += 1
	identityMatrix = sp.eye(23)
	S =  identityMatrix.row(rowsToKeep)
	
	### Second part: restoring rho, ux, uy, uz
	nullspace = mfu.T.nullspace()
	def nonzero_order_sum(v):
		key = 0
		for i, val in enumerate(v):
			if val != 0:
				key += meqOrder[i]
		return key
	sorted_nullspace = sorted(nullspace, key=nonzero_order_sum)
	C = sp.Matrix.hstack(*sorted_nullspace)
	rowsToKeep = [0]
	i = 1
	while len(rowsToKeep) < 4:
		if i >= len(C.T * Q):
			print("nullspace Error: Ran out of independent columns")
			return None
		
		previousRank = (C.T * Q).row(rowsToKeep).rank()
		attemptedRank =  (C.T * Q).row(rowsToKeep + [i]).rank()
		if attemptedRank > previousRank:
			rowsToKeep.append(i)
		i += 1
	identityMatrix = sp.eye(len((C.T * Q).tolist()))
	Sq =  identityMatrix.row(rowsToKeep)
	
	return fk, fu, mfk, mfu, S, Sq, C

os.makedirs("results", exist_ok=True)

allNormals = [	[1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1],
				[0,1,1], [1,0,1], [1,1,0], [0,-1,-1], [-1,0,-1], [-1,-1,0],
				[0,1,-1], [1,0,-1], [1,-1,0], [0,-1,1], [-1,0,1], [-1,1,0],
				[1,1,1], [-1,-1,-1], 
				[-1,1,1], [1,-1,1], [1,1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1] ]

"""
#### applyMBBC
allLines = []
for i, normal in enumerate(allNormals):
	fk, fu, mLabels, mfk, mfu, chosenMoments, K, U, UInv, meqDict, meqLinDict = getBCdata(normal)
	allLines += applyMBBC(i, normal, fk, fu, mfk, mfu, chosenMoments, K, U, UInv, meqDict)

with open("results/applyVelocityInlet.hpp", "w") as file:
	file.write("\n".join(allLines))
"""

#### Latex Tables
normal = [1, -1, 1]
fk, fu, mfk, mfu, S, Sq, C = getBCdata(normal)
latexCode = getLatex(f, normal, fk, fu, mfk, mfu, S, Sq, C, Q)

