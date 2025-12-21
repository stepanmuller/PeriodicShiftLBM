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
		
	mfk = np.array(mfk)
	mfu = np.array(mfu)

	def getSelector():
		selector = [0]
		numberOfEquations = len(fu)
		i = 1
		while len(selector) < numberOfEquations:
			if i >= len(mfu):
				print("getChosenMoments Error: Ran out of independent groups")
				return None
			previousRank = np.linalg.matrix_rank(mfu[selector])
			attemptedRank = np.linalg.matrix_rank(mfu[(selector + [i])])
			if attemptedRank > previousRank:
				selector.append(i)
			i += 1
		return selector

	selector = getSelector()
	
	K = mfk[selector]
	U = mfu[selector]
	
	USympy = sp.Matrix(U)
	UInvSympy = USympy.inv()
	UInv = UInvSympy.tolist()
	
	### Now second part - restoring macroscopic quantities
	
	mfuSympy = sp.Matrix(mfu)
	nullspace = mfuSympy.T.nullspace()
	def nonzero_index_sum(v):
		return sum(i for i, val in enumerate(v) if val != 0)
	sorted_nullspace = sorted(nullspace, key=nonzero_index_sum)
	mct = sp.Matrix.hstack(*sorted_nullspace).T
	
	mctq = mct * Q
	
	def getSelector2():
		selector2 = [0]
		i = 1
		while len(selector2) < 4:
			if i >= 4:
				print("nullspace Error: Ran out of independent columns")
				return None
			
			previousRank = mctq.row(selector2).rank()
			attemptedRank =  mctq.row(selector2 + [i]).rank()
			if attemptedRank > previousRank:
				selector2.append(i)
			i += 1
		return selector2
	selector2 = getSelector2()
	
	C = (mct * sp.Matrix(mfk)).row(selector2)
	D = sp.Matrix(mctq).row(selector2)
	
	DInv = D.inv()
	
	DInvC = DInv * C
	
	DInv = DInv.tolist()
	mc = mct.T.tolist()
	
	return fk, fu, mfk, mfu, selector, K, U, UInv, selector2, mc, DInv

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
fk, fu, mfk, mfu, selector, K, U, UInv, selector2, mc, DInv = getBCdata(normal)
latexCode = getLatex(f, normal, fk, fu, mfk, mfu, K, U, UInv, mc, DInv)

