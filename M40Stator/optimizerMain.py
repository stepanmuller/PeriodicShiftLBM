import os
import math
import numpy as np
import subprocess
from matplotlib import pyplot as plt
from stl import mesh
from bladeGenerator import *

maximumIterations = 1000

-3.5; 15.625; -14.75; 0.5; 0.0; -19.0; 26.875; -41.5; 0.5; 0.0; -3.75; 24.6875; -37.5; 0.5; 0.0

parameterNames = [ "gamma0", "c0", "beta0", "btoc0", "lnw0", "gamma1", "c1", "beta1", "btoc1", "lnw1", "gamma2", "c2", "beta2", "btoc2", "lnw2" ]
parameterDefaults = [	
-3.5, 15.625, -14.75, 0.5, 0, 
-19, 26.875, -41.5, 0.5, 0, 
-3.75, 24.6875, -37.5, 0.5, 0
]
parameterEnables = [
True, True, True, False, False, 
True, True, True, False, False, 
True, True, True, False, False
]
parameterMins = [
-90, 3.5, -90, 0.1, -3, 
-90, 5.5, -90, 0.1, -3, 
-90, 7.5, -90, 0.1, -3, 
]
parameterMaxs = [
0, 21, 0, 0.9, 3, 
0, 33, 0, 0.9, 3, 
0, 45, 0, 0.9, 3, 
]

parametersFile = "optimizerParameters.txt"
resultsFile = "optimizerResults.txt"

activeIndexList = [i for i, enabled in enumerate(parameterEnables) if enabled]
activeIndexCount = len(activeIndexList)

def launchCase(caseID, parameters):
	print(f"Launching case {caseID} with parameters: {parameters}")
	
	gamma0, c0, beta0, btoc0, lnw0, gamma1, c1, beta1, btoc1, lnw1, gamma2, c2, beta2, btoc2, lnw2 = parameters
	gamma = [gamma0, gamma1, gamma2]
	c = [c0, c1, c2]
	beta = [beta0, beta1, beta2]
	b = [c0 * btoc0, c1 * btoc1, c2 * btoc2]
	w = [np.exp(lnw0), np.exp(lnw1), np.exp(lnw2)]
	parameterMatrixList = [gamma, c, beta, b, w]
	parameterMatrix = np.array(parameterMatrixList).T
	generateBlade(parameterMatrix)

	command = ["./statorMain", str(caseID)]
	subprocess.run(
		command, 
		check=True, 
		stdout=subprocess.DEVNULL, 
		stderr=subprocess.DEVNULL
	)

def readDelimitedFile(filepath):
	"""Reads a semicolon-delimited file and returns a list of data rows."""
	data = []
	if not os.path.exists(filepath):
		return data
	with open(filepath, "r") as f:
		lines = f.readlines()[1:] # Skip header
		for line in lines:
			if not line.strip(): continue
			parts = [float(x.strip()) for x in line.strip().split(";")]
			data.append(parts)
	return data

def isDuplicate(newParameters, allCases, stepSize):
	"""Checks if a parameter set has already been evaluated."""
	for case in allCases:
		oldParameters = case[1:] # Skip caseID
		isMatch = True
		for i in range(len(newParameters)):
			tolerance = 0.51 * stepSize * (parameterMaxs[i] - parameterMins[i])
			if not math.isclose(newParameters[i], oldParameters[i], abs_tol=tolerance):
				isMatch = False
				break
		if isMatch:
			return True
	return False

# Main Optimization Loop
for iteration in range(maximumIterations):
	
	# 1. Initialize parameters file if it doesn't exist
	if not os.path.exists(parametersFile):
		print(f"Creating {parametersFile}")
		with open(parametersFile, "w") as f:
			header = "caseID; " + "; ".join(parameterNames)
			f.write(header + "\n")
			
			# Write Case 0 (Defaults)
			case0 = [0] + parameterDefaults
			f.write("; ".join(map(str, case0)) + "\n")

	# 2. Initialize results file if it doesn't exist
	if not os.path.exists(resultsFile):
		with open(resultsFile, "w") as f:
			f.write("caseID; F\n")

	# 3. Read current state
	allCases = readDelimitedFile(parametersFile)
	allResults = readDelimitedFile(resultsFile)
	
	resultDict = {int(row[0]): row[1] for row in allResults}
	
	# 4. Check for missing results and launch
	missingResultFound = False
	for case in allCases:
		caseID = int(case[0])        
		if caseID not in resultDict:
			parameters = case[1:]
			launchCase(caseID, parameters)
			missingResultFound = True
			break # Break cases loop to restart outermost loop     
	if missingResultFound:
		continue # Restarts the outermost loop

	# 5. Calculate next point via Coordinate Descent
	bestCaseID = None
	bestF = float('inf')
	bestParameters = None
	for case in allCases:
		caseID = int(case[0])
		F = resultDict[caseID]
		if F < bestF:
			bestF = F
			bestParameters = case[1:]
			bestCaseID = caseID
				
	print(f"Best case ID = {bestCaseID}, (F = {bestF})")
	
	# 6. Generate Candidates
	stepSize = 0.25
	newParameters = None
	
	while newParameters is None:
		for index in activeIndexList:
			delta = stepSize * (parameterMaxs[index] - parameterMins[index])
			
			# Try MINUS direction
			candidateMinus = list(bestParameters)
			candidateMinus[index] -= delta
			if candidateMinus[index] >= parameterMins[index]:
				if not isDuplicate(candidateMinus, allCases, stepSize):
					newParameters = candidateMinus
					print(f"New parameters found at step size {stepSize}, changing {parameterNames[index]} from {bestParameters[index]} to {newParameters[index]}")
					break
					
			# Try PLUS direction
			candidatePlus = list(bestParameters)
			candidatePlus[index] += delta
			if candidatePlus[index] <= parameterMaxs[index]:
				if not isDuplicate(candidatePlus, allCases, stepSize):
					newParameters = candidatePlus
					print(f"New parameters found at step size {stepSize}, changing {parameterNames[index]} from {bestParameters[index]} to {newParameters[index]}")
					break
					
		# If no candidates were valid/unique at this scale, halve it
		if newParameters is None:
			stepSize *= 0.5

	# 7. Write and launch new case		
	newCaseID = int(allCases[-1][0]) + 1
	
	with open(parametersFile, "a") as f:
		caseLine = [newCaseID] + newParameters
		f.write("; ".join(map(str, caseLine)) + "\n")
		
	launchCase(newCaseID, newParameters)

print("Can't believe it, but we reached the maximumIterations limit. Terminating.")
