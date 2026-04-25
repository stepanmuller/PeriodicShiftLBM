import os
import math
import numpy as np
import subprocess
from matplotlib import pyplot as plt
from stl import mesh
from bladeGenerator import *

maximumIterations = 10000

-3.5; 15.625; -14.75; 0.5; 0.0; -19.0; 26.875; -41.5; 0.5; 0.0; -3.75; 24.6875; -37.5; 0.5; 0.0

parameterNames = [ "gamma0", "c0", "beta0", "btoc0", "lnw0", "gamma1", "c1", "beta1", "btoc1", "lnw1", "gamma2", "c2", "beta2", "btoc2", "lnw2" ]
parameterDefaults = [	
-3.5, 15.625, -14.75, 0.5, 0, 
-19, 26.875, -41.5, 0.5, 0, 
-3.75, 24.6875, -37.5, 0.5, 0
]
parameterEnables = [
True, True, True, True, False, 
True, True, True, True, False, 
True, True, True, True, False
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

def getSmartCandidate(allCases, resultDict, bestParameters, activeIndexList, parameterMins, parameterMaxs, parameterNames):
    stepSize = 0.125
    
    while True:
        candidates = []
        
        # 1. Generate and score all possible candidates for the current stepSize
        for index in activeIndexList:
            delta = stepSize * (parameterMaxs[index] - parameterMins[index])
            
            for direction in [-1, 1]: # -1 for MINUS, 1 for PLUS
                candidate = list(bestParameters)
                candidate[index] += direction * delta
                
                # Check bounds
                if candidate[index] < parameterMins[index] or candidate[index] > parameterMaxs[index]:
                    continue
                    
                # Check uniqueness (assuming your isDuplicate function exists)
                if isDuplicate(candidate, allCases, stepSize):
                    continue
                
                # 2. Evaluate historical success for this specific parameter change
                pCandidate = candidate[index]
                threshold = 0.51 * delta
                
                historical_Fs = []
                for past_case in allCases:
                    past_caseID = int(past_case[0])
                    # Ensure the case actually finished and has an F value
                    if past_caseID not in resultDict:
                        continue 
                        
                    past_p = float(past_case[1+index])
                    
                    # If the past record had this parameter near pCandidate
                    if abs(past_p - pCandidate) <= threshold:
                        historical_Fs.append(float(resultDict[past_caseID]))
                
                # 3. Assign priority score
                if not historical_Fs:
                    # Priority 0: No historical data at all -> Test this first
                    score = (0, 0.0) 
                else:
                    # Priority 1: Has history -> Rank by the best (minimum) F found
                    score = (1, min(historical_Fs))
                    
                candidates.append((score, candidate, index, direction))
        
        # 4. Sort and select the best candidate
        if candidates:
            # Python sorts tuples element-by-element. 
            # It will group by Priority (0 then 1), then sort by the lowest F within Priority 1.
            candidates.sort(key=lambda x: x[0])
            
            best_score, best_candidate, best_index, best_dir = candidates[0]
            
            dir_str = "PLUS" if best_dir == 1 else "MINUS"
            history_str = "Unexplored" if best_score[0] == 0 else f"Explored (Best past F: {best_score[1]})"
            
            print(f"New parameters found at step size {stepSize}, changing {parameterNames[best_index]} ({dir_str}) to {best_candidate[best_index]}.")
            print(f" -> Reason: {history_str}")
            
            return best_candidate
            
        # 5. If NO candidates were valid/unique across all active indices, halve the scale
        stepSize *= 0.5
        
        # Optional: Safety break to prevent infinite loops if the space is completely exhausted
        if stepSize < 1e-6:
            print("Step size diminished to zero. Optimization exhausted.")
            return None

# Main Optimization Loop
for iteration in range(maximumIterations):
	print()
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
	
	# 6. Generate Candidates (Smart Coordinate Descent)
	newParameters = getSmartCandidate(
		allCases=allCases, 
		resultDict=resultDict, 
		bestParameters=bestParameters, 
		activeIndexList=activeIndexList, 
		parameterMins=parameterMins, 
		parameterMaxs=parameterMaxs, 
		parameterNames=parameterNames
	)

	# 7. Write and launch new case		
	newCaseID = int(allCases[-1][0]) + 1
	
	with open(parametersFile, "a") as f:
		caseLine = [newCaseID] + newParameters
		f.write("; ".join(map(str, caseLine)) + "\n")
		
	launchCase(newCaseID, newParameters)

print("Can't believe it, but we reached the maximumIterations limit. Terminating.")
