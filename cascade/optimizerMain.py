import os
import math
import numpy as np
import subprocess
from matplotlib import pyplot as plt
from stl import mesh
from vaneGenerator import *

parameterNames = ["h", "a", "b", "l", "w", "t", "s"]
parameterDefaultValues = [10.0, 10.0, 5.0, 15.0, 1.0, 1.0, 15.0]
parameterEnables = [False, True, True, True, True, False, True]
parameterDeltas = [1.0, 10.0, 3.0, 10.0, 5.0, 1.0, 15.0]

# --- ADDED: Bounds for the Coordinate Descent ---
parameterMinimums = [10.0, 1.0, 1.0, 1.0, 0.01, 1.0, 5.0]
parameterMaximums = [10.0, 60.0, 9.0, 80.0, 20.0, 1.0, 80.0]

maximumIterations = 100
parametersFile = "optimizerParameters.txt"
resultsFile = "optimizerResults.txt"

activeIndexList = [i for i, enabled in enumerate(parameterEnables) if enabled]
activeIndexCount = len(activeIndexList)

# --- Helper Functions ---
def launchCase(caseID, parameters):
    # Fixed 'params' to 'parameters' in the print statement
    print(f"Launching case {caseID} with params: {parameters}")
    generateVane(parameters)
    h = parameters[0]
    s = parameters[6]
    command = ["./cascadeMain", str(caseID), str(h), str(s)]
    #subprocess.run(command, check=True)
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

def isDuplicate(newParameters, allCases, tolerance=1e-5):
    """Checks if a parameter set has already been evaluated."""
    for caseData in allCases:
        oldParams = caseData[1:] # Skip caseID
        isMatch = True
        for i in activeIndexList:
            if not math.isclose(newParameters[i], oldParams[i], abs_tol=tolerance):
                isMatch = False
                break
        if isMatch:
            return True
    return False

# --- Main Optimization Loop ---
for iteration in range(maximumIterations):
    print(f"\n--- Iteration {iteration} ---")
    
    # 1. Initialize parameters file if it doesn't exist
    if not os.path.exists(parametersFile):
        print(f"Creating {parametersFile} and initializing simplex...")
        with open(parametersFile, "w") as f:
            header = "caseID; " + "; ".join(parameterNames)
            f.write(header + "\n")
            
            # Write Case 0 (Defaults)
            case0 = [0] + parameterDefaultValues
            f.write("; ".join(map(str, case0)) + "\n")
            
            # Write Cases 1 to N (One altered parameter each)
            currentCaseID = 1
            for activeIdx in activeIndexList:
                alteredParams = list(parameterDefaultValues)
                alteredParams[activeIdx] += parameterDeltas[activeIdx]
                caseLine = [currentCaseID] + alteredParams
                f.write("; ".join(map(str, caseLine)) + "\n")
                currentCaseID += 1

    # 2. Initialize results file if it doesn't exist
    if not os.path.exists(resultsFile):
        with open(resultsFile, "w") as f:
            f.write("caseID; p\n")

    # 3. Read current state
    allCases = readDelimitedFile(parametersFile)
    allResults = readDelimitedFile(resultsFile)
    
    resultDict = {int(row[0]): row[1] for row in allResults}
    
    # 4. Check for missing results and launch
    missingResultFound = False
    for caseData in allCases:
        caseID = int(caseData[0])
        params = caseData[1:]
        
        if caseID not in resultDict:
            launchCase(caseID, params)
            missingResultFound = True
            break # Break cases loop to restart outermost loop
            
    if missingResultFound:
        continue # Restarts the outermost loop

    # 5. Calculate next point via Coordinate Descent
    print("All cases evaluated. Finding next point...")
    
    bestCaseID = None
    bestP = float('inf')
    bestParams = None
    
    for caseData in allCases:
        cID = int(caseData[0])
        if cID in resultDict:
            p = resultDict[cID]
            if p < bestP:
                bestP = p
                bestParams = caseData[1:]
                bestCaseID = cID
                
    print(f"Anchoring to best case so far: Case {bestCaseID} (p = {bestP})")
    
    # 6. Generate Candidates
    changeScale = 1.0
    minChangeScale = 1e-4
    newParameters = None
    
    while changeScale > minChangeScale and newParameters is None:
        for activeIdx in activeIndexList:
            
            # Try MINUS direction
            candidateMinus = list(bestParams)
            candidateMinus[activeIdx] -= (parameterDeltas[activeIdx] * changeScale)
            
            if candidateMinus[activeIdx] >= parameterMinimums[activeIdx]:
                if not isDuplicate(candidateMinus, allCases):
                    newParameters = candidateMinus
                    print(f"Selected MINUS step on '{parameterNames[activeIdx]}' at scale {changeScale}")
                    break
                    
            # Try PLUS direction
            candidatePlus = list(bestParams)
            candidatePlus[activeIdx] += (parameterDeltas[activeIdx] * changeScale)
            
            if candidatePlus[activeIdx] <= parameterMaximums[activeIdx]:
                if not isDuplicate(candidatePlus, allCases):
                    newParameters = candidatePlus
                    print(f"Selected PLUS step on '{parameterNames[activeIdx]}' at scale {changeScale}")
                    break
                    
        # If no candidates were valid/unique at this scale, halve it
        if newParameters is None:
            print(f"No unique valid moves found at scale {changeScale}. Halving scale...")
            changeScale *= 0.5

    # 7. Write and launch new case
    if newParameters is None:
        print("Optimization converged! No new unique points can be generated.")
        break # Exit the maximumIterations loop entirely
        
    newCaseID = int(allCases[-1][0]) + 1
    
    with open(parametersFile, "a") as f:
        caseLine = [newCaseID] + newParameters
        f.write("; ".join(map(str, caseLine)) + "\n")
        
    launchCase(newCaseID, newParameters)
    
print("Optimization loop finished.")
