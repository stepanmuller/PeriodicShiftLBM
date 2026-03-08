import sys
import numpy as np
from matplotlib import pyplot as plt
import copy
from stl import mesh
from STLVisualizer import show_stl
from vaneGenerator import *
	
parameterNames = ["h", "a", "b", "l", "w", "t", "s"]
parameterDefaultValues = [10, 10, 5, 15, 1, 1, 15]
parameterEnables = [False, True, True, True, True, False, True]

while True:
	1) load optimizerParameters.txt, create it if it does not exist and add header "caseID; h; a; b; l; w; t; s"
	2) load optimizerResults.txt, create it if it does not exist and add header "caseID; p"
	3) for every case in optimizerParameters.txt, look if there is a matching record in optimizerResults. if not, launch the case and continue
	4) If there is a matching record for all cases, browse them and find if there is at least 2 different values for every parameter. 
	
