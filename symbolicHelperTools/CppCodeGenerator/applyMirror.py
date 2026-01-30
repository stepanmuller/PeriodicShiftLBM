import re
import numpy as np

def applyMirror(index, normal, fu, fuIndexList, cx, cy, cz):
	lines = []
	if index == 0:
		lines.append("__host__ __device__ void applyMirror(")
		lines.append("	const int &outerNormalX, const int &outerNormalY, const int &outerNormalZ,")
		lines.append("	float (&f)[27]")
		lines.append(")")
		lines.append("{")
		lines.append("	if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
	else:
		lines.append("	else if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
	lines.append("	{")
	for i, fuIndex in enumerate(fuIndexList):
		myCx = cx[fuIndex]
		if myCx == -normal[0]:
			myCx = -myCx
		myCy = cy[fuIndex]
		if myCy == -normal[1]:
			myCy = -myCy
		myCz = cz[fuIndex]
		if myCz == -normal[2]:
			myCz = -myCz
		result = 0
		for j in range(len(cx)):
			if cx[j] == myCx and cy[j] == myCy and cz[j] == myCz:
				result = j
				break
		lines.append("	" + fu[i] + " = f[" + str(result) + "];" )
	lines.append("	}")
	
	return lines
