import re
import numpy as np

def restoreRho(index, normal, fk, SCQrho, SCMfk, SCQu):
	lines = []
	if index == 0:
		lines.append("__host__ __device__ void restoreRho(")
		lines.append("	const int &outerNormalX, const int &outerNormalY, const int &outerNormalZ,")
		lines.append("	const float &rho, float &ux, float &uy, float &uz,")
		lines.append("	const float (&f)[27]")
		lines.append(")")
		lines.append("{")
		lines.append("	if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
	else:
		lines.append("	else if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
	lines.append("	{")
	
	directions = ["x", "y", "z"]
	line = "		const float scqu ="
	for j, multiplier in enumerate(SCQu[0]):
		if multiplier == "0":
			continue
		elif multiplier == "1":
			line += " + u" + directions[j]
			continue
		elif multiplier == "-1":
			line += " - u" + directions[j]
			continue
		else:
			line += " + "
			line += multiplier
			line += " * u"
			line += directions[j]
	if line[-1] == "=":
		line += " 0.f"
	line += ";"
	lines.append(line)
		
	line = "		const float s = " + SCQrho[0][0] + " + scqu;"
	lines.append(line)
		
	line = "		const float scmf ="
	for j, multiplier in enumerate(SCMfk[0]):
		if multiplier == "0":
			continue
		elif multiplier == "1":
			line += " + " + fk[j]
			continue
		elif multiplier == "-1":
			line += " - " + fk[j]
			continue
		else:
			line += " + "
			line += multiplier
			line += " * "
			line += fk[j]
	if line[-1] == "=":
		line += " 0.f"
	line += ";"
	lines.append(line)
		
	line = "		rho = scmf / s;"
	lines.append(line)
	
	lines.append("	}")
	
	return lines
