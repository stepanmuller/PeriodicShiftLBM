import re
import numpy as np

def stringToNumber(s):
	match = re.search(r"\d+", s)
	if match:
		return int(match.group())
	else:
		return None

def applyMBBC(i, normal, fk, fu, mfkRows, mfuRows, chosenMoments, K, U, UInv, meqDict):
	lines = []
	if i == 0:
		lines.append("// Applies moment based boundary condition")
		lines.append("// Calculates unknown distributions from rho, ux, uy, uz")
		lines.append("")
		lines.append("// Source paper: Pavel Eichler, Radek Fucik, and Pavel Strachota.")
		lines.append("// Investigation of mesoscopic boundary conditions for lattice boltzmann method in laminar flow problems.")
		lines.append("// Computers & Mathematics with Applications, 173:87–101, 2024.")
		lines.append("")
		lines.append("if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
	else:
		lines.append("else if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
	lines.append("{")
	lines.append("	// Multiply K fk")
	for i, row in enumerate(K):
		line = "	const float kf" + str(i) + " ="
		for j, multiplier in enumerate(row):
			if multiplier == 0:
				continue
			if multiplier == 1:
				number = stringToNumber(fk[j])
				line += " + f" + str(number)
			elif multiplier == -1:
				number = stringToNumber(fk[j])
				line += " - f" + str(number)
		if line[-1] == "=":
			line += " 0.f;"
		else:
			line += ";"
		lines.append(line)
	lines.append("	// Calculate equilibrium moments")
	for i, moment in enumerate(chosenMoments):
		line = "	const float m" + str(i) + " = "
		line += meqDict[moment]
		line += ";"
		lines.append(line)
	lines.append("	// Subtract m - Kfk")
	for i, moment in enumerate(chosenMoments):
		line = "	const float s" + str(i) + " = "
		line += "m" + str(i) + " - kf" + str(i) + ";"
		lines.append(line)
	lines.append("	// Multiply U^-1 * (m - Kfk) to get unknown distributions")
	for i, unknown in enumerate(fu):
		number = stringToNumber(unknown)
		line = "	f" + str(number) + " ="
		for j, multiplier in enumerate(UInv[i]):
			if multiplier == 0:
				continue
			elif multiplier == 1:
				line += " + s" + str(j)
				continue
			elif multiplier == -1:
				line += " - s" + str(j)
				continue
			elif multiplier > 0:
				line += " + "
			elif multiplier < 0:
				line += " - "
			multiplier = abs(multiplier)
			multiplier = str(float(multiplier)) + "f"
			line += multiplier
			line += " * "
			line += "s" + str(j)
		if line[-1] == "=":
			line += " 0.f"
		line += ";"
		lines.append(line)
	lines.append("}")

	return lines
