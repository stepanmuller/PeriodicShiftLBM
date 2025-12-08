import re
import numpy as np

def stringToNumber(s):
	match = re.search(r"\d+", s)
	if match:
		return int(match.group())
	else:
		return None

def applyVelocityInlet(i, normal, fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict):
	lines = []
	flag = 1000 + (100 * (normal[0] + 5) + 10 * (normal[1] + 5) + (normal[2] + 5))
	if i == 0:
		lines.append("// Applies velocity inlet moment based boundary condition on a face cell")
		lines.append("// 1. Reads known distributions and prescribed ux, uy, uz")
		lines.append("// 2. Calculates rho from consistency condition")
		lines.append("// 3. Applies general moment based boundary condition to find unknown distributions")
		lines.append("")
		lines.append("// Source paper: Pavel Eichler, Radek Fucik, and Pavel Strachota.")
		lines.append("// Investigation of mesoscopic boundary conditions for lattice boltzmann method in laminar flow problems.")
		lines.append("// Computers & Mathematics with Applications, 173:87–101, 2024.")
		lines.append("")
		lines.append("// Reading prescribed velocity inlet ux, uy, uz")
		lines.append("float ux = uxArrayView[cell];")
		lines.append("float uy = uyArrayView[cell];")
		lines.append("float uz = uzArrayView[cell];")
		lines.append("// Declaring rho")
		lines.append("float rho = 1.f;")
		lines.append("")
		lines.append("if (flag == " + str(flag) + ") // outer normal " + str(normal))
	else:
		lines.append("else if (flag == " + str(flag) + ") // outer normal " + str(normal))
	lines.append("{")
	
	lines.append("	// Reading known distributions fk")
	for k in fk:
		number = stringToNumber(k)
		line = "	f"
		line += str(number)
		line += " = f"
		line += str(number)
		line += "ArrayView[shiftedIndex["
		line += str(number)
		line += "]];"
		lines.append(line)
		
	lines.append("	// Applying consistency condition to find rho")
	momentIndex = np.where(np.array(normal) != 0)[0][0] + 1
	momentSign = np.sum(np.array(normal))
	line = "	const float fkProduct ="
	for j, coeff2 in enumerate(mfkRows[momentIndex]):
		multiplier = 1 + momentSign * coeff2
		if multiplier == 0:
			continue
		if multiplier == 1:
			number = stringToNumber(fk[j])
			line += " + f" + str(number)
		elif multiplier == -1:
			number = stringToNumber(fk[j])
			line += " - f" + str(number)
		elif multiplier > 0:
			number = stringToNumber(fk[j])
			line += " + " + str(int(multiplier)) + ".f * f" + str(number)
		elif multiplier < 0:
			number = stringToNumber(fk[j])
			line += " - " + str(int(abs(multiplier))) + ".f * f" + str(number)
	line += ";"
	lines.append(line)
	velocity = "ux"
	if momentIndex == 2:
		velocity = "uy"
	elif momentIndex == 3:
		velocity = "uz"
	if momentSign < 0:
		lines.append("    rho = fkProduct / (1.f - " + velocity + ");")
	else:
		lines.append("    rho = fkProduct / (1.f + " + velocity + ");")
	
	lines.append("	// At this point rho, ux, uy, uz are known")
	lines.append("	// Applying general MBBC")
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
