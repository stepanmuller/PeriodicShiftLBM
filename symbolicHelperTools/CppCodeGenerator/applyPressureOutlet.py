import re
import numpy as np

def stringToNumber(s):
	match = re.search(r"\d+", s)
	if match:
		return int(match.group())
	else:
		return None

def applyPressureOutlet(i, normal, fk, fu, mLabels, mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict):
	lines = []
	if i == 0:
		lines.append("// Applies pressure outlet moment based boundary condition on a face cell")
		lines.append("// 1. Reads known distributions and prescribed rho")
		lines.append("// 2. Calculates normal u component from consistency condition")
		lines.append("// 3. Uses moments to calculate tangential u components")
		lines.append("// 4. Applies general moment based boundary condition to find unknown distributions")
		lines.append("")
		lines.append("// Source paper: Pavel Eichler, Radek Fucik, and Pavel Strachota.")
		lines.append("// Investigation of mesoscopic boundary conditions for lattice boltzmann method in laminar flow problems.")
		lines.append("// Computers & Mathematics with Applications, 173:87–101, 2024.")
		lines.append("")
		lines.append("// Reading prescribed pressure outlet rho")
		lines.append("const float rho = rhoArrayView[cell];")
		lines.append("// Declaring ux, uy, uz")
		lines.append("float ux = 0.f;")
		lines.append("float uy = 0.f;")
		lines.append("float uz = 0.f;")
		lines.append("")
		lines.append("if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
	else:
		lines.append("else if ( outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + " )")
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
		
	lines.append("	// Applying consistency condition to find normal u component")
	momentIndexNormal = np.where(np.array(normal) != 0)[0][0] + 1
	momentSign = np.sum(np.array(normal))
	
	line = "	const float fkProduct ="
	for j, coeff2 in enumerate(mfkRows[momentIndexNormal]):
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
	if momentIndexNormal == 2:
		velocity = "uy"
	elif momentIndexNormal == 3:
		velocity = "uz"
	if momentSign < 0:
		lines.append("    " + velocity + " = 1.f - fkProduct / rho;")
	else:
		lines.append("    " + velocity + " = - 1.f + fkProduct / rho;")
		
	lines.append("	// Using moments to find tangential u components in a local way")
	if momentIndexNormal == 1: #x is normal direction
		indexLin1 = mLabels.index("m_{010}") #to get uy
		indexCub1 = mLabels.index("m_{210}")
		indexLin2 = mLabels.index("m_{001}") #to get uz
		indexCub2 = mLabels.index("m_{201}")
	elif momentIndexNormal == 2: #y is normal direction
		indexLin1 = mLabels.index("m_{100}") #to get ux
		indexCub1 = mLabels.index("m_{120}")
		indexLin2 = mLabels.index("m_{001}") #to get uz
		indexCub2 = mLabels.index("m_{021}")
	else: #z is normal direction
		indexLin1 = mLabels.index("m_{100}") #to get ux
		indexCub1 = mLabels.index("m_{102}")
		indexLin2 = mLabels.index("m_{010}") #to get uy
		indexCub2 = mLabels.index("m_{012}")
	
	line = "	const float fkTangential1 ="
	for j in range(len(fk)):
		multiplier = mfkRows[indexLin1][j] - mfkRows[indexCub1][j]
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
	
	line = "	const float fkTangential2 ="
	for j in range(len(fk)):
		multiplier = mfkRows[indexLin2][j] - mfkRows[indexCub2][j]
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
	
	if momentIndexNormal == 1: #x is normal direction
		velocity = "ux * ux"
		lines.append("    uy = fkTangential1 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uy
		lines.append("    uz = fkTangential2 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uz
	elif momentIndexNormal == 2: #y is normal direction
		velocity = "uy * uy"
		lines.append("    ux = fkTangential1 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get ux
		lines.append("    uz = fkTangential2 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uz
	else: #z is normal direction
		velocity = "uz * uz"
		lines.append("    ux = fkTangential1 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get ux
		lines.append("    uy = fkTangential2 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uy
	
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
