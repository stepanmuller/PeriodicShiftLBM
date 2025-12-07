import re
import numpy as np

def stringToNumber(s):
	match = re.search(r"\d+", s)
	if match:
		return int(match.group())
	else:
		return None

def pOutletConsistencyCpp(i, normal, fk, mLabels, mfkRows):
	lines = []
	if i == 0:
		lines.append("// Known distributions must be read above")
		lines.append("// Reads prescribed float rho from local cell")
		lines.append("// Calculates normal float u_ from consistency condition")
		lines.append("// Calculates two tangential u_, u_ from moments")
		lines.append("// Works for face cells only. Outer normal points in the direction where no neighbours are")
		lines.append("")
		lines.append("const float rho = rhoArrayView[cell];")
		lines.append("")
		lines.append("if (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	else:
		lines.append("elif (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	lines.append("{")
	
	momentIndexNormal = np.where(np.array(normal) != 0)[0][0] + 1
	momentSign = np.sum(np.array(normal))
	
	line = "	const float fkProduct ="
	for j, coeff2 in enumerate(mfkRows[momentIndexNormal]):
		multiplier = 1 - coeff2
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
	lines.append("    float " + velocity + " = 1.f - fkProduct / rho;")
	
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
		lines.append("    float uy = fkTangential1 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uy
		lines.append("    float uz = fkTangential2 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uz
	elif momentIndexNormal == 2: #y is normal direction
		velocity = "uy * uy"
		lines.append("    float ux = fkTangential1 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get ux
		lines.append("    float uz = fkTangential2 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uz
	else: #z is normal direction
		velocity = "uz * uz"
		lines.append("    float ux = fkTangential1 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get ux
		lines.append("    float uy = fkTangential2 / ( (2.f/3.f - " + str(velocity) + " ) * rho );") #to get uy
	
	lines.append("}")
	return lines
