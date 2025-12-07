import re
import numpy as np

def stringToNumber(s):
	match = re.search(r"\d+", s)
	if match:
		return int(match.group())
	else:
		return None

def uInletConsistencyCpp(i, normal, fk, mfkRows):
	lines = []
	if i == 0:
		lines.append("// Known distributions must be read above")
		lines.append("// Reads prescribed floats ux, uy, uz from local cell")
		lines.append("// Calculates float rho from consistency condition")
		lines.append("// Works for face cells only. Outer normal points in the direction where no neighbours are")
		lines.append("")
		lines.append("float ux = uxArrayView[cell];")
		lines.append("float uy = uyArrayView[cell];")
		lines.append("float uz = uzArrayView[cell];")
		lines.append("")
		lines.append("if (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	else:
		lines.append("else if (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	lines.append("{")
	
	momentIndex = np.where(np.array(normal) != 0)[0][0] + 1
	momentSign = np.sum(np.array(normal))
	
	line = "	const float fkProduct ="
	for j, coeff2 in enumerate(mfkRows[momentIndex]):
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
	if momentIndex == 2:
		velocity = "uy"
	elif momentIndex == 3:
		velocity = "uz"
	lines.append("    const float rho = fkProduct / (1.f - " + velocity + ");")
	lines.append("}")
	return lines
