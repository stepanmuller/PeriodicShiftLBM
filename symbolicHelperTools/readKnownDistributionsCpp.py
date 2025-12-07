import re

def stringToNumber(s):
	match = re.search(r"\d+", s)
	if match:
		return int(match.group())
	else:
		return None

def readKnownDistributionsCpp(i, normal, fk, fu, chosenMoments, K, UInv, meqDict):
	lines = []
	if i == 0:
		lines.append("// Reads known distributions")
		lines.append("// Outer normal identifies cell type - face cell / edge cell / corner cell and its orientation")
		lines.append("// Example: cell has outer normal [1, -1, 0]")
		lines.append("// -> Edge cell, has neighbours in both Z directions but no neighbours in +X and -Y direction")
		lines.append("")
		lines.append("if (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	else:
		lines.append("else if (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	lines.append("{")
	lines.append("	// Reading known distributions")
	for k in fk:
		number = stringToNumber(k)
		line = "	float f"
		line += str(number)
		line += " = f"
		line += str(number)
		line += "ArrayView[shiftedIndex["
		line += str(number)
		line += "]];"
		lines.append(line)
	lines.append("}")
	return lines
