import re

def stringToNumber(s):
	match = re.search(r"\d+", s)
	if match:
		return int(match.group())
	else:
		return None

def periodicBCCpp(i, normal, fu):
	lines = []
	if i == 0:
		lines.append("// Periodic BC, can also work as zero gradient depending how source cell is chosen")
		lines.append("// Reads unknown distributions from source cell")
		lines.append("// Writes them into local cell")
		lines.append("// Outer normal identifies cell type - face cell / edge cell / corner cell and its orientation")
		lines.append("// Example: cell has outer normal [1, -1, 0]")
		lines.append("// -> Edge cell, has neighbours in both Z directions but no neighbours in +X and -Y direction")
		lines.append("// In this case it is only used to identify unknown distributions and might not match the geometrical meaning")
		lines.append("")
		lines.append("if (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	else:
		lines.append("elif (outerNormalX == " + str(normal[0]) + " && outerNormalY == " + str(normal[1]) + " && outerNormalZ == " + str(normal[2]) + ")")
	lines.append("{")
	lines.append("	// Reading unknown distributions from source cell")
	for k in fu:
		number = stringToNumber(k)
		line = "	const float f"
		line += str(number)
		line += " = f"
		line += str(number)
		line += "ArrayView[shiftedSourceIndex["
		line += str(number)
		line += "]];"
		lines.append(line)
	lines.append("	// Writing them into local cell")
	for k in fu:
		number = stringToNumber(k)
		line = "	ArrayView[shiftedIndex["
		line += str(number)
		line += "]] = f"
		line += str(number)
		line += ";"
		lines.append(line)
	lines.append("}")
	return lines
