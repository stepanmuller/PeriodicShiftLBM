import numpy as np

def getLatex(f, cx, cy, cz, w, normal, fk, fu, mLabels,	mfkRows, mfuRows, uniqueMfuRows, mGroups, chosenMoments, K, U, UInv, meqDict):
	velocitySet(f, cx, cy, cz, w)
	fkfu(normal, fk, fu)
	moments(normal, mLabels, mfkRows, mfuRows)
	groups(normal, uniqueMfuRows, mGroups)
	functionChosenMoments(normal, chosenMoments)
	Matrix(U, "U")
	Matrix(K, "K")
	Matrix(UInv, "U^{-1}")
	consistencyCondition(normal, mLabels, mfkRows, mfuRows)
	
	
def velocitySet(f, cx, cy, cz, w):
	lines = []
	lines.append(r"\begin{table}[H]")
	lines.append(r"\centering")
	lines.append(r"\caption{D3Q27 velocity set}")
	lines.append(r"\setlength{\tabcolsep}{1.6pt}")
	lines.append(r"\renewcommand{\arraystretch}{1.3}")
	lines.append(r"\resizebox{\textwidth}{!}{%")
	lines.append(r"\begin{tabular}{c|*{27}{c}} ")
	line = ""
	for el in f:
		line += "& $" + el + "$ "
	line += "\\\\"
	lines.append(line)
	lines.append(r"\hline")
	line = "$c_x$ "
	for el in cx:
		line += "& $" + str(el) + "$ "
	line += "\\\\"
	lines.append(line)
	line = "$c_y$ "
	for el in cy:
		line += "& $" + str(el) + "$ "
	line += "\\\\"
	lines.append(line)
	line = "$c_z$ "
	for el in cz:
		line += "& $" + str(el) + "$ "
	line += "\\\\"
	lines.append(line)
	lines.append(r"\hline")
	line = "$w$ "
	for el in w:
		num, den = el.strip().split("/")
		line += f"& $\\frac{{{num}}}{{{den}}}$ "
	line += "\\\\"
	lines.append(line)
	lines.append(r"\end{tabular}")
	lines.append(r"}")
	lines.append(r"\end{table}")
	lines.append("")
	print("\n".join(lines))

def fkfu(normal, fk, fu):
	knownString = ", ".join(f"{d}" for d in fk)
	unknownString = ", ".join(f"{d}" for d in fu)
	normalString = ", ".join(f"{d}" for d in normal)
	text = rf"""
\begin{{table}}[H]
\centering
\caption{{Known and unknown distributions for a cell with (${normalString}$) outer normal}}
\begin{{tabular}}{{l|l}}
\textbf{{Known distributions}} & $\uline{{f_k}} = ({knownString})^T$ \\
\textbf{{Unknown distributions}} & $\uline{{f_u}} = ({unknownString})^T$ \\
\end{{tabular}}
\end{{table}}

	"""
	print(text)
	
def moments(normal, mLabels, mfkRows, mfuRows):
	normalString = ", ".join(f"{d}" for d in normal)
	
	lines = []

	# Header
	lines.append(r"\begin{table}[H]")
	lines.append(r"\centering")
	lines.append(rf"\caption{{Moments for a cell with $({normalString})$ outer normal}}")
	lines.append(r"\setlength{\tabcolsep}{3pt}")
	lines.append(r"\renewcommand{\arraystretch}{1.3}")
	lines.append(r"\newcommand{\num}[1]{%")
	lines.append(r"\ifnum#1<0 #1\else \phantom{-}#1\fi")
	lines.append(r"}")
	lines.append(r"\resizebox{\textwidth}{!}{%")
	lines.append(r"\begin{tabular}{l|l|l}")
	lines.append(r"\textbf{Moment} & \textbf{Known distributions} & \textbf{Unknown distributions} \\")
	lines.append(r"\hline")

	# Each moment produces a two–row block (multirow)
	for label, knownList, unknownList in zip(mLabels, mfkRows, mfuRows):

		# Join terms into latex inline math and remove quotes
		knownStr = "("
		for el in knownList:
			knownStr += r"\num{" + str(el) + "}, "
		knownStr = knownStr[:-2]
		knownStr += ")^T"
		unknownStr = "("
		for el in unknownList:
			unknownStr += r"\num{" + str(el) + "}, "
		unknownStr = unknownStr[:-2]
		unknownStr += ")^T"
		momentLabelLatex = f"${label}$"
		lines.append(rf"{momentLabelLatex} & $\uline{{f_k}} \cdot {knownStr}$ & $ \uline{{f_u}} \cdot {unknownStr}$ \\")
		lines.append(r"\hline")

	lines.append(r"\end{tabular}")
	lines.append(r"}")
	lines.append(r"\end{table}")
	lines.append("")
	print("\n".join(lines))
	
def groups(normal, uniqueMfuRows, mGroups):
	normalString = ", ".join(f"{d}" for d in normal)
	
	lines = []

	lines.append(r"\begin{table}[H]")
	lines.append(r"\centering")
	lines.append(rf"\caption{{Unknown distribution groups for a cell with $({normalString})$ outer normal}}")
	lines.append(r"\setlength{\tabcolsep}{3pt}")
	lines.append(r"\renewcommand{\arraystretch}{1.3}")
	lines.append(r"\newcommand{\num}[1]{%")
	lines.append(r"\ifnum#1<0 #1\else \phantom{-}#1\fi")
	lines.append(r"}")
	lines.append(r"\begin{tabular}{l|l}")
	lines.append(r"\textbf{Unknown distribution group} & \textbf{Moments} \\")
	lines.append(r"\hline")

	for group, moments in zip(uniqueMfuRows, mGroups):
		momentStr = ", ".join(f"${m}$" for m in moments)
		unknownStr = "("
		for el in group:
			unknownStr += r"\num{" + str(el) + "}, "
		unknownStr = unknownStr[:-2]
		unknownStr += ")^T"
		lines.append(rf"$ \uline{{f_u}} \cdot {unknownStr}$ & {momentStr}\\")
		lines.append(r"\hline")
	lines.append(r"\end{tabular}")
	lines.append(r"\end{table}")

	print("\n".join(lines))
	
def functionChosenMoments(normal, chosenMoments):
    momentString = ", ".join(str(d) for d in chosenMoments)
    normalString = ", ".join(str(d) for d in normal)

    text = rf"""
\begin{{table}}[H]
\centering
\caption{{Chosen moments for a cell with $({normalString})$ outer normal}}
\setlength{{\tabcolsep}}{{3pt}}
\renewcommand{{\arraystretch}}{{1.3}}
\resizebox{{\textwidth}}{{!}}{{%
\begin{{tabular}}{{l|l}}
\textbf{{Chosen moments}} & $\uline{{m}} = ({momentString})^T$ \\
\end{{tabular}}
}}
\end{{table}}
"""

    print(text)

def Matrix(U, name):
    n_cols = len(U[0])
    col_format = "c" * n_cols  # 'c' for each column
    lines = [rf"\uuline{{{name}}} = ", r"\left[", rf"\begin{{array}}{{{col_format}}}"]
    for row in U:
        row_str = " & ".join(str(val) for val in row)
        lines.append(f"{row_str} \\\\")
    lines.append(r"\end{array}")
    lines.append(r"\right]")
    lines.append("")
    print("\n".join(lines))
    
def consistencyCondition(normal, mLabels, mfkRows, mfuRows):
	normalnp = np.array(normal)
	sign = np.sum(normalnp)
	normalnpabs = normalnp * normalnp
	if np.sum(normalnpabs) > 1:
		return
	which = np.where(normalnpabs >0)[0][0]
	which += 1

	normalString = ", ".join(f"{d}" for d in normal)
	
	lines = []
	
	label = mLabels[0]
	knownList = mfkRows[0]
	unknownList = mfuRows[0]

	knownStr = "("
	for el in knownList:
		knownStr += r"\num{" + str(el) + "}, "
	knownStr = knownStr[:-2]
	knownStr += ")^T"
	unknownStr = "("
	for el in unknownList:
		unknownStr += r"\num{" + str(el) + "}, "
	unknownStr = unknownStr[:-2]
	unknownStr += ")^T"
	momentLabelLatex = f"${label}$"
	lines.append(rf"From {momentLabelLatex} it yields that")
	lines.append(r"\begin{equation}")
	lines.append(r"\newcommand{\num}[1]{%")
	lines.append(r"\ifnum#1<0 #1\else \phantom{-}#1\fi")
	lines.append(r"}")
	lines.append(r"\begin{gathered}")
	lines.append(rf"\rho = \uline{{f_k}} \cdot {knownStr} \\ + \uline{{f_u}} \cdot {unknownStr}")
	lines.append(r"\end{gathered}")
	lines.append(r"\label{eq:BCConsistency1}")
	lines.append(r"\end{equation}")

	label = mLabels[which]
	knownList = mfkRows[which]
	unknownList = mfuRows[which]
	
	if which == 1:
		rhouxyz = r"\rho u_x"
	elif which == 2:
		rhouxyz = r"\rho u_y"
	elif which == 3:
		rhouxyz = r"\rho u_z"
	
	knownStr2 = "("
	for el in knownList:
		knownStr2 += r"\num{" + str(el) + "}, "
	knownStr2 = knownStr2[:-2]
	knownStr2 += ")^T"
	unknownStr2 = "("
	for el in unknownList:
		unknownStr2 += r"\num{" + str(el) + "}, "
	unknownStr2 = unknownStr2[:-2]
	unknownStr2 += ")^T"
	momentLabelLatex = f"${label}$"
	lines.append(rf"\noindent Similarly, from {momentLabelLatex} it yields that")
	lines.append(r"\begin{equation}")
	lines.append(r"\newcommand{\num}[1]{%")
	lines.append(r"\ifnum#1<0 #1\else \phantom{-}#1\fi")
	lines.append(r"}")
	lines.append(r"\begin{gathered}")
	lines.append(rf"{rhouxyz} = \uline{{f_k}} \cdot {knownStr2} \\ + \uline{{f_u}} \cdot {unknownStr2}")
	lines.append(r"\end{gathered}")
	lines.append(r"\label{eq:BCConsistency2}")
	lines.append(r"\end{equation}")
	lines.append("")
	
	finalStr = "("
	for i, el in enumerate(mfkRows[0]):
		el += sign * mfkRows[which][i]
		finalStr += r"\num{" + str(el) + "}, "
	finalStr = finalStr[:-2]
	finalStr += ")^T"
	
	lines.append(r"\begin{equation}")
	lines.append(r"\newcommand{\num}[1]{%")
	lines.append(r"\ifnum#1<0 #1\else \phantom{-}#1\fi")
	lines.append(r"}")
	lines.append(r"\begin{gathered}")
	lines.append(rf"\rho + ({sign}) {rhouxyz} \\ = \\ \uline{{f_k}} \cdot {knownStr} \\ + ({sign})\uline{{f_k}} \cdot {knownStr2} \\ = \\ \uline{{f_k}} \cdot {finalStr}")
	lines.append(r"\end{gathered}")
	lines.append(r"\end{equation}")
	lines.append("")
	
	print("\n".join(lines))
