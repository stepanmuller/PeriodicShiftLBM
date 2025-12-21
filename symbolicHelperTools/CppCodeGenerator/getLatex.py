import numpy as np

def getLatex(f, normal, fk, fu, mfk, mfu, K, U, UInv, mc, DInv):
	fkfu(normal, fk, fu)
	Matrix(mfk, "M_{fk}")
	Matrix(mfu, "M_{fu}")
	Matrix(U, "U")
	Matrix(K, "K")
	Matrix(UInv, "U^{-1}")
	Matrix(mc, "M_{c}")
	Matrix(DInv, "D^{-1}")

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
