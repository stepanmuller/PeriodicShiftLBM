import numpy as np
import sympy as sp

def getLatex(f, normal, fk, fu, mfk, mfu, S, Sq, C, Q):
	printFkfu(normal, fk, fu)
	printMatrix(mfk.tolist(), rf"\uuline{{M_{{fk}}}}")
	printMatrix(mfu.tolist(), rf"\uuline{{M_{{fu}}}}")
	printMatrix(S.tolist(), rf"\uuline{{S}}")
	printMatrix((S*mfu).inv().tolist(), rf"\left( \uuline{{S}} \ \uuline{{M_{{fu}}}} \right)^{{-1}}")
	printMatrix(Sq.tolist(), rf"\uuline{{S_q}}")
	printMatrix((Sq*C.T*Q).inv().tolist(), rf"\left( \uuline{{S_q}} \uuline{{C^T}} \uuline{{Q}} \right)^{{-1}}")

def printFkfu(normal, fk, fu):
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

def printMatrix(U, name):
    n_cols = len(U[0])
    col_format = "c" * n_cols  # 'c' for each column
    lines = [rf"{name} = ", r"\left[", rf"\begin{{array}}{{{col_format}}}"]
    for row in U:
        row_str = " & ".join(str(val) for val in row)
        lines.append(f"{row_str} \\\\")
    lines.append(r"\end{array}")
    lines.append(r"\right]")
    lines.append("")
    print("\n".join(lines))
