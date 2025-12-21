import sympy as sp
from math import gcd
from functools import reduce

# ---- helper: normalize integer vectors ----
def normalize_integer_vector(v):
    nums = [int(x) for x in v if x != 0]
    if not nums:
        return v
    g = reduce(gcd, nums)
    return [int(x // g) for x in v]

# ---- helper: support of a vector ----
def support(v):
    return {i for i, x in enumerate(v) if x != 0}

# ---- your b-vectors (rows) ----
B = sp.Matrix([
    [ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1],
    [-1,  0,  0,  1, -1, -1, -1,  0,  0, -1,  1,  0, -1,  1, -1,  1,  1, -1, -1],
    [ 0,  0, -1,  0,  0,  0, -1,  1, -1,  1, -1, -1,  1, -1, -1,  1, -1,  1, -1],
    [ 0, -1,  0, -1,  1, -1,  0, -1,  1,  0,  0, -1, -1,  1,  1, -1, -1,  1, -1],
    [ 1,  0,  0,  1,  1,  1,  1,  0,  0,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1],
    [ 0,  0,  1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1],
    [ 0,  1,  0,  1,  1,  1,  0,  1,  1,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1],
    [ 0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  1, -1, -1, -1, -1,  1,  1,  1],
    [ 0,  0,  0, -1, -1,  1,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1],
    [ 0,  0,  0,  0,  0,  0,  1,  0,  0, -1, -1,  0, -1, -1,  1,  1, -1, -1,  1],
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1, -1],
])

# ---- compute nullspace of B^T ----
nullspace = B.T.nullspace()

print(f"Nullspace dimension: {len(nullspace)}\n")

# ---- collect circuits (minimal supports) ----
circuits = []

for v in nullspace:
    # clear denominators
    v = sp.lcm([x.q for x in v]) * v
    v = normalize_integer_vector(v)
    supp = support(v)

    # check minimality
    is_minimal = True
    for c in circuits:
        if supp > c["support"]:
            is_minimal = False
            break

    if is_minimal:
        circuits.append({
            "coeffs": v,
            "support": supp
        })

# ---- print results ----
for i, c in enumerate(circuits, 1):
    print(f"Circuit {i}:")
    print("  rows involved:", sorted(j+1 for j in c["support"]))
    print("  coefficients :", c["coeffs"])
    print()
