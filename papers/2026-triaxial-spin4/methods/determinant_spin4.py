#!/usr/bin/env python3
"""Independent spin-4 pure-Q determinant and instability root check.

Synthetic finite matrix only; NOT validation of a physical background.
Needs numpy and sympy. No previous OMF code imported.
"""
import math
import numpy as np


def spin4_operator(a, b, c):
    if not (a >= b > 0 and c >= 0):
        raise ValueError("need a >= b > 0 and c >= 0")
    s = 4
    jp = np.zeros((2*s+1, 2*s+1), dtype=complex)
    for i, m in enumerate(range(-s, s)):
        jp[i+1, i] = math.sqrt((s-m)*(s+m+1))
    jm = jp.conj().T
    j1 = (jp+jm)/2
    j2 = (jp-jm)/(2j)
    j3 = np.diag(np.arange(-s, s+1)).astype(complex)
    return a*(j1@j1) + b*(j2@j2) - c*(j3@j3)


def factor_polynomials(a, b, c):
    f1 = 2*(a+b)*c*c-(2*a*a+13*a*b+2*b*b)*c+2*a*b*(a+b)
    f2 = 64*c*c-136*(a+b)*c+(9*a*a+46*a*b+9*b*b)
    f3 = 9*c*c-(46*a+136*b)*c+(9*a*a+136*a*b+64*b*b)
    f4 = 9*c*c-(136*a+46*b)*c+(64*a*a+136*a*b+9*b*b)
    return f1, f2, f3, f4


def first_root(a, b):
    B = 2*a*a + 13*a*b + 2*b*b
    A = 2*(a+b)
    C = 2*a*b*(a+b)
    D = B*B-4*A*C
    return 2*C/(B+math.sqrt(D))


def numerical_checks():
    records = []
    for a, b in [(1,1), (2,1), (3,2), (5,1), (10,9)]:
        c4 = first_root(a,b)
        for c in [0, c4*.5, c4*.9, c4*1.1, max(1,2*c4)]:
            m = spin4_operator(a,b,c)
            assert np.linalg.norm(m-m.conj().T) < 1e-12
            direct = float(np.linalg.det(m).real)
            factors = 320*math.prod(factor_polynomials(a,b,c))
            assert math.isclose(direct, factors, rel_tol=2e-9, abs_tol=2e-4), (a,b,c,direct,factors)
        below = float(np.linalg.eigvalsh(spin4_operator(a,b,c4*.9))[0])
        at = float(np.linalg.eigvalsh(spin4_operator(a,b,c4))[0])
        above = float(np.linalg.eigvalsh(spin4_operator(a,b,c4*1.1))[0])
        assert below > 0 and abs(at) < 1e-10 and above < 0, (a,b,c4,below,at,above)
        records.append((a,b,c4,below,at,above))
    assert math.isclose(first_root(1,1), .25, abs_tol=1e-12)
    print("PASS: 25 determinant checks and 5 independent eigenvalue crossings")
    for row in records: print(row)


def symbolic_check():
    """Exact independent symbolic 9x9 determinant, not numeric interpolation."""
    import sympy as sp
    a,b,c = sp.symbols("a b c", positive=True)
    spin = 4
    jp = sp.zeros(2*spin+1)
    for col, m in enumerate(range(-spin,spin)):
        jp[col+1,col] = sp.sqrt((spin-m)*(spin+m+1))
    jm = jp.T
    j1 = (jp+jm)/2
    j2 = (jp-jm)/(2*sp.I)
    j3 = sp.diag(*range(-spin,spin+1))
    q = a*j1**2+b*j2**2-c*j3**2
    actual = sp.factor(q.det())
    claimed = 320*sp.prod(factor_polynomials(a,b,c))
    assert sp.simplify(actual-claimed) == 0
    print("PASS: exact symbolic det(Q4) = 320 F1 F2 F3 F4")


if __name__ == "__main__":
    numerical_checks()
    symbolic_check()
