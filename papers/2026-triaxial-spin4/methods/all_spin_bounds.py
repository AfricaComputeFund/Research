#!/usr/bin/env python3
"""Independent all-spin Rayleigh bounds for pure SU(2) quadratic Q_s.
This checks a conditional mathematical lemma, not any IKKT physical spectrum.
Requires numpy, scipy and sympy.
"""
import numpy as np
import sympy as sp
from scipy.linalg import eigvalsh
from scipy.optimize import brentq

s, a, b, c = sp.symbols('s a b c', positive=True)
assert sp.simplify(b*s*(s+1) - (b+b/s)*s*s) == 0
assert sp.simplify((a+b)*s/2 - ((a+b)/(2*s))*s*s) == 0

def matrices(spin):
    m = np.arange(-spin, spin + 1, dtype=float)
    p = np.zeros((2*spin+1, 2*spin+1), dtype=complex)
    for i, mi in enumerate(m[:-1]):
        p[i+1, i] = np.sqrt((spin-mi)*(spin+mi+1))
    jx = (p+p.conj().T)/2
    jy = (p-p.conj().T)/(2j)
    jz = np.diag(m)
    assert np.allclose(jx@jx+jy@jy+jz@jz, spin*(spin+1)*np.eye(len(m)), atol=2e-12)
    return jx,jy,jz

for aa,bb in [(1,1),(2,1),(5,1),(10,9)]:
    for spin in (1,2,3,4,6,8):
        jx,jy,jz = matrices(spin)
        base=aa*(jx@jx)+bb*(jy@jy)
        zz=jz@jz
        low=bb/spin
        high=(aa+bb)/(2*spin)
        assert eigvalsh(base-(1-1e-6)*low*zz)[0] > -2e-11
        top = np.zeros(2*spin+1); top[-1]=1
        assert abs(np.real(top@(base-high*zz)@top)) < 2e-11
        root=brentq(lambda x: eigvalsh(base-x*zz)[0],0,2*high+1,xtol=1e-13)
        assert low-1e-10 <= root <= high+1e-10, (aa,bb,spin,root,low,high)
        if aa==bb:
            assert abs(root-bb/spin)<2e-10
print('PASS: 2 exact symbolic identities, 24 spin-matrix Casimir controls, 24 spectral bound/threshold checks, isotropic equalities')
