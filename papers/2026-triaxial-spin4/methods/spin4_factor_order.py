#!/usr/bin/env python3
"""Exact global spin-4 determinant-factor root ordering.

Source-conditional algebra for the pure-Q spin-four operator; no claim
about lower-spin ordering or the physical fluctuation operator.
Requires SymPy.
"""
import sympy as sp

r, x = sp.symbols("r x", real=True)
F1 = 2*(r+1)*x**2-(2*r**2+13*r+2)*x+2*r*(r+1)
F2 = 64*x**2-136*(r+1)*x+9*r**2+46*r+9
F3 = 9*x**2-(46*r+136)*x+9*r**2+136*r+64
F4 = 9*x**2-(136*r+46)*x+64*r**2+136*r+9
F = [F1,F2,F3,F4]
discriminants = [
    4*r**4+36*r**3+145*r**2+36*r+4,
    64*(253*r**2+394*r+253),
    64*(28*r**2+119*r+253),
    64*(253*r**2+119*r+28),
]
resultants = [
    -2268*(r-1)**6,
    -63*(r+4)**4*(4*r+1)**2,
    -63*(r+4)**2*(4*r+1)**4,
]
for f,expected in zip(F,discriminants):
    assert sp.expand(sp.discriminant(f,x)-expected)==0
    A,B,C = sp.Poly(f,x).all_coeffs()
    # For every r>=1: A>0, B<0, C>0, discriminant>0.
    # Thus both roots are distinct, real and strictly positive.
    assert A.subs(r,1)>0 and B.subs(r,1)<0 and C.subs(r,1)>0
for f,expected in zip(F[1:],resultants):
    assert sp.expand(sp.resultant(F1,f,x)-expected)==0

# At r=2 disjoint rational intervals isolate the smaller roots.
intervals = [
    (sp.Rational(354,1000),sp.Rational(355,1000)),
    (sp.Rational(355,1000),sp.Rational(356,1000)),
    (sp.Rational(1752,1000),sp.Rational(1754,1000)),
    (sp.Rational(1777,1000),sp.Rational(1779,1000)),
]
for f,(lo,hi) in zip(F,intervals):
    f2=f.subs(r,2)
    A,B,C=sp.Poly(f2,x).all_coeffs()
    assert f2.subs(x,lo)>0 and f2.subs(x,hi)<0
    assert hi < -B/(2*A)  # root isolated on decreasing branch
assert intervals[0][1]==intervals[1][0]
assert intervals[0][1]<intervals[2][0]
assert intervals[0][1]<intervals[3][0]

# Nonzero resultants prohibit crossing for r>1. Continuous smaller
# roots keep the exact order established at r=2 throughout r>1.
# The r=1 endpoint has an exact F1/F2 first-root tie.
for f,roots in zip(F,
                   [(sp.Rational(1,4),sp.Integer(4))]*2
                   +[(sp.Rational(11,9),sp.Integer(19))]*2):
    A=sp.Poly(f.subs(r,1),x).LC()
    assert sp.expand(f.subs(r,1)-A*(x-roots[0])*(x-roots[1]))==0
print("PASS: 4 discriminants, 3 noncrossing resultants,")
print("4 rational isolating intervals and isotropic endpoint.")
print("For a>=b>0 the smallest positive spin-4 factor root is F1;")
print("strict versus F2/F3/F4 when a>b, tied with F2 at a=b.")
print("This script proves only within-spin ordering; see cross_spin_order.py for spins 1..4.")
