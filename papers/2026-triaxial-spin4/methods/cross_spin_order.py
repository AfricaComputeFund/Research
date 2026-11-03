#!/usr/bin/env python3
"""Exact pure-Q spin 1..4 threshold-order certificate; not physical theory."""
import sympy as S

r, x = S.symbols("r x", real=True)
a, b, c = S.symbols("a b c", real=True)

def matrix(j):
    ms = list(range(-j, j+1))
    U = S.zeros(2*j+1)
    for k, m in enumerate(ms[:-1]):
        U[k+1, k] = S.sqrt((j-m)*(j+m+1))
    X = (U+U.T)/2
    Y = (U-U.T)/(2*S.I)
    Z = S.diag(*ms)
    return a*X**2+b*Y**2-c*Z**2

F1 = 2*(r+1)*x*x-(2*r*r+13*r+2)*x+2*r*(r+1)
g1 = r+1-x
g2 = 3*r*r+8*r-(26*r+8)*x+3*x*x
g3 = 3*r*r+26*r+3-8*(r+1)*x
g4 = 8*r+3-(8*r+26)*x+3*x*x
spin3 = [g1,g2,g3,g4]
spin2 = [r+1-4*x,r+4-x,4*r+1-x,r-(r+1)*x]
spin1 = [r-x,1-x]

# Reconstruct lower-spin determinants directly from SU(2) ladder matrices.
for j,factors,k in [(1,spin1,r+1),(2,spin2,12),(3,spin3,108)]:
    det = S.factor(matrix(j).det().subs({a:r,b:1,c:x}))
    assert S.expand(det-k*S.prod(factors)) == 0, (j,det)

F2 = 64*x*x-136*(r+1)*x+9*r*r+46*r+9
F3 = 9*x*x-(46*r+136)*x+9*r*r+136*r+64
F4 = 9*x*x-(136*r+46)*x+64*r*r+136*r+9
assert S.expand(matrix(4).det().subs({a:r,b:1,c:x})-320*F1*F2*F3*F4)==0

# Every resultant is nonzero for r>=1: no spin-4 F1 root crosses
# any spin-3 root. At r=1, F1 first root 1/4; spin-3 first 1/3.
expected = [
    -7*r*(r+1),
    -r*r*(240*r**4+120*r**3+1703*r*r+928*r+320),
    -2*(r+1)*(15*r**4+144*r**3+578*r*r+144*r+15),
    -(320*r**4+928*r**3+1703*r*r+120*r+240),
]
for g,e in zip(spin3,expected):
    assert S.factor(S.resultant(F1,g,x)-e)==0
assert S.factor(F1.subs(r,1))==(x-4)*(4*x-1)
assert S.factor(g2.subs(r,1))==(x-11)*(3*x-1)
assert S.factor(g4.subs(r,1))==(x-11)*(3*x-1)
assert S.solve(g1.subs(r,1),x)==[2]
assert S.solve(g3.subs(r,1),x)==[2]
assert S.expand(S.discriminant(g2,x)-64*(10*r*r+5*r+1))==0
assert S.expand(S.discriminant(g4,x)-64*(r*r+5*r+10))==0

# At c2=r/(r+1), g4 is negative while g4(0)>0.
# Therefore the first spin-3 zero lies strictly below c2.
c2 = r/(r+1)
assert S.factor(g4.subs(x,c2)+3*(4*r*r+4*r-1)/(r+1)**2)==0
# Spin-2 threshold c2=min{(r+1)/4,r+4,4r+1,r/(r+1)}
# by (r-1)^2>=0; spin-1 threshold is 1.
print("PASS: independent exact spin-1..4 determinant reconstruction,")
print("four cross-spin resultants, isotropic anchor, discriminants,")
print("and spin-3 sign witness at c2.")
print("For a>=b>0: c4<c3<c2<c1 in the pure-Q model.")
