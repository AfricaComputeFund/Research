#!/usr/bin/env python3
"""Self-contained symbolic and numerical certificate for the pure SU(2) Q_s paper.

This checks exact identities with SymPy, finite-matrix spectral signs with NumPy,
and root ordering by algebraic resultants plus rational sign isolations.
No Lean proof or physical embedding is claimed.
"""
import json
import math
import platform
from pathlib import Path

import numpy as np
import sympy as S

r,x=S.symbols('r x',real=True)
a,b,c=S.symbols('a b c',real=True)

F1=2*(r+1)*x*x-(2*r*r+13*r+2)*x+2*r*(r+1)
F2=64*x*x-136*(r+1)*x+(9*r*r+46*r+9)
F3=9*x*x-(46*r+136)*x+(9*r*r+136*r+64)
F4=9*x*x-(136*r+46)*x+(64*r*r+136*r+9)
F=(F1,F2,F3,F4)
g1=r+1-x
g2=3*r*r+8*r-(26*r+8)*x+3*x*x
g3=3*r*r+26*r+3-8*(r+1)*x
g4=8*r+3-(8*r+26)*x+3*x*x
g=(g1,g2,g3,g4)


def matrices_sym(s):
    U=S.zeros(2*s+1)
    for i,m in enumerate(range(-s,s)):
        U[i+1,i]=S.sqrt((s-m)*(s+m+1))
    return (U+U.T)/2,(U-U.T)/(2*S.I),S.diag(*range(-s,s+1))


def matrix_num(s,aa,bb,cc):
    U=np.zeros((2*s+1,2*s+1),dtype=np.complex128)
    for i,m in enumerate(range(-s,s)):
        U[i+1,i]=np.sqrt((s-m)*(s+m+1))
    j1=(U+U.conj().T)/2
    j2=(U-U.conj().T)/(2j)
    j3=np.diag(np.arange(-s,s+1).astype(float))
    return aa*(j1@j1)+bb*(j2@j2)-cc*(j3@j3)


def roots_quad_small(poly,R):
    z=S.Poly(poly.subs(r,R),x)
    aa,bb,cc=map(float,z.all_coeffs())
    return 2*cc/(-bb+math.sqrt(bb*bb-4*aa*cc))


def c4(aa,bb):
    B=2*aa*aa+13*aa*bb+2*bb*bb
    D=B*B-16*aa*bb*(aa+bb)*(aa+bb)
    return 4*aa*bb*(aa+bb)/(B+math.sqrt(D))


def c3(aa,bb):
    R=aa/bb
    return bb*(8*R+3)/(4*R+13+4*math.sqrt(R*R+5*R+10))


def exact_assert(expr):
    if S.factor(expr)!=0:
        raise AssertionError(f'Nonzero exact polynomial: {S.factor(expr)}')


def main():
    exact_checks=0
    # Reconstruct all determinants from ladder operators (rather than trusting factor transcript).
    expected=[(r+1)*(r-x)*(1-x),
              12*(r+1-4*x)*(r+4-x)*(4*r+1-x)*(r-(r+1)*x),
              108*S.prod(g),320*S.prod(F)]
    for s,expect in enumerate(expected,start=1):
        j1,j2,j3=matrices_sym(s)
        q=r*j1**2+j2**2-x*j3**2
        det=S.factor(q.det(method='berkowitz'))
        exact_assert(det-expect)
        exact_checks+=1
        print(f'PASS exact spin-{s} determinant (dimension {2*s+1})',flush=True)

    crossF=[-2268*(r-1)**6,
            -63*(r+4)**4*(4*r+1)**2,
            -63*(r+4)**2*(4*r+1)**4]
    for f,expect in zip(F[1:],crossF):
        exact_assert(S.resultant(F1,f,x)-expect)
        exact_checks+=1
    print('PASS spin-4 factor resultants',flush=True)

    crossG=[-7*r*(r+1),
            -r*r*(240*r**4+120*r**3+1703*r*r+928*r+320),
            -2*(r+1)*(15*r**4+144*r**3+578*r*r+144*r+15),
            -(320*r**4+928*r**3+1703*r*r+120*r+240)]
    for f,expect in zip(g,crossG):
        exact_assert(S.resultant(F1,f,x)-expect)
        exact_checks+=1
    print('PASS spin-4/spin-3 resultants',flush=True)

    crossG3=[-5*(r+2)**2,-1215*(r-1)**4,-15*(r+3)**3*(11*r+1)]
    for f,expect in zip(g[:3],crossG3):
        exact_assert(S.resultant(g4,f,x)-expect)
        exact_checks+=1
    print('PASS spin-3 global first-factor resultants',flush=True)

    D=[4*r**4+36*r**3+145*r**2+36*r+4,
       64*(253*r*r+394*r+253),
       64*(28*r*r+119*r+253),
       64*(253*r*r+119*r+28)]
    for f,d in zip(F,D):
        exact_assert(S.discriminant(f,x)-d)
        exact_checks+=1
    for f,d in zip((g2,g4),(64*(10*r*r+5*r+1),64*(r*r+5*r+10))):
        exact_assert(S.discriminant(f,x)-d)
        exact_checks+=1

    # The rational intervals certify all first-root comparisons at r=2.
    intervals=[(S.Rational(354,1000),S.Rational(355,1000)),
               (S.Rational(355,1000),S.Rational(356,1000)),
               (S.Rational(1752,1000),S.Rational(1754,1000)),
               (S.Rational(1777,1000),S.Rational(1779,1000))]
    for f,(lo,hi) in zip(F,intervals):
        poly=f.subs(r,2)
        A,B,C=S.Poly(poly,x).all_coeffs()
        assert poly.subs(x,lo)>0 and poly.subs(x,hi)<0 and hi<-B/(2*A)
        exact_checks+=1
    assert intervals[0][1]==intervals[1][0]

    # Spin-three g4 is first: compare r=2 and nonzero resultants r>1.
    # Exact rational isolations, not floating-point root guesses.
    # g4 first root in (468/1000,469/1000), g2 in (478/1000,479/1000).
    for polynomial,lo,hi in [(g4,S.Rational(468,1000),S.Rational(469,1000)),
                              (g2,S.Rational(478,1000),S.Rational(479,1000))]:
        f=polynomial.subs(r,2)
        A,B,C=S.Poly(f,x).all_coeffs()
        assert f.subs(x,lo)>0 and f.subs(x,hi)<0 and hi < -B/(2*A)
        exact_checks+=1
    assert S.solve(g3.subs(r,2),x)==[S.Rational(67,24)]
    assert S.solve(g1.subs(r,2),x)==[S.Integer(3)]
    exact_checks+=2
    assert S.factor(g4.subs(r,1))==(x-11)*(3*x-1)
    assert S.factor(g2.subs(r,1))==(x-11)*(3*x-1)
    assert S.factor(F1.subs(r,1))==(x-4)*(4*x-1)
    exact_assert(F2.subs(r,1)-16*(x-4)*(4*x-1))
    exact_checks+=4

    # Show spin3 crossing spin2 before its first zero, and strict spin2<spin1.
    x2=r/(r+1)
    exact_assert(g4.subs(x,x2)+3*(4*r*r+4*r-1)/(r+1)**2)
    exact_checks+=1

    # Stable c4 formula is exactly the low root of F1 (algebraic check).
    B=2*r*r+13*r+2
    discriminant=B*B-16*r*(r+1)**2
    exact_assert(S.expand(discriminant-D[0]))
    exact_checks+=1
    root=(B-S.sqrt(discriminant))/(4*(r+1))
    exact_assert(F1.subs(x,root))
    exact_checks+=1

    # All-spin inequalities at two endpoints are exact Casimir/Rayleigh evaluations.
    ss,aa,bb=S.symbols('s a b',positive=True)
    exact_assert(bb*ss*(ss+1)-(bb+bb/ss)*ss**2)
    exact_assert(ss*(aa+bb)/2-((aa+bb)/(2*ss))*ss**2)
    exact_checks+=2

    spectral=[]
    pairs=[(1,1),(2,1),(3,2),(5,1),(10,9)]
    for aa,bb in pairs:
        tc4=c4(aa,bb)
        tc3=c3(aa,bb)
        tc2=aa*bb/(aa+bb)
        tc1=bb
        assert 0<tc4<tc3<tc2<tc1
        for s,thr in [(1,tc1),(2,tc2),(3,tc3),(4,tc4)]:
            below=np.linalg.eigvalsh(matrix_num(s,aa,bb,thr*0.99))[0]
            at=np.linalg.eigvalsh(matrix_num(s,aa,bb,thr))[0]
            above=np.linalg.eigvalsh(matrix_num(s,aa,bb,thr*1.01))[0]
            assert below>1e-8 and abs(at)<1e-8 and above<-1e-8,(aa,bb,s,thr,below,at,above)
            spectral.append(dict(a=aa,b=bb,spin=s,threshold=thr,lambda_below=below,lambda_at=at,lambda_above=above))
    print('PASS 20 spectral sign checks across s=1..4',flush=True)

    allspin=0
    for aa,bb in [(1,1),(2,1),(5,1),(10,9)]:
        for s in (1,2,3,4,6,8):
            lower=bb/s
            upper=(aa+bb)/(2*s)
            lowE=float(np.linalg.eigvalsh(matrix_num(s,aa,bb,lower*(1-1e-6)))[0])
            hiE=float(np.linalg.eigvalsh(matrix_num(s,aa,bb,upper*(1+1e-5)))[0])
            assert lowE>0 and hiE<0,(aa,bb,s,lowE,hiE)
            allspin+=1
    print('PASS 24 all-spin variational-bound spectral tests',flush=True)

    data={"date":"2026-11-03","scope":"pure Q_s integer spins; not a physical fluctuation spectrum",
          "python":platform.python_version(),"sympy":S.__version__,"numpy":np.__version__,
          "exact_identity_checks":exact_checks,"crossing_checks":len(spectral),
          "all_spin_bound_checks":allspin,"spin_1_to_4_examples":spectral,
          "formal_lean":"not_run","independent_external_peer_review":"not_performed",
          "result":"PASS"}
    Path('verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print('Wrote verification.json: PASS',flush=True)

if __name__=='__main__':
    main()
