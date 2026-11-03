# Exact Triaxial Spin-4 Stability Thresholds for an Anisotropic SU(2) Quadratic Operator

**Africa Compute Fund · Research Programme**  
**Technical preprint · 3 November 2026**

## Abstract

For the integer-spin irreducible representations of $\mathfrak{su}(2)$, we study the positive-definiteness threshold of the Hermitian quadratic operator $Q_s(c)=aJ_1^2+bJ_2^2-cJ_3^2$, with $a\geq b>0$ and $c\geq0$. We factor the spin-four determinant exactly, identify its first positive zero over the anisotropy domain, and obtain a closed radical expression for the threshold $c_4$. A further factor analysis gives a closed formula for $c_3$. Exact resultant and sign arguments establish $c_4<c_3<c_2<c_1$ for every $a\geq b>0$, including isotropy. At arbitrary integer spin, a variational estimate proves $b/s\leq c_s\leq(a+b)/(2s)$, so fixed $c>0$ cannot preserve positive definiteness through an unbounded pure-quadratic tower. Determinant identities and comparison certificates are independently reproducible using symbolic and numerical methods.

---

**Keywords:** angular momentum; anisotropic quadratic operators; finite-dimensional spectral theory; exact algebraic thresholds; symbolic verification; matrix models.

# Introduction

Anisotropic quadratic forms in angular-momentum generators appear in the spectral theory of quantum asymmetric tops, collective-spin systems and related matrix models. Their algebraic structure is well established; see, for example, the asymmetric-top treatment of Breev and Gitman , the collective-spin perspective of Opatrný *et al.* , and the exact-spectrum analysis of Ribeiro, Vidal and Mosseri . A negative coefficient in such a quadratic form changes a precise mathematical question: when does its spectrum first reach zero?

We address that question for a fixed family of finite-dimensional Hermitian matrices and a specified coefficient domain. The primary object is the first loss of *strict positive definiteness*, not a dynamical growth exponent or an experimentally demonstrated instability. The distinction is essential: a physical system can have additional Hessian terms, constraints, gauge or ghost sectors, and spin-dependent coefficients absent from the operator studied here.

The concrete results are: (i) an explicit spin-four determinant and its globally minimal factor root; (ii) a closed spin-three formula and a strict ordering of all four integer-spin thresholds; and (iii) an all-spin variational envelope. Proofs use exact polynomial identities, nonvanishing resultants, rational root isolation and the Rayleigh–Ritz principle. We make no claim of priority for the general operator class, the standard variational technique, or an exhaustively established novelty assessment.

# Operator, notation and meaning of the threshold

For integer $s\geq1$, let $J_1,J_2,J_3$ be the Hermitian generators of the spin-$s$ irreducible representation, acting on $\mathcal H_s\cong\mathbb C^{2s+1}$. Our conventions are $$
[J_i,J_j]=i\varepsilon_{ijk}J_k,\qquad
J_1^2+J_2^2+J_3^2=s(s+1)I.$$ In the $J_3$ eigenbasis, $m=-s,-s+1,\ldots,s$, $$\begin{aligned}
J_3|s,m\rangle&=m|s,m\rangle,\\
J_\pm|s,m\rangle&=\sqrt{(s\mp m)(s\pm m+1)}
\,|s,m\pm1\rangle,
\end{aligned}$$ where $J_\pm=J_1\pm iJ_2$.

Throughout the paper, $$
Q_s(c):=aJ_1^2+bJ_2^2-cJ_3^2,
\qquad a\geq b>0,\quad c\geq0.$$ We write $A\succ0$ for strict positive definiteness and use the following unambiguous definition: $$
c_s(a,b):=\inf\{c\geq0:Q_s(c)\not\succ0\}.$$ The spin-zero block is excluded because $Q_0(c)=0$ identically. By homogeneity, $c_s(a,b)=b\,x_s(r)$ for $$
r:=a/b\geq1,\qquad x:=c/b.$$ A zero of $\det Q_s(c)$ is a threshold candidate; identifying the *first* zero and proving positivity on its left are separate steps.



**Lemma 1** (Initial positivity and monotonicity). *For every integer $s\geq1$, $Q_s(0)\succ0$ and $\lambda_{\min}(Q_s(c))$ is nonincreasing in $c$. Consequently, $c_s$ is the smallest nonnegative determinant zero and $Q_s(c)\succ0$ exactly when $0\leq c<c_s$.*





*Proof.* For any $v$, $\langle v,Q_s(0)v\rangle\geq
b(\|J_1v\|^2+\|J_2v\|^2)$. If this vanishes, then $J_1v=J_2v=0$ and hence $J_3v=(1/i)[J_1,J_2]v=0$. The Casimir identity <a href="#eq:su2" data-reference-type="eqref" data-reference="eq:su2">[eq:su2]</a> then forces $s(s+1)v=0$, so $v=0$. For $c'\geq c$, $Q_s(c')-Q_s(c)=-(c'-c)J_3^2\preceq0$. The Rayleigh–Ritz characterization of the least eigenvalue proves monotonicity. Since the highest-weight expectation tends to $-\infty$ as $c\to\infty$, a threshold exists and is positive; continuity and monotonicity complete the argument. ◻



Expanding <a href="#eq:Q" data-reference-type="eqref" data-reference="eq:Q">[eq:Q]</a> in ladder operators gives a useful sparse form valid for every integer $s$: $$\begin{aligned}
\langle s,m|Q_s(c)|s,m\rangle
 &=\frac{a+b}{2}\,[s(s+1)-m^2]-cm^2,\\
\langle s,m+2|Q_s(c)|s,m\rangle
 &=\frac{a-b}{4}
\sqrt{(s-m)(s+m+1)(s-m-1)(s+m+2)}.
\end{aligned}$$ All other off-diagonal entries vanish. Thus the even-$m$ and odd-$m$ subspaces are invariant, and the determinants factor into small polynomial blocks.

# Exact spin-four threshold

Let $s=4$, so that $\dim\mathcal H_4=9$. An exact determinant evaluation using <a href="#eq:diag" data-reference-type="eqref" data-reference="eq:diag">[eq:diag]</a>–<a href="#eq:offdiag" data-reference-type="eqref" data-reference="eq:offdiag">[eq:offdiag]</a> gives $$
\det Q_4(c)=320F_1F_2F_3F_4,$$ where $$\begin{aligned}
F_1&=2(a+b)c^2-(2a^2+13ab+2b^2)c+2ab(a+b),\\
F_2&=64c^2-136(a+b)c+9a^2+46ab+9b^2,\\
F_3&=9c^2-(46a+136b)c+9a^2+136ab+64b^2,\\
F_4&=9c^2-(136a+46b)c+64a^2+136ab+9b^2.
\end{aligned}$$ These are exact polynomial identities, not numerical interpolations. The $c$-degree is eight rather than nine because the $m=0$ basis vector has zero $J_3^2$ expectation.



**Theorem 2** (Global spin-four first crossing). *For every $a\geq b>0$, the first positive determinant zero of $Q_4$ is the smaller root of $F_1$. In particular, with $$
B=2a^2+13ab+2b^2,
\qquad \Delta=B^2-16ab(a+b)^2,$$ its positive-definiteness threshold is $$
\boxed{
 c_4=\frac{B-\sqrt\Delta}{4(a+b)}
 =\frac{4ab(a+b)}{B+\sqrt\Delta}.}$$ For $a>b$ this root is strictly smaller than the positive roots of $F_2,F_3,F_4$. When $a=b$, the least roots of $F_1$ and $F_2$ coincide.*





*Proof.* Use <a href="#eq:scaling" data-reference-type="eqref" data-reference="eq:scaling">[eq:scaling]</a> and denote by $F_i(r,x)$ the polynomials <a href="#eq:F1" data-reference-type="eqref" data-reference="eq:F1">[eq:F1]</a>–<a href="#eq:F4" data-reference-type="eqref" data-reference="eq:F4">[eq:F4]</a> with $a=r,b=1,c=x$. Their quadratic discriminants are respectively $$\begin{aligned}
D_1&=4r^4+36r^3+145r^2+36r+4,\\
D_2&=64(253r^2+394r+253),\\
D_3&=64(28r^2+119r+253),\\
D_4&=64(253r^2+119r+28).
\end{aligned}$$ For $r\geq1$ these are strictly positive; each quadratic has positive leading and constant terms and a negative linear coefficient. Both roots of each factor are therefore distinct and positive, with smaller branches continuous in $r$.

The exact resultants, eliminating $x$, are $$\begin{aligned}
\operatorname{Res}_x(F_1,F_2)&=-2268(r-1)^6,\nonumber\\
\operatorname{Res}_x(F_1,F_3)&=-63(r+4)^4(4r+1)^2,\\
\operatorname{Res}_x(F_1,F_4)&=-63(r+4)^2(4r+1)^4.\nonumber
\end{aligned}$$ They never vanish for $r>1$; no root of $F_1$ can meet a root of another factor there. At the exact rational point $r=2$, positive-to-negative sign changes below each quadratic vertex isolate the smaller positive roots in $$
\begin{aligned}
x_{1,-}&\in(354/1000,355/1000),&
x_{2,-}&\in(355/1000,356/1000),\\
x_{3,-}&\in(1752/1000,1754/1000),&
x_{4,-}&\in(1777/1000,1779/1000).
\end{aligned}$$ The $F_1$ smaller root is strictly left of the others. Continuity and <a href="#eq:res4" data-reference-type="eqref" data-reference="eq:res4">[eq:res4]</a> preserve that order for all $r>1$. At $r=1$, $F_1$ and $F_2$ have roots $\{1/4,4\}$, while $F_3$ and $F_4$ have roots $\{11/9,19\}$. The earliest zero is therefore always the smaller $F_1$ root. Solving <a href="#eq:F1" data-reference-type="eqref" data-reference="eq:F1">[eq:F1]</a> and rationalizing gives <a href="#eq:c4" data-reference-type="eqref" data-reference="eq:c4">[eq:c4]</a>; <a href="#lem:monotone" data-reference-type="ref" data-reference="lem:monotone">1</a> identifies it as the threshold. ◻



The rationalized expression in <a href="#eq:c4" data-reference-type="eqref" data-reference="eq:c4">[eq:c4]</a> is numerically preferable when $a/b$ is large because it avoids cancellation in $B-\sqrt\Delta$.

# Exact thresholds at spins one, two and three

The lower-spin determinants may also be built directly from the ladder matrices. With $r=a/b$ and $x=c/b$, $$\begin{aligned}
\frac{\det Q_1}{b^3}
 &=(r+1)(r-x)(1-x),\\
\frac{\det Q_2}{b^5}
 &=12(r+1-4x)(r+4-x)(4r+1-x)
       [r-(r+1)x],\\
\frac{\det Q_3}{b^7}
 &=108g_1g_2g_3g_4,
\end{aligned}$$ with $$\begin{aligned}
g_1&=r+1-x,\\
g_2&=3r^2+8r-(26r+8)x+3x^2,\\
g_3&=3r^2+26r+3-8(r+1)x,\\
g_4&=8r+3-(8r+26)x+3x^2.
\end{aligned}$$ It follows immediately that $$
\boxed{c_1=b,\qquad c_2=\frac{ab}{a+b}.}$$ Indeed $r/(r+1)\leq(r+1)/4$ because $(r-1)^2\geq0$, and the other spin-two roots are larger.



**Proposition 3** (Closed spin-three threshold). *For every $a\geq b>0$, the first positive spin-three zero is the smaller root of $g_4$, namely $$
\boxed{c_3=b\,\frac{4r+13-4\sqrt{r^2+5r+10}}{3}
 =b\,\frac{8r+3}{4r+13+4\sqrt{r^2+5r+10}}},
\qquad r=a/b.$$ The competing smaller root of $g_2$ is equal to this value at $a=b$, but strictly larger when $a>b$.*





*Proof.* The two roots of each of $g_2,g_4$ are positive and distinct for $r\geq1$: their respective discriminants are $64(10r^2+5r+1)$ and $64(r^2+5r+10)$. The roots of the two linear factors $g_1,g_3$ are positive. The exact resultants are $$\begin{aligned}
\operatorname{Res}_x(g_4,g_1)&=-5(r+2)^2,\nonumber\\
\operatorname{Res}_x(g_4,g_2)&=-1215(r-1)^4,\\
\operatorname{Res}_x(g_4,g_3)&=-15(r+3)^3(11r+1).\nonumber
\end{aligned}$$ They are nonzero on $r>1$. At $r=2$, exact rational sign changes below the quadratic vertices isolate the least roots: $$x_{g_4,-}\in(468/1000,469/1000),\qquad
x_{g_2,-}\in(478/1000,479/1000).$$ The linear factors $g_3$ and $g_1$ vanish at $67/24$ and $3$. Thus $g_4$ has the earliest zero. As none of these roots can cross those of $g_4$ for $r>1$, its smaller root is globally first. At $r=1$, the least roots of $g_2$ and $g_4$ coincide at $1/3$; the other factors first vanish at $2$. The quadratic formula, with its rationalized variant, gives <a href="#eq:c3" data-reference-type="eqref" data-reference="eq:c3">[eq:c3]</a>. ◻



# Global threshold ordering for the finite tower



**Theorem 4** (Strict order through spin four). *For every $a\geq b>0$, the positive-definiteness thresholds of the four integer-spin blocks satisfy $$
\boxed{0<c_4<c_3<c_2<c_1.}$$ Consequently, the direct sum $\bigoplus_{s=1}^4Q_s(c)$ is strictly positive definite if and only if $0\leq c<c_4$.*





*Proof.* By <a href="#thm:spin4,prop:spin3" data-reference-type="ref" data-reference="thm:spin4,prop:spin3">[thm:spin4,prop:spin3]</a>, the relevant first roots are the smaller roots of $F_1$ and $g_4$. The exact resultant $$
\operatorname{Res}_x(F_1,g_4)
=-(320r^4+928r^3+1703r^2+120r+240)$$ is strictly negative for every $r\geq1$. Both smaller roots are continuous in $r$, and at $r=1$ they are $1/4$ and $1/3$, respectively. They cannot cross, so $c_4<c_3$ for the entire domain.

Next, with $x_2=r/(r+1)$, $$
g_4(0)=8r+3>0,
\qquad
g_4(x_2)=-\frac{3(4r^2+4r-1)}{(r+1)^2}<0.$$ By continuity, the first $g_4$ zero precedes $x_2$, so $c_3<c_2$. Finally $r/(r+1)<1$, giving $c_2<c_1$. Initial positivity and <a href="#lem:monotone" data-reference-type="ref" data-reference="lem:monotone">1</a> establish the direct-sum criterion. ◻



<figure id="fig:thresholds">
<embed src="thresholds.pdf" />
<figcaption>Dimensionless thresholds <span class="math inline"><em>x</em><sub><em>s</em></sub> = <em>c</em><sub><em>s</em></sub>/<em>b</em></span> versus <span class="math inline"><em>r</em> = <em>a</em>/<em>b</em></span> for the exact integer-spin sectors <span class="math inline"><em>s</em> = 1, 2, 3, 4</span>. The plotted curves are evaluations of the analytical formulas, not evidence for a theorem beyond these four sectors.</figcaption>
</figure>

<div id="tab:examples">

| $a$ | $b$ |    $c_1$ |    $c_2$ |    $c_3$ |    $c_4$ |
|:---:|:---:|---------:|---------:|---------:|---------:|
|  1  |  1  | 1.000000 | 0.500000 | 0.333333 | 0.250000 |
|  2  |  1  | 1.000000 | 0.666667 | 0.468027 | 0.354249 |
|  3  |  2  | 2.000000 | 1.200000 | 0.815741 | 0.613062 |
|  5  |  1  | 1.000000 | 0.833333 | 0.672044 | 0.543069 |
| 10  |  9  | 9.000000 | 4.736842 | 3.162430 | 2.371925 |

Representative thresholds from the closed formulas. Rows test the isotropic boundary and several anisotropic points.



# An all-spin variational envelope

The proof of <a href="#thm:order" data-reference-type="ref" data-reference="thm:order">4</a> stops at spin four. A different, elementary argument yields nontrivial information for every integer spin without assuming that adjacent-spin thresholds remain ordered.



**Theorem 5** (All-spin threshold bounds). *For any integer $s\geq1$ and fixed $a\geq b>0$, $$
\boxed{\frac b s\leq c_s(a,b)\leq\frac{a+b}{2s}.}$$ The bounds coincide at $a=b$, where $c_s=b/s$. In particular, $c_s(a,b)\to0$ as $s\to\infty$. For any fixed $c>0$, the collection of all integer-spin pure quadratic blocks $\{Q_s(c):s\geq1\}$ cannot remain positive definite.*





*Proof.* Since $a\geq b$, the Casimir identity yields $$\begin{aligned}
Q_s(c)&\succeq b(J_1^2+J_2^2)-cJ_3^2\\
&=bs(s+1)I-(b+c)J_3^2\\
&\succeq [bs(s+1)-(b+c)s^2]I,
\end{aligned}$$ using $J_3^2\preceq s^2I$. The last scalar is strictly positive when $c<b/s$, proving the lower bound. On the normalized highest-weight vector $|s,s\rangle$, $\langle J_1^2\rangle=\langle J_2^2\rangle=s/2$ and $\langle J_3^2\rangle=s^2$. Hence $$\langle s,s|Q_s(c)|s,s\rangle
=\frac{s(a+b)}2-cs^2.$$ This expectation is negative for $c>(a+b)/(2s)$; the Rayleigh–Ritz principle gives the upper bound. At $a=b$, $Q_s(c)=bs(s+1)I-(b+c)J_3^2$ is diagonal and its first zero occurs at $|m|=s$ and $c=b/s$. The limit and nonuniform-positivity conclusion follow from <a href="#eq:envelope" data-reference-type="eqref" data-reference="eq:envelope">[eq:envelope]</a>. ◻



<div class="remark">

**Remark 6**. *does *not* establish $c_{s+1}<c_s$ for all $s$ and says nothing about a physical unbounded tower with extra operators, gauge constraints, finite cutoffs, or spin-dependent coefficients. It is a conditional statement about exactly <a href="#eq:Q" data-reference-type="eqref" data-reference="eq:Q">[eq:Q]</a>.*



# Computer-assisted verification and reproducibility

A self-contained Python script, `verify.py`, reproduces the finite-spin claims without importing the original internal calculation. It constructs each ladder matrix, builds $Q_s$ independently and checks the following:

1.  Exact determinant identities <a href="#eq:det1" data-reference-type="eqref" data-reference="eq:det1">[eq:det1]</a>–<a href="#eq:det4" data-reference-type="eqref" data-reference="eq:det4">[eq:det4]</a> for spins $1$ through $4$ in SymPy.

2.  All resultants used in <a href="#eq:res4" data-reference-type="eqref" data-reference="eq:res4">[eq:res4]</a>, <a href="#eq:res3" data-reference-type="eqref" data-reference="eq:res3">[eq:res3]</a> and <a href="#eq:res43" data-reference-type="eqref" data-reference="eq:res43">[eq:res43]</a>, discriminants, rational root-isolating intervals and endpoint identities.

3.  Hermitian smallest-eigenvalue signs below, at and above each of the four closed-form thresholds for five anisotropy pairs.

4.  The variational endpoints for 24 additional finite-spin cases at $s\in\{1,2,3,4,6,8\}$.

The verification record distributed alongside this paper gives runtime versions, the checks performed and their numerical outputs. The accompanying `make_figure.py` regenerates <a href="#fig:thresholds" data-reference-type="ref" data-reference="fig:thresholds">1</a>. With Python, SymPy, NumPy and Matplotlib installed, the reproduction commands are

    python verify.py
    python make_figure.py
    pdflatex -interaction=nonstopmode -halt-on-error paper.tex
    pdflatex -interaction=nonstopmode -halt-on-error paper.tex

All exact assertions are symbolic polynomial equalities; floating-point eigenvalues are used as adversarial numerical controls, not substitutes for the algebraic arguments. Numerical determinant comparisons directly at an exact zero are deliberately avoided as a relative-error metric because subtractive cancellation makes that metric ill-conditioned. The threshold is instead tested against the least eigenvalue.

The algebraic arguments in the text are human-readable proofs, backed by computer algebra. They are *not* a Lean-checked formal proof: no such formalization is claimed. Likewise, this release is not an externally peer-reviewed article. Independent expert scrutiny, additional prior-art checking and formalization are appropriate follow-on work rather than fictitious prerequisites already completed.

# Interpretation, relationship to prior work and limits

The operator in <a href="#eq:Q" data-reference-type="eqref" data-reference="eq:Q">[eq:Q]</a> belongs to the established family of quadratic angular-momentum Hamiltonians. The quantum asymmetric top literature treats the more conventional positive principal-axis coefficients ; collective-spin and generalized Lipkin–Meshkov–Glick literature demonstrates the breadth of related quadratic forms . These sources supply mathematical context; they do not by themselves establish priority or novelty of the particular first-zero formulas and comparison arguments derived here.

This paper’s contribution is specifically the explicit, fully triaxial pure-quadratic *positive-definiteness* calculation through spin four and its precise domain-wide comparison, supplemented by an all-spin bound. The phrase “stability threshold” is used in this limited spectral sense. A zero mode of $Q_s$ must not be conflated with a physical tachyon, instability of a Lorentzian vacuum, or a mode surviving physical constraints.

A physical application would require a derivation from the relevant parent action, the correct gauge-fixed Hessian and constraints, and a demonstrated identification of its physical spin sectors with <a href="#eq:Q" data-reference-type="eqref" data-reference="eq:Q">[eq:Q]</a>. In the absence of those additional steps, the finite-dimensional results cannot be promoted to claims about a particular physical background. No such inference is made here.

# Conclusion

For $a\geq b>0$ we have identified the exact first spin-four zero of the Hermitian operator $aJ_1^2+bJ_2^2-cJ_3^2$ and proved that the strict positive-definiteness thresholds obey $c_4<c_3<c_2<c_1$. The spin-three and spin-four formulas are explicit, algebraically certified, and independently tested against reconstructed finite matrices. A separate all-spin Rayleigh envelope shows that an unbounded family with fixed coefficients cannot remain positive definite for fixed $c>0$. These are finite-dimensional spectral results with transparent assumptions and reproducible methods; any proposed physical application must be established separately.

# Code and data availability

The source LaTeX, verification script, figure generator and verification JSON are distributed together with this preprint. They permit independent reconstruction of the finite matrices, algebraic identities and illustrative numerical examples. A versioned public repository URL may be added upon release. No repository URL or DOI is asserted in advance.

# Review status

This technical preprint provides a mathematically scoped result and an independently executable local audit. It has not been subject to external peer review or Lean formalization.

<div class="thebibliography">

9 A. Breev and D. Gitman, “New Exact Solutions Describing Quantum Asymmetric Top,” *Symmetry* **15** (2023), 503. [doi:10.3390/sym15020503](https://doi.org/10.3390/sym15020503). T. Opatrný, L. Richterek and M. Opatrný, “Analogies of the classical Euler top with a rotor to spin squeezing and quantum phase transitions in a generalized Lipkin–Meshkov–Glick model,” *Scientific Reports* **8** (2018), 1984. [doi:10.1038/s41598-018-20486-y](https://doi.org/10.1038/s41598-018-20486-y). P. Ribeiro, J. Vidal and R. Mosseri, “Exact spectrum of the Lipkin–Meshkov–Glick model in the thermodynamic limit and finite-size corrections,” *Physical Review E* **78** (2008), 021106. [doi:10.1103/PhysRevE.78.021106](https://doi.org/10.1103/PhysRevE.78.021106).

