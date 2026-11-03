# Exact Triaxial Spin-4 Stability Thresholds for an Anisotropic SU(2) Quadratic Operator

**Africa Compute Fund — Research Programme**  
**Technical preprint, 3 November 2026**

This is a self-contained mathematical preprint and reproduction package. It studies the Hermitian finite-spin family

\[
Q_s(c)=aJ_1^2+bJ_2^2-cJ_3^2,\qquad a\geq b>0,\quad c\geq0,
\]

and the first value of \(c\) at which strict positive definiteness is lost. The paper gives the global exact spin-4 threshold, closed formulas for spins 1–4, a strict spin-threshold ordering for integer spins 1–4, and a variational envelope valid at every integer spin.

**Scope:** These are mathematical claims about this specified matrix family. They do not identify a physical fluctuation spectrum or prove the stability or instability of any physical background.

## Contents

| File | Purpose |
| --- | --- |
| `paper.pdf` | Typeset, self-contained technical paper (8 pages). |
| `paper.tex` | Editable LaTeX manuscript. |
| `verify.py` | Exact symbolic determinant/resultant checks and independent numerical eigenvalue controls. |
| `verification.json` | Actual recorded verification results, numerical values and runtime versions. |
| `make_figure.py` | Rebuilds the threshold comparison plot. |
| `thresholds.pdf` / `thresholds.png` | Reproducible threshold figure. |
| `requirements.txt` | Python package dependencies. |
| `CITATION.cff` | GitHub citation metadata. |
| `SHA256SUMS` | File-integrity hashes (excluding itself). |

## Reproduce

With Python 3 installed:

```bash
python -m pip install -r requirements.txt
python verify.py
python make_figure.py
```

The verification command prints the exact determinants, resultants and spectral tests as they pass, and writes `verification.json`. To rebuild the PDF, install a LaTeX distribution with `pdflatex`, `amsmath`, `amsthm`, `mathtools`, `microtype`, `geometry`, `booktabs`, `graphicx`, `enumitem`, `hyperref`, and `cleveref` then run:

```bash
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
```

The original release-generation environment used Python 3.11+, SymPy 1.14.0, NumPy 2.3.5 and Matplotlib 3.10.8. See `verification.json` for the actual verification runtime.

## Verification and review status

- **Exact algebra:** finite-spin determinant identities for spins 1–4, spin-4 cross-factor resultants, spin-3 cross-factor resultants, the spin-4/spin-3 resultant, and supporting discriminant and endpoint identities passed under SymPy.
- **Numerical controls:** 20 threshold-crossing checks for spins 1–4 across five anisotropy pairs, plus 24 all-spin spectral-bound checks passed under NumPy.
- **Proof status:** detailed, self-contained mathematical arguments and exact computer algebra; **not** a Lean-formalized proof and **not** externally peer-reviewed.
- **Prior art:** the article cites relevant published literature on asymmetric-top and collective-spin quadratic Hamiltonians, but does not claim a completed novelty or priority audit.
- **Physical scope:** the operator is not asserted to be a physically derived gauge-fixed fluctuation operator.

## Publication guidance

This folder is the intended **public subset**. It intentionally excludes internal project identifiers, private source-commit references, broader research repositories, and unpublished physical-model work. Publishing this folder does not require making the underlying private corpus public.

After publishing it in GitHub, use the versioned URL to `paper.pdf` in the announcement or link to the repository itself. Record the GitHub release tag/commit so that future revisions remain traceable. Choose an appropriate license **before** offering reuse rights; no license has been assigned here because the owner has not selected one.

Authorship is currently recorded institutionally as **Africa Compute Fund**. The owner should review the title-page authorship and citation metadata before public release.
