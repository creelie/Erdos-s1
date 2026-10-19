# AI provenance: New lower bounds for Erdős–Szekeres numbers in three and four dimensions

Authors: Deep Bhattacharjee (corresponding author, ORCID 0000-0003-0466-750X)
and Ushashi Bhattacharya (ORCID 0009-0002-2254-3914).

- Paper: https://github.com/creelie/Erdos-s1/blob/main/preprintErdosSzekeresHighDim/ErdosSzekeresHighDim.pdf
- Repository: https://github.com/creelie/Erdos-s1
- Archive: release v1.0.0 on Zenodo

## Result

Let g⁽ᵈ⁾(k) be the least number such that every set of that many points in
general position in ℝᵈ contains k points in convex position. Scheucher (2023)
proved g⁽³⁾(7) = 13, g⁽⁴⁾(8) ≤ 13 and g⁽⁵⁾(9) ≤ 13 by satisfiability solving;
his twelve-element witnesses in rank five and six are not realizable by points.
The lower bounds g⁽⁴⁾(8) ≥ 11 and g⁽³⁾(8) ≥ 14 were the best known. The paper
gives eleven integer points in ℝ⁴ and fifteen in ℝ³, in general position,
with no eight in convex position, so g⁽⁴⁾(8) ∈ {12, 13} and g⁽³⁾(8) ≥ 16.

## Who did what

Deep Bhattacharjee chose the problem, directed the work, set the format and
the target journal, and ran the computations. Claude (Anthropic), used
through Claude Code under Deep Bhattacharjee's direction, produced the
following:

- the search programs (a graded score for the k-subsets in convex position,
  simulated annealing, repair moves) and the runs that found both sets;
- the polishing and rounding to small integer coordinates;
- the certificate generator and its minimisation;
- the proofs of the certificate criterion and of the layer structure;
- the checks in Python, C, C++, Julia and Lean;
- the figures;
- the LaTeX text of the paper.

The BLLM engine Anaxia (release v3.0.5, doi:10.5281/zenodo.23469389;
github.com/creelie/eBLLM) supplied the status check of the problem from its
indexes of DeepMind's formal-conjectures collection and AlphaEvolve problem
repository, of the OpenRefute archive and of OpenAI's collection of
mathematical manuscripts; these hold only the planar question. The paper says
so in a footnote of Section V, in the data availability section and in its
bibliography. The point sets and the certificates were not produced by the
engine; the search methods were added to it afterwards as construction rules.
No other AI system took part.

## How the result is checked

- **Python** (`checks/python/check_all.py`, 980 checks): every computational
  claim of the paper in exact integer arithmetic, including a direct
  convex-position test of every 8-subset that does not use the certificate.
- **C** (`checks/c/check.c`): brute force over all 8-subsets, without the
  certificate.
- **C++** (`checks/cpp/check.cpp`): the certificates, including the Cramer
  numerators of every relation and the covering condition over bitmasks.
- **Julia** (`checks/julia/check.jl`): relations, covering and brute force in
  BigInt.
- **Lean 4** (`checks/lean/ErdosSzekeresCertificates.lean`, core Lean
  v4.24.0, no Mathlib): the kernel proves general position, every relation
  and the covering condition; six theorems use no axiom and the two covering
  theorems only `propext`.

## What is not claimed

The paper does not determine g⁽⁴⁾(8) or g⁽³⁾(8). A twelve-point set in ℝ⁴
without an 8-gon would give g⁽⁴⁾(8) = 13; the searches for it ended with sets
that still have four 8-gons. The authors checked on 19 October 2026, in the
literature citing Scheucher (2023), on arXiv and in the indexes named above,
that no earlier work gives these lower bounds.

## Text and references

The bibliography contains, besides the works on the Erdős–Szekeres problem
that the paper builds on, three self-citations of the authors. Every entry
was checked against its DOI record or publisher page, and none comes from
TechRxiv, Qeios, Authorea or EasyChair.
