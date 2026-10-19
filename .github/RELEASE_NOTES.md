**New lower bounds for Erdős–Szekeres numbers in three and four dimensions**, by Deep Bhattacharjee (corresponding author) and Ushashi Bhattacharya.

Let g⁽ᵈ⁾(k) be the least number such that every set of that many points in general position in ℝᵈ contains k points in convex position. Scheucher (2023) proved g⁽³⁾(7) = 13, g⁽⁴⁾(8) ≤ 13 and g⁽⁵⁾(9) ≤ 13; the best lower bounds for g⁽⁴⁾(8) and g⁽³⁾(8) were 11 and 14.

**Theorem A.** Eleven integer points in general position in ℝ⁴ (coordinates 0 to 40) with no eight in convex position. Hence g⁽⁴⁾(8) ≥ 12 and h⁽⁴⁾(8) ≥ 12, so both are 12 or 13.

**Theorem B.** Fifteen integer points in general position in ℝ³ (coordinates 0 to 180) with no eight in convex position. Hence g⁽³⁾(8) ≥ 16, h⁽³⁾(8) ≥ 16 and g⁽³⁾(9) ≥ 17.

| Set | Relations | 8-subsets | Convex 8-subsets |
|---|---|---|---|
| 11 points in ℝ⁴ | 28 | 165 | 0 |
| 15 points in ℝ³ | 161 | 6435 | 0 |

Each relation Σⱼ mⱼ x_{sⱼ} = M x_p, with positive integer weights, places x_p strictly inside a simplex of other points; every eight points contain one. The exact values of g⁽⁴⁾(8) and g⁽³⁾(8) remain open.

Checks: `checks/python/check_all.py` (980 exact checks), C (brute force over all subsets), C++ (certificates), Julia (BigInt), and a Lean 4 file in which the kernel proves general position, every relation and the covering condition (no Mathlib; six theorems use no axiom, the two covering theorems only `propext`).

Files:
- `ErdosSzekeresHighDim_final2.pdf`: the paper
- `ErdosSzekeresHighDim_final2-tex.zip`: LaTeX source (`main.tex`, tables, figures as PNG)
- `ErdosSzekeresHighDim_final2-arxiv.tar.gz`: the same source, ready for arXiv

The repository also holds `preprintErdosSzekeresHighDim/ErdosSzekeresHighDim.pdf` and `MasterErdosSzekeresHighDim_masterzipfinal2.zip`, which holds every file of the repository. Run `scripts/run_all.sh` to repeat the checks and `scripts/build_paper.sh` to rebuild the files above.
