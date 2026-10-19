#!/usr/bin/env python3
"""Write the coordinate and certificate tables of the paper from the witness and
certificate files in highdim/.  Usage: python3 paper/make_tables.py"""
import json, os
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
HD = os.path.join(HERE, "..", "highdim")
OUT = os.path.join(HERE, "tables")


def det(rows):
    m = [[Fraction(x) for x in r] for r in rows]
    n, s = len(m), Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            m[c], m[p] = m[p], m[c]; s = -s
        s *= m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for j in range(c, n):
                m[r][j] -= f * m[c][j]
    return int(s)


def check(P, rels, d, k):
    """Recheck every relation and the covering condition; return (max M, list of g)."""
    n = len(P)
    gs, Ms = [], []
    for r in rels:
        p, S, m = r["p"], r["simplex"], r["weights"]
        M = sum(m)
        assert all(x > 0 for x in m)
        for j in range(d):
            assert sum(mi * P[v][j] for mi, v in zip(m, S)) == M * P[p][j]
        D = det([[1] + P[v] for v in S])
        assert D != 0 and D % M == 0
        Ms.append(M); gs.append(abs(D) // M)
    sets = [frozenset([r["p"]] + r["simplex"]) for r in rels]
    for T in combinations(range(n), k):
        T = frozenset(T)
        assert any(s <= T for s in sets), T
    return max(Ms), gs


def coords(P, name, caption, label):
    d = len(P[0])
    axes = "xyzw"[:d] if d <= 4 else [f"x_{j+1}" for j in range(d)]
    head = " & ".join(f"${a}$" for a in axes)
    lines = [r"\begin{table}[tb]", r"\caption{" + caption + r"}\label{" + label + "}",
             r"\begin{ruledtabular}", r"\begin{tabular}{r" + "r" * d + "}",
             r"point & " + head + r" \\ \hline"]
    for i, x in enumerate(P):
        lines.append(f"$x_{{{i+1}}}$ & " + " & ".join(str(c) for c in x) + r" \\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    open(os.path.join(OUT, name), "w").write("\n".join(lines) + "\n")


def idx(S):
    return ",".join(str(v + 1) for v in S)


def cert_r4(P, rels, gs, k=8, thm="thm:A", label="tab:certA", fname="cert_r4.tex", count="twenty-eight", note="", size=""):
    n = len(P)
    t = n - k
    lines = [r"\begin{table*}[tp]"] + ([size] if size else []) + [
             r"\caption{The " + count + r" relations behind Theorem~\ref{" + thm + "}" + note + r". In each row the point $x_p$ satisfies "
             r"$\sum_j m_j x_{s_j}=M x_p$ with $\sigma=\{s_1,\dots,s_5\}$ listed in increasing order, the weights $m_j$ "
             r"are positive integers with sum $M$, and $|\det\sigma|=gM$ for the $5\times5$ determinant of the rows $(1,x_s)$, "
             r"$s\in\sigma$. The last column lists the " + ["", "", "", "", "", "five", "six", "seven"][n - 6] +
             r" points not used by the relation; every " + {3: "triple", 4: "four-element subset"}[t] + r" of "
             r"$\{1,\dots," + str(n) + r"\}$ lies in one of these sets.}\label{" + label + "}",
             r"\begin{ruledtabular}", r"\begin{tabular}{rlrrrrrrrl}",
             r"$p$ & $\sigma$ & $m_1$ & $m_2$ & $m_3$ & $m_4$ & $m_5$ & $M$ & $g$ & unused points \\ \hline"]
    for r, g in zip(rels, gs):
        rest = sorted(set(range(n)) - {r["p"]} - set(r["simplex"]))
        lines.append(f"{r['p']+1} & {idx(r['simplex'])} & " + " & ".join(str(m) for m in r["weights"]) +
                     f" & {sum(r['weights'])} & {g} & {idx(rest)} \\\\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}"]
    open(os.path.join(OUT, fname), "w").write("\n".join(lines) + "\n")


def cert_r3(P, rels, gs, per=41, thm="thm:B", label="tab:certB", fname="cert_r3.tex", note=""):
    rows = [f"{r['p']+1} & {idx(r['simplex'])} & " + " & ".join(str(m) for m in r["weights"]) +
            f" & {sum(r['weights'])} & {g}" for r, g in zip(rels, gs)]
    head = r"$p$ & $\sigma$ & $m_1$ & $m_2$ & $m_3$ & $m_4$ & $M$ & $g$"
    chunks = [rows[i:i + 2 * per] for i in range(0, len(rows), 2 * per)]
    out = []
    for c, ch in enumerate(chunks):
        half = min(per, (len(ch) + 1) // 2)
        a, b = ch[:half], ch[half:]
        b += [" & " * 7] * (len(a) - len(b))
        first = (r"The " + str(len(rels)) + r" relations behind Theorem~\ref{" + thm + "}" + note if c == 0
                 else r"The relations behind Theorem~\ref{" + thm + "}" + note + ", continued")
        cap = (first + r". Each row gives $\sum_j m_j x_{s_j}=M x_p$ with $\sigma=\{s_1,\dots,s_4\}$ in increasing order, "
               r"positive integer weights $m_j$ with sum $M$, and $|\det\sigma|=gM$ for the $4\times4$ determinant of the rows $(1,x_s)$.")
        lab = r"\label{" + label + "}" if c == 0 else r"\label{" + label + "abcdef"[c] + "}"
        out += [r"\begin{table*}[tp]", r"\scriptsize", r"\caption{" + cap + "}" + lab,
                r"\begin{ruledtabular}", r"\begin{tabular}{rlrrrrrr@{\hspace{2.2em}}rlrrrrrr}",
                head + " & " + head + r" \\ \hline"]
        for x, y in zip(a, b):
            out.append(x + " & " + y + r" \\")
        out += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table*}", ""]
    open(os.path.join(OUT, fname), "w").write("\n".join(out))


def orient(P, S):
    return det([[1] + P[v] for v in S])


def hull_layers(P, d):
    """Convex layers of P (lists of indices) and the facets of each layer with more than d+1 points."""
    def inside(p, S):
        o = orient(P, S)
        for j in range(d + 1):
            T = list(S); T[j] = p
            l = orient(P, T)
            if l == 0 or (l > 0) != (o > 0):
                return False
        return o != 0
    rem, out = list(range(len(P))), []
    while rem:
        v = rem[:] if len(rem) <= d + 1 else [p for p in rem if not any(
            inside(p, S) for S in combinations([q for q in rem if q != p], d + 1))]
        F = []
        if len(v) > d + 1:
            for S in combinations(v, d):
                sg = {orient(P, list(S) + [q]) > 0 for q in v if q not in S}
                if len(sg) == 1:
                    F.append(S)
        out.append((v, F)); rem = [p for p in rem if p not in v]
    return out


def structure_r3(P):
    L = hull_layers(P, 3)
    assert [len(v) for v, _ in L] == [5, 6, 4], [len(v) for v, _ in L]
    (v1, F1), (v2, F2), (v3, _) = L
    assert len(F1) == 6 and len(F2) == 8
    apex = [i for i in v1 if sum(i in f for f in F1) == 3]
    base = [i for i in v1 if i not in apex]
    edges = {frozenset(e) for f in F2 for e in combinations(f, 2)}
    pairs = sorted(tuple(sorted(e)) for e in combinations(v2, 2) if frozenset(e) not in edges)
    assert len(apex) == 2 and len(pairs) == 3
    x = lambda i: f"x_{{{i+1}}}"
    return {"LAYERONEAPEX": f"${x(apex[0])}$ and ${x(apex[1])}$",
            "LAYERONEBASE": "$" + "".join(x(i) for i in base) + "$",
            "LAYERTWOPAIRS": ", ".join(f"${x(a)},{x(b)}$" for a, b in pairs[:-1]) + f" and ${x(pairs[-1][0])},{x(pairs[-1][1])}$",
            "layers": [[i + 1 for i in v] for v, _ in L]}


def main():
    os.makedirs(OUT, exist_ok=True)
    P4 = json.load(open(os.path.join(HD, "witness_r4_11.json")))
    R4 = json.load(open(os.path.join(HD, "cert_r4_11.json")))
    P3 = json.load(open(os.path.join(HD, "witness_r3_15.json")))
    R3 = json.load(open(os.path.join(HD, "cert_r3_15.json")))
    M4, g4 = check(P4, R4, 4, 8)
    M3, g3 = check(P3, R3, 3, 8)
    coords(P4, "coords_r4.tex", r"The eleven points of Theorem~\ref{thm:A} in $\mathbb{R}^4$.", "tab:coordsA")
    coords(P3, "coords_r3.tex", r"The fifteen points of Theorem~\ref{thm:B} in $\mathbb{R}^3$.", "tab:coordsB")
    cert_r4(P4, R4, g4)
    cert_r3(P3, R3, g3)
    mins = []
    for P, d in ((P4, 4), (P3, 3)):
        mins.append(min(abs(det([[1] + P[v] for v in S])) for S in combinations(range(len(P)), d + 1)))
    nums = {"relA": len(R4), "relB": len(R3), "maxMA": f"{M4:,}".replace(",", "{,}"), "maxMB": f"{M3:,}".replace(",", "{,}"),
            "maxcoordA": max(map(max, P4)), "maxcoordB": max(map(max, P3)), "mindetA": mins[0], "mindetB": mins[1],
            "relBinner": sum(1 for r in R3 if r["p"] >= 11)}
    st = structure_r3(P3)
    for key in ("LAYERONEAPEX", "LAYERONEBASE", "LAYERTWOPAIRS"):
        nums[key] = st[key]
    assert st["layers"] == [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11], [12, 13, 14, 15]], st["layers"]
    L4 = hull_layers(P4, 4)
    assert [sorted(v) for v, _ in L4] == [list(range(6)), list(range(6, 11))] and len(L4[0][1]) == 8
    # the nine-point sets of Theorem C
    Q4 = json.load(open(os.path.join(HD, "witness_r4_13_k9.json")))
    S4 = json.load(open(os.path.join(HD, "cert_r4_13_k9.json")))
    Q3 = json.load(open(os.path.join(HD, "witness_r3_18_k9.json")))
    S3 = json.load(open(os.path.join(HD, "cert_r3_18_k9.json")))
    N4, h4 = check(Q4, S4, 4, 9)
    N3, h3 = check(Q3, S3, 3, 9)
    coords(Q4, "coords_r4_k9.tex", r"The thirteen points of Theorem~\ref{thm:C} in $\mathbb{R}^4$.", "tab:coordsCfour")
    coords(Q3, "coords_r3_k9.tex", r"The eighteen points of Theorem~\ref{thm:C} in $\mathbb{R}^3$.", "tab:coordsCthree")
    cert_r4(Q4, S4, h4, k=9, thm="thm:C", label="tab:certCfour", fname="cert_r4_k9.tex", count="forty-six",
            note=r" in $\mathbb{R}^4$", size=r"\footnotesize")
    cert_r3(Q3, S3, h3, per=52, thm="thm:C", label="tab:certCthree", fname="cert_r3_k9.tex", note=r" in $\mathbb{R}^3$")
    mC = [min(abs(det([[1] + P[v] for v in S])) for S in combinations(range(len(P)), d + 1)) for P, d in ((Q4, 4), (Q3, 3))]
    nums.update({"relCfour": len(S4), "relCthree": len(S3), "maxMCfour": f"{N4:,}".replace(",", "{,}"),
                 "maxMCthree": f"{N3:,}".replace(",", "{,}"), "maxcoordCfour": max(map(max, Q4)),
                 "maxcoordCthree": max(map(max, Q3)), "mindetCfour": mC[0], "mindetCthree": mC[1],
                 "relCthreeinner": sum(1 for r in S3 if r["p"] >= 14)})
    # structure stated in Propositions 6.1 and 6.2
    L = hull_layers(Q4, 4)
    assert [sorted(v) for v, _ in L] == [list(range(7)), list(range(7, 13))]
    F = {tuple(f) for f in L[0][1]}
    caps = {(1, 2, 3, 4), (1, 2, 3, 5), (1, 2, 4, 5), (1, 3, 4, 6), (1, 3, 5, 6), (1, 4, 5, 6), (2, 3, 4, 6),
            (2, 3, 5, 7), (2, 3, 6, 7), (2, 4, 5, 7), (2, 4, 6, 7), (3, 5, 6, 7), (4, 5, 6, 7)}
    assert F == {tuple(i - 1 for i in f) for f in caps}
    E = {frozenset(e) for f in L[1][1] for e in combinations(f, 2)}
    assert len(L[1][1]) == 9 and len(E) == 15
    assert all(r["p"] >= 7 for r in S4)
    L = hull_layers(Q3, 3)
    assert [sorted(v) for v, _ in L] == [list(range(7)), list(range(7, 14)), list(range(14, 18))]
    F1 = {tuple(i + 1 for i in f) for f in L[0][1]}
    assert F1 == {(1, 2, 3), (1, 2, 7), (1, 3, 7), (2, 3, 4), (2, 4, 5), (2, 5, 7), (3, 4, 6), (3, 6, 7), (4, 5, 6), (5, 6, 7)}
    F2 = {tuple(i + 1 for i in f) for f in L[1][1]}
    ring = [8, 10, 14, 13, 12]
    assert F2 == {tuple(sorted((a, ring[j], ring[(j + 1) % 5]))) for a in (9, 11) for j in range(5)}
    assert all(r["p"] >= 7 for r in S3)
    with open(os.path.join(OUT, "numbers.tex"), "w") as f:
        for k, v in nums.items():
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    print(nums)


if __name__ == "__main__":
    main()
