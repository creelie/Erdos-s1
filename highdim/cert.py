#!/usr/bin/env python3
"""Hand-checkable certificate that a point set has no k points in convex position
(or, with --hole, no empty convex k-polytope).

A relation is a point p and d+1 points v_0..v_d of the set with positive integers
m_0..m_d such that m_0 v_0 + ... + m_d v_d = (m_0 + ... + m_d) p and
det[1 v_0; ...; 1 v_d] != 0.  Then p lies in the interior of the simplex conv(v_i),
so no set containing p and the v_i is in convex position, and this stays true
after any small enough perturbation of the points.  The script picks, greedily,
a short list of relations such that every k-subset contains one of them (for
--hole: every k-subset T either contains a relation, or contains the v_i of a
relation whose p is outside T).  It prints the list and re-checks it.
Usage: cert.py points.json d k [--hole] [--tex out.tex]"""
import sys, json, itertools
from fractions import Fraction
from math import gcd
sys.path.insert(0, __import__('os').path.dirname(__file__))
from verify_exact import orient

def relations(P, d):
    n = len(P); rel = []
    for S in itertools.combinations(range(n), d + 1):
        o = orient(P, S)
        if o == 0: continue
        for p in range(n):
            if p in S: continue
            lam = []
            for j in range(d + 1):
                T = list(S); T[j] = p
                lam.append(orient(P, T))
            if all((l > 0) == (o > 0) and l != 0 for l in lam):
                g = 0
                for l in lam: g = gcd(g, abs(l))
                m = [abs(l) // g for l in lam]
                rel.append((p, S, m, o))
    return rel

def check(P, d, p, S, m):
    M = sum(m)
    return all(sum(m[j] * P[S[j]][c] for j in range(d + 1)) == M * P[p][c] for c in range(d)) and \
        all(x > 0 for x in m) and orient(P, S) != 0

def cover(P, d, k, hole=False):
    n = len(P); rel = relations(P, d)
    subsets = list(itertools.combinations(range(n), k))
    index = {T: i for i, T in enumerate(subsets)}
    def mask(r):
        # bitmask of the k-subsets this relation rules out
        p, S, _, _ = r
        base = set(S) if hole else set(S) | {p}
        rest = [i for i in range(n) if i not in base]
        b = 0
        for E in itertools.combinations(rest, k - len(base)):
            T = tuple(sorted(base | set(E)))
            if hole or p in T: b |= 1 << index[T]
        return b
    masks = [mask(r) for r in rel]
    left = (1 << len(subsets)) - 1; chosen = []
    while left:
        j = max(range(len(rel)), key=lambda i: (masks[i] & left).bit_count())
        if masks[j] & left == 0:
            return None, [subsets[i] for i in range(len(subsets)) if left >> i & 1]
        chosen.append(rel[j]); left &= ~masks[j]
    return chosen, []

if __name__ == "__main__":
    a = sys.argv
    P = json.load(open(a[1])); d, k = int(a[2]), int(a[3]); hole = "--hole" in a
    chosen, left = cover(P, d, k, hole)
    if chosen is None:
        print(f"no certificate: {len(left)} {k}-subsets are not covered, e.g. {left[0]}"); sys.exit(1)
    n = len(P)
    ok = all(check(P, d, p, S, m) for p, S, m, _ in chosen)
    cov = all(any(set(S) <= set(T) and (p in T or hole) for p, S, m, _ in chosen)
              for T in itertools.combinations(range(n), k))
    print(f"n={n} d={d} k={k}{' holes' if hole else ''}: {len(chosen)} relations; identities {'OK' if ok else 'FAIL'}; cover {'OK' if cov else 'FAIL'}")
    for p, S, m, o in chosen:
        print(f"  {p+1:>2} in ({' '.join(str(s+1) for s in S)})  weights {m}")
    if "--tex" in a:
        out = a[a.index("--tex") + 1]
        with open(out, "w") as f:
            for p, S, m, o in chosen:
                f.write(f"${p+1}$ & ${', '.join(str(s+1) for s in S)}$ & ${', '.join(map(str, m))}$ \\\\\n")
