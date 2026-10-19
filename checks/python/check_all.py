#!/usr/bin/env python3
"""Exact checks of every computational claim of the paper
"New Lower Bounds for Erdos-Szekeres Numbers in Three and Four Dimensions".

All arithmetic is in Python integers. For each of the two point sets
(eleven points in R^4, fifteen points in R^3) the script checks

  (A) the coordinates: integers in [0, C] with C = 40 and C = 180;
  (B) general position: every d+1 points span a full simplex, with the least
      absolute determinant 11 and 6;
  (C) every relation of the certificate: positive weights with gcd 1, the identity
      sum m_j x_{s_j} = M x_p, the stored determinant, and M dividing it;
  (D) the covering condition: every 8-subset contains the support of a relation;
  (E) independently of the certificate, an exact convex-position test of every
      8-subset: no 8 points are in convex position;
  (F) the convex-layer structure of Propositions III.1 and IV.1, and the facts
      about which points occur as p in the relations;
  (G) Schoenheim's bounds 18 and 74 of Remark II.6.

It also checks the convex layers of Scheucher's twelve points in R^3 without a
7-gon (sizes 5, 6, 1: a triangular bipyramid and an octahedron).

Usage: python3 checks/python/check_all.py    (from any directory)
"""
import json
import os
import sys
from itertools import combinations
from math import comb, gcd
from functools import reduce

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
HD = os.path.join(ROOT, "highdim")

passed = failed = 0


def ok(cond, what):
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        print("FAIL:", what)


def det(M):
    """Bareiss fraction-free determinant of a square integer matrix."""
    M = [list(r) for r in M]
    n, sign, prev = len(M), 1, 1
    for i in range(n - 1):
        if M[i][i] == 0:
            for r in range(i + 1, n):
                if M[r][i] != 0:
                    M[i], M[r] = M[r], M[i]
                    sign = -sign
                    break
            else:
                return 0
        for r in range(i + 1, n):
            for c in range(i + 1, n):
                M[r][c] = (M[r][c] * M[i][i] - M[r][i] * M[i][c]) // prev
        prev = M[i][i]
    return sign * M[-1][-1]


def orient(P, S):
    return det([[1] + list(P[v]) for v in S])


def in_open_simplex(P, p, S, o=None):
    """p lies in the interior of conv S (S of d+1 points) iff replacing any vertex by p
    keeps the nonzero orientation sign."""
    o = orient(P, S) if o is None else o
    if o == 0:
        return False
    for j in range(len(S)):
        T = list(S)
        T[j] = p
        t = orient(P, T)
        if t == 0 or (t > 0) != (o > 0):
            return False
    return True


def convex_layers(P, d):
    rem, layers = list(range(len(P))), []
    while rem:
        if len(rem) <= d + 1:
            layers.append(rem[:])
            break
        v = [p for p in rem if not any(in_open_simplex(P, p, S)
                                       for S in combinations([q for q in rem if q != p], d + 1))]
        layers.append(v)
        rem = [p for p in rem if p not in v]
    return layers


def facets(P, V, d):
    """d-subsets of V whose hyperplane leaves all other points of V strictly on one side."""
    F = []
    for S in combinations(V, d):
        sg = {orient(P, list(S) + [q]) > 0 for q in V if q not in S}
        if len(sg) == 1:
            F.append(S)
    return F


def check_set(name, P, rels, d, k, cmax, mindet, nrel):
    n = len(P)
    print(f"== {name}: {n} points in R^{d}, {len(rels)} relations, k = {k}")
    # (A) coordinates
    ok(all(len(x) == d and all(isinstance(c, int) and 0 <= c <= cmax for c in x) for x in P),
       f"{name}: coordinates are integers in [0, {cmax}]")
    ok(max(map(max, P)) == cmax, f"{name}: largest coordinate {cmax}")
    # (B) general position
    dets = {S: orient(P, S) for S in combinations(range(n), d + 1)}
    ok(all(v != 0 for v in dets.values()), f"{name}: general position")
    ok(min(abs(v) for v in dets.values()) == mindet, f"{name}: least |det| = {mindet}")
    # (C) relations
    ok(len(rels) == nrel, f"{name}: {nrel} relations")
    for r in rels:
        p, S, m = r["p"], r["simplex"], r["weights"]
        M = sum(m)
        D = orient(P, S)
        ok(len(S) == d + 1 and len(set(S)) == d + 1 and p not in S and S == sorted(S),
           f"{name}: relation {r} has d+1 distinct sorted vertices other than p")
        ok(all(x > 0 for x in m) and reduce(gcd, m) == 1, f"{name}: weights {m} positive with gcd 1")
        ok(all(sum(mi * P[v][j] for mi, v in zip(m, S)) == M * P[p][j] for j in range(d)),
           f"{name}: identity of relation {r}")
        ok(D == r["orient"] and D != 0 and D % M == 0, f"{name}: determinant of relation {r}")
        ok(in_open_simplex(P, p, S, D), f"{name}: x_{p+1} inside the simplex of relation {r}")
    # (D) covering
    supports = [frozenset([r["p"]] + r["simplex"]) for r in rels]
    subsets = list(combinations(range(n), k))
    ok(len(subsets) == comb(n, k), f"{name}: number of {k}-subsets")
    ok(all(any(s <= frozenset(T) for s in supports) for T in subsets),
       f"{name}: every {k}-subset contains a support")
    # (E) independent convex-position test of every k-subset
    sign = {S: v > 0 for S, v in dets.items()}

    def inside(p, S):
        s0 = sign[S]
        for j in range(d + 1):
            T = tuple(sorted(S[:j] + (p,) + S[j + 1:]))
            # orientation of the tuple with p in position j, from the sorted one
            perm = list(S[:j] + (p,) + S[j + 1:])
            inv = sum(1 for a in range(d + 1) for b in range(a + 1, d + 1) if perm[a] > perm[b])
            if (sign[T] ^ (inv % 2 == 1)) != s0:
                return False
        return True

    def count_convex(size):
        return sum(1 for T in combinations(range(n), size)
                   if not any(inside(p, S) for p in T for S in combinations([q for q in T if q != p], d + 1)))

    convex = count_convex(k)
    ok(convex == 0, f"{name}: no {k} points in convex position (found {convex})")
    # control: the same test finds (k-1)-gons, and agrees with the direct sign test on them
    c7 = count_convex(k - 1)
    direct = sum(1 for T in combinations(range(n), k - 1)
                 if not any(in_open_simplex(P, p, S) for p in T for S in combinations([q for q in T if q != p], d + 1)))
    ok(c7 > 0 and c7 == direct, f"{name}: control, {c7} convex {k-1}-subsets by both tests ({direct})")
    print(f"   {k-1}-gons: {c7}; {k}-gons: {convex}")
    return dets


def main():
    P4 = json.load(open(os.path.join(HD, "witness_r4_11.json")))
    R4 = json.load(open(os.path.join(HD, "cert_r4_11.json")))
    P3 = json.load(open(os.path.join(HD, "witness_r3_15.json")))
    R3 = json.load(open(os.path.join(HD, "cert_r3_15.json")))
    PS = json.load(open(os.path.join(HD, "scheucher_r3_12.json")))

    check_set("X4", P4, R4, 4, 8, 40, 11, 28)
    check_set("X3", P3, R3, 3, 8, 180, 6, 161)

    print("== structure of X4 (Proposition III.1)")
    L4 = convex_layers(P4, 4)
    ok([sorted(v) for v in L4] == [list(range(6)), list(range(6, 11))], "X4: hull vertices x1..x6, x7..x11 interior")
    F4 = facets(P4, list(range(6)), 4)
    want = sorted(tuple(sorted((a,) + t)) for a in (0, 5) for t in combinations(range(1, 5), 3))
    ok(sorted(F4) == want, "X4: the eight facets are x1 or x6 with three of x2..x5")
    ok(all(r["p"] in range(6, 11) for r in R4), "X4: every relation has p among x7..x11")
    ok(max(sum(r["weights"]) for r in R4) == 136733, "X4: largest M = 136733")
    comp = [frozenset(range(11)) - frozenset([r["p"]] + r["simplex"]) for r in R4]
    ok(all(any(frozenset(t) <= c for c in comp) for t in combinations(range(11), 3)),
       "X4: every triple lies in the complement of a support")

    print("== structure of X3 (Proposition IV.1)")
    L3 = convex_layers(P3, 3)
    ok([sorted(v) for v in L3] == [list(range(5)), list(range(5, 11)), list(range(11, 15))], "X3: layers 5, 6, 4")
    F1 = facets(P3, list(range(5)), 3)
    ok(len(F1) == 6 and all(sum(i in f for f in F1) == 3 for i in (1, 2))
       and all(sum(i in f for f in F1) == 4 for i in (0, 3, 4)),
       "X3: outer layer is a bipyramid with apexes x2, x3 over x1x4x5")
    F2 = facets(P3, list(range(5, 11)), 3)
    edges = {frozenset(e) for f in F2 for e in combinations(f, 2)}
    nonedges = sorted(tuple(sorted(e)) for e in combinations(range(5, 11), 2) if frozenset(e) not in edges)
    ok(len(F2) == 8 and nonedges == [(5, 8), (6, 10), (7, 9)],
       "X3: second layer is an octahedron with opposite pairs x6x9, x7x11, x8x10")
    ok(all(orient(P3, list(f) + [q]) != 0 and
           (orient(P3, list(f) + [q]) > 0) == (orient(P3, list(f) + [next(i for i in range(5) if i not in f)]) > 0)
           for f in F1 for q in range(5, 15)), "X3: the inner ten points lie strictly inside the bipyramid")
    ps = {r["p"] for r in R3}
    ok(ps == set(range(5, 15)), "X3: the points p run through x6..x15 and never through x1..x5")
    ok(sum(1 for r in R3 if r["p"] >= 11) == 90, "X3: 90 relations have p in the tetrahedron")
    ok(max(sum(r["weights"]) for r in R3) == 261355, "X3: largest M = 261355")

    print("== Scheucher's twelve points in R^3 without a 7-gon")
    LS = convex_layers(PS, 3)
    ok([len(v) for v in LS] == [5, 6, 1], "Scheucher: convex layers of sizes 5, 6, 1")
    FS1 = facets(PS, LS[0], 3)
    FS2 = facets(PS, LS[1], 3)
    ok(len(FS1) == 6 and sorted(sum(i in f for f in FS1) for i in LS[0]) == [3, 3, 4, 4, 4],
       "Scheucher: outer layer is a triangular bipyramid")
    ok(len(FS2) == 8 and all(sum(i in f for f in FS2) == 4 for i in LS[1]),
       "Scheucher: second layer is an octahedron")

    print("== Schoenheim's bound (Remark II.6)")

    def schoenheim(n, s, t):
        b = 1
        for i in range(t - 1, -1, -1):
            b = -((-(n - i) * b) // (s - i))
        return b
    ok(schoenheim(11, 5, 3) == 18, "Schoenheim bound for (11, 5, 3) is 18")
    ok(schoenheim(15, 10, 7) == 74, "Schoenheim bound for (15, 10, 7) is 74")

    print(f"{passed} checks passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
