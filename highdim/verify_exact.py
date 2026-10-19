#!/usr/bin/env python3
"""Exact check, in integer arithmetic, that a set of integer points in R^d is
in general position and has no k points in convex position.
Usage: verify_exact.py points.json d k   (points.json: list of integer lists)"""
import sys, json, itertools

def det(M):
    # Bareiss fraction-free determinant of an integer matrix
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for i in range(n - 1):
        if M[i][i] == 0:
            for r in range(i + 1, n):
                if M[r][i] != 0:
                    M[i], M[r] = M[r], M[i]; sign = -sign; break
            else:
                return 0
        for r in range(i + 1, n):
            for c in range(i + 1, n):
                M[r][c] = (M[r][c] * M[i][i] - M[r][i] * M[i][c]) // prev
        prev = M[i][i]
    return sign * M[-1][-1]

def orient(P, idx):
    return det([[1] * len(idx)] + [[P[i][j] for i in idx] for j in range(len(P[0]))])

def main(P, d, k):
    n = len(P)
    sgn = {}
    for S in itertools.combinations(range(n), d + 1):
        o = orient(P, S)
        if o == 0:
            return f"not in general position: {S}"
        sgn[S] = o > 0
    def inside(p, S):
        # p in conv(S) iff replacing each vertex of S by p keeps the orientation sign
        s0 = sgn[S]
        for j in range(d + 1):
            T = list(S); T[j] = p
            o = orient(P, T)
            if (o > 0) != s0:
                return False
        return True
    wit = [(1 << p) | sum(1 << s for s in S)
           for p in range(n) for S in itertools.combinations([i for i in range(n) if i != p], d + 1)
           if inside(p, S)]
    bad = [T for T in itertools.combinations(range(n), k)
           if not any((w & sum(1 << t for t in T)) == w for w in wit)]
    return f"n={n} d={d} k={k}: general position OK; convex {k}-subsets: {len(bad)}" + (f" e.g. {bad[0]}" if bad else "")

if __name__ == "__main__":
    P = json.load(open(sys.argv[1])); d, k = int(sys.argv[2]), int(sys.argv[3])
    print(main(P, d, k))
