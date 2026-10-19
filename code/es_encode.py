#!/usr/bin/env python3
"""CNF encoding of: n points in general position (rank-3 signotope, x-sorted)
with no k points in convex position.  UNSAT for (n,k) proves ES(k) <= n.

Usage: es_encode.py n k [--nosym] > out.cnf
Variables: O(a,b,c), a<b<c, true iff (a,b,c) is a cup (b below segment ac).
"""
import sys, itertools

def build(n, k, sym=True):
    nv = 0
    def new():
        nonlocal nv
        nv += 1
        return nv
    O = {}
    for t in itertools.combinations(range(n), 3):
        O[t] = new()
    cl = []
    # signotope axioms: for a<b<c<d the sequence abc, abd, acd, bcd is monotone
    for a, b, c, d in itertools.combinations(range(n), 4):
        s = [O[(a, b, c)], O[(a, b, d)], O[(a, c, d)], O[(b, c, d)]]
        for i, j, l in itertools.combinations(range(4), 3):
            cl.append([-s[i], s[j], -s[l]])
            cl.append([s[i], -s[j], s[l]])
    if sym:
        # point 0 is extremal, others sorted around it: all (0,b,c) are cups
        for b, c in itertools.combinations(range(1, n), 2):
            cl.append([O[(0, b, c)]])
        # reflection a -> n-a of points 1..n-1 maps O(a,b,c) to O(n-c,n-b,n-a).
        # Lex-leader on the consecutive triples, read from the middle outward
        # (Heule-Scheucher, predicate (20)).
        L, R = [], []
        left = [a for a in range(1, n - 2) if a < n - 2 - a]
        left.sort(reverse=True)
        for a in left:
            L.append(O[(a, a + 1, a + 2)])
            b = n - 2 - a
            R.append(O[(b, b + 1, b + 2)])
        # encode L <=_lex R (false < true)
        eq = None  # literal: prefix equal so far (None = true)
        for x, y in zip(L, R):
            pre = [] if eq is None else [-eq]
            cl.append(pre + [-x, y])          # prefix equal -> x <= y
            e = new()
            # e -> prefix equal and x == y  (only this direction needed)
            if eq is not None:
                cl.append([-e, eq])
            cl.append([-e, -x, y]); cl.append([-e, x, -y])
            # prefix equal and x == y -> e
            cl.append(pre + [x, y, e]); cl.append(pre + [-x, -y, e])
            eq = e
    # chain variables: U[s,m,b,c]: an m-cup starting at s, last two points b<c
    # (m>=3); V likewise for caps. Only the forcing direction is needed.
    U, V = {}, {}
    E, F = {}, {}   # E[s,m,t]: m-cup from s to t ; F caps
    for s in range(n):
        for t in range(s + 1, n):
            for m in range(2, k + 1):
                E[(s, m, t)] = new()
                F[(s, m, t)] = new()
            cl.append([E[(s, 2, t)]])
            cl.append([F[(s, 2, t)]])
    for s in range(n):
        for b in range(s + 1, n):
            for c in range(b + 1, n):
                for m in range(3, k + 1):
                    U[(s, m, b, c)] = new()
                    V[(s, m, b, c)] = new()
    for s in range(n):
        for b in range(s + 1, n):
            for c in range(b + 1, n):
                # 3-chains
                cl.append([-O[(s, b, c)], U[(s, 3, b, c)]])
                cl.append([O[(s, b, c)], V[(s, 3, b, c)]])
                for m in range(3, k + 1):
                    cl.append([-U[(s, m, b, c)], E[(s, m, c)]])
                    cl.append([-V[(s, m, b, c)], F[(s, m, c)]])
                for d in range(c + 1, n):
                    for m in range(3, k):
                        cl.append([-U[(s, m, b, c)], -O[(b, c, d)], U[(s, m + 1, c, d)]])
                        cl.append([-V[(s, m, b, c)], O[(b, c, d)], V[(s, m + 1, c, d)]])
    # no convex k-gon: an i-cap and a j-cup from s to t with i+j = k+2
    for s in range(n):
        for t in range(s + 1, n):
            for j in range(2, k + 1):
                i = k + 2 - j
                if 2 <= i <= k:
                    cl.append([-E[(s, j, t)], -F[(s, i, t)]])
    return nv, cl

if __name__ == "__main__":
    n, k = int(sys.argv[1]), int(sys.argv[2])
    nv, cl = build(n, k, sym="--nosym" not in sys.argv)
    out = sys.stdout
    out.write(f"p cnf {nv} {len(cl)}\n")
    out.write("".join(" ".join(map(str, c)) + " 0\n" for c in cl))
