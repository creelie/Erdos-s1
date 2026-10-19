#!/usr/bin/env python3
"""Search for n points in R^d in general position with no k points in convex
position (lower-bound witnesses for the Erdos-Szekeres numbers g^(d)(k)).

A k-subset is in convex position iff no member lies in the simplex spanned by
d+1 other members (Caratheodory).  We precompute, for every point p and every
(d+1)-subset S not containing p, whether p lies in conv(S), and count convex
k-subsets.  Simulated annealing on the coordinates drives the count to 0.

Usage: search.py n d k seed [iters]
"""
import sys, itertools, math, random
import numpy as np

def inside_table(P, d):
    n = len(P)
    subs = list(itertools.combinations(range(n), d + 1))
    A = np.ones((len(subs), d + 1, d + 1))
    for i, S in enumerate(subs):
        A[i, 1:, :] = P[list(S)].T
    det0 = np.linalg.det(A)
    ins = {}
    for p in range(n):
        sel = [i for i, S in enumerate(subs) if p not in S]
        # barycentric signs: replace column j by p
        ok = np.ones(len(sel), dtype=bool)
        B = A[sel].copy()
        for j in range(d + 1):
            C = B.copy()
            C[:, 1:, j] = P[p]
            ok &= np.sign(np.linalg.det(C)) == np.sign(det0[sel])
        for i, o in zip(sel, ok):
            if o:
                ins.setdefault(p, []).append(subs[i])
    return ins

def convex_subsets(P, d, k, ins=None):
    n = len(P)
    if ins is None:
        ins = inside_table(P, d)
    witness = [(p, set(S)) for p in ins for S in ins[p]]
    # bitmasks
    wm = [(1 << p) | sum(1 << s for s in S) for p, S in witness]
    cnt = 0
    for T in itertools.combinations(range(n), k):
        m = sum(1 << t for t in T)
        if not any((w & m) == w for w in wm):
            cnt += 1
    return cnt

def cost(P, d, k):
    return convex_subsets(P, d, k)

if __name__ == "__main__":
    n, d, k, seed = map(int, sys.argv[1:5])
    iters = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
    rng = np.random.default_rng(seed)
    P = rng.normal(size=(n, d))
    c = cost(P, d, k); best = c; T = 1.0
    for it in range(iters):
        i = rng.integers(n)
        Q = P.copy()
        scale = rng.choice([0.03, 0.1, 0.3, 1.0])
        Q[i] += rng.normal(size=d) * scale
        if rng.random() < 0.05:
            Q[i] = rng.normal(size=d) * rng.choice([0.2, 1, 3])
        c2 = cost(Q, d, k)
        if c2 <= c or rng.random() < math.exp((c - c2) / T):
            P, c = Q, c2
        T = max(0.05, T * 0.9995)
        if c < best:
            best = c
            print(f"it {it} cost {c}", flush=True)
        if c == 0:
            np.save(f"witness_n{n}_d{d}_k{k}_s{seed}.npy", P)
            print("FOUND", flush=True)
            print(repr(P.tolist()))
            break
