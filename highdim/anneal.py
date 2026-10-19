#!/usr/bin/env python3
"""Simulated annealing with reheats for search.py's cost.
Usage: anneal.py n d k seed seconds [init.npy]"""
import sys, time, math
import numpy as np
import search

n, d, k, seed, secs = map(int, sys.argv[1:6])
rng = np.random.default_rng(seed)
P = np.load(sys.argv[6]) if len(sys.argv) > 6 else rng.normal(size=(n, d))
if len(P) < n:   # extend a smaller witness by points near random members
    P = np.vstack([P] + [P[rng.integers(len(P))] + rng.normal(size=d) * 0.05 * P.std() for _ in range(n - len(P))])
c = search.cost(P, d, k); best = c; bestP = P.copy()
T = 2.0; t0 = time.time(); it = 0; last_imp = 0
while time.time() - t0 < secs and best > 0:
    it += 1
    Q = P.copy()
    r = rng.random()
    if r < 0.8:
        i = rng.integers(n); Q[i] += rng.normal(size=d) * rng.choice([0.01, 0.03, 0.1, 0.3])
    elif r < 0.95:
        i = rng.integers(n); Q[i] = Q[rng.integers(n)] + rng.normal(size=d) * 0.2
    else:
        i = rng.integers(n); Q[i] = rng.normal(size=d) * rng.choice([0.5, 1, 3])
    c2 = search.cost(Q, d, k)
    if c2 <= c or rng.random() < math.exp((c - c2) / T):
        P, c = Q, c2
    if c < best:
        best, bestP, last_imp = c, P.copy(), it
        print(f"{time.time()-t0:.0f}s it {it} best {best}", flush=True)
        np.save(f"best_n{n}_d{d}_k{k}_s{seed}.npy", bestP)
    T = max(0.1, T * 0.999)
    if it - last_imp > 4000:          # reheat from the best point set
        P, c, T, last_imp = bestP.copy(), best, 1.5, it
if best == 0:
    print("FOUND", flush=True); print(bestP.tolist(), flush=True)
print(f"# end best {best} it {it}", flush=True)
