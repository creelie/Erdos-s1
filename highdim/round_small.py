#!/usr/bin/env python3
"""Round a point set to the coarsest integer grid on which it stays in general
position and keeps a certificate (cert.py) that no k points are in convex position.
Usage: round_small.py points.txt d k out.json [scales]   (scales: comma separated)"""
import sys, json, random
sys.path.insert(0, __import__('os').path.dirname(__file__))
from cert import cover
from verify_exact import main as check
TRIES = int(__import__("os").environ.get("TRIES", "20"))
def attempt(X, d, k, s, shift):
    lo = [min(x[j] for x in X) for j in range(d)]
    span = max(max(x[j] for x in X) - lo[j] for j in range(d))
    P = [[round((x[j] - lo[j]) / span * s + shift[j]) for j in range(d)] for x in X]
    m = [min(p[j] for p in P) for j in range(d)]
    return [[p[j] - m[j] for j in range(d)] for p in P]
if __name__ == "__main__":
    a = sys.argv
    X = [list(map(float, l.split())) for l in open(a[1]) if l.strip()]
    d, k, out = int(a[2]), int(a[3]), a[4]
    scales = [int(t) for t in a[5].split(",")] if len(a) > 5 else [10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500, 1000]
    random.seed(1)
    for s in scales:
        for t in range(TRIES):
            shift = [0.0] * d if t == 0 else [random.random() for _ in range(d)]
            P = attempt(X, d, k, s, shift)
            if len(set(map(tuple, P))) < len(P): continue
            r = check(P, d, k)
            if "general position OK" not in r or not r.endswith(": 0"): continue
            ch, left = cover(P, d, k)
            if ch:
                json.dump(P, open(out, "w"))
                print(f"scale {s} (try {t}): {r}; certificate with {len(ch)} relations; max coordinate {max(max(p) for p in P)}")
                sys.exit(0)
        print(f"scale {s}: no success in {TRIES} tries")
