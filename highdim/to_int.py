#!/usr/bin/env python3
"""Round a floating point set (one point per line) to integer coordinates,
increasing the scale until the exact check reports general position and no
k points in convex position (or, with --hole, no empty convex k-polytope).
Usage: to_int.py points.txt d k out.json [--hole]"""
import sys, json
sys.path.insert(0, __import__('os').path.dirname(__file__))
from verify_exact import main as check
def run(fn, d, k, out, hole=False):
    X = [list(map(float, l.split())) for l in open(fn) if l.strip()]
    lo = [min(x[j] for x in X) for j in range(d)]; hi = [max(x[j] for x in X) for j in range(d)]
    span = max(hi[j] - lo[j] for j in range(d))
    for scale in [10**e for e in range(2, 9)]:
        P = [[round((x[j] - lo[j]) / span * scale) for j in range(d)] for x in X]
        r = check(P, d, k, hole=hole)
        print(scale, r)
        if r.endswith(": 0"):
            json.dump(P, open(out, "w")); return P
    return None
if __name__ == "__main__":
    a = sys.argv
    run(a[1], int(a[2]), int(a[3]), a[4], "--hole" in a)
