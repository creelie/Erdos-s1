#!/usr/bin/env python3
"""Relabel a point set layer by layer: first the vertices of the convex hull, then the
vertices of the hull of the rest, and so on (within a layer, by the first coordinate).
Usage: relabel.py points.json d out.json"""
import sys, json, itertools
sys.path.insert(0, __import__('os').path.dirname(__file__))
from verify_exact import orient
def layers(P, d):
    def inside(p, S):
        o = orient(P, S)
        if o == 0: return False
        for j in range(d + 1):
            T = list(S); T[j] = p
            l = orient(P, T)
            if l == 0 or (l > 0) != (o > 0): return False
        return True
    rem = list(range(len(P))); out = []
    while rem:
        if len(rem) <= d + 1: v = rem[:]
        else: v = [p for p in rem if not any(inside(p, S) for S in itertools.combinations([q for q in rem if q != p], d + 1))]
        out.append(sorted(v, key=lambda i: P[i])); rem = [p for p in rem if p not in v]
    return out
if __name__ == "__main__":
    P = json.load(open(sys.argv[1])); d = int(sys.argv[2])
    L = layers(P, d); order = [i for l in L for i in l]
    json.dump([P[i] for i in order], open(sys.argv[3], "w"))
    print("layer sizes", [len(l) for l in L])
