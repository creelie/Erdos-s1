#!/usr/bin/env python3
"""Cubes for the ES(k) formula of es_encode.py (Heule-Scheucher split).

Assign the consecutive triples O(a,a+1,a+2) for a window of length l around
the middle.  Since all (0,a,b) are cups, r consecutive true values give an
(r+3)-cup and r consecutive false values give an (r+2)-cap, so runs of k-3
trues or k-2 falses are excluded.  Each occurrence of the alternating pattern
+,-,+ starting at b is further split on O(b,b+2,b+4).
Prints one cube per line as DIMACS literals of es_encode variable numbering.
"""
import sys, itertools

def tvar(n, a, b, c):
    # variable number of O(a,b,c) in es_encode (combinations order, 1-based)
    idx = 0
    # rank of (a,b,c) among 3-combinations of range(n) in lexicographic order
    from math import comb
    r = 0
    for x in range(a):
        r += comb(n - 1 - x, 2)
    for y in range(a + 1, b):
        r += n - 1 - y
    r += c - b - 1
    return r + 1

def lex_ok(n, win, bits):
    """Mirror symmetry (es_encode): reading pairs (a, n-2-a) from the middle
    outward, the first unequal pair must have O(a)=false, O(n-2-a)=true."""
    val = dict(zip(win, bits))
    for a in sorted((a for a in win if a < n - 2 - a), reverse=True):
        b = n - 2 - a
        if b not in val:
            break
        if val[a] != val[b]:
            return val[a] < val[b]
    return True

def cubes(n, k, l, split_alt=True):
    starts = list(range(1, n - 2))          # a = 1..n-3 (0 is fixed)
    mid = (len(starts) - l) // 2
    win = starts[mid:mid + l]
    tmax, fmax = k - 4, k - 3               # allowed run lengths
    out = []
    for bits in itertools.product((1, 0), repeat=l):
        ok = True; rt = rf = 0
        for x in bits:
            rt = rt + 1 if x else 0
            rf = 0 if x else rf + 1
            if rt > tmax or rf > fmax:
                ok = False; break
        if not ok:
            continue
        if not lex_ok(n, win, bits):
            continue
        lits = [tvar(n, a, a + 1, a + 2) * (1 if x else -1) for a, x in zip(win, bits)]
        extra = [[]]
        if split_alt:
            for i in range(l - 2):
                if bits[i] and not bits[i + 1] and bits[i + 2]:
                    b = win[i]
                    v = tvar(n, b, b + 2, b + 4)
                    extra = [e + [s * v] for e in extra for s in (1, -1)]
        for e in extra:
            out.append(lits + e)
    return out

if __name__ == "__main__":
    n, k, l = map(int, sys.argv[1:4])
    for c in cubes(n, k, l):
        print(" ".join(map(str, c)))
