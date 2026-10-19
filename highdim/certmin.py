#!/usr/bin/env python3
"""Shortest certificate found by randomized greedy covering with pruning, among the
relations of cert.py whose total weight is at most the least threshold that allows a cover.
Usage: certmin.py points.json d k rounds out.json [--hole]"""
import sys, json, random, itertools
sys.path.insert(0, __import__('os').path.dirname(__file__))
from cert import relations, check
def main(P, d, k, rounds, hole=False, seed=1):
    n = len(P); rel = relations(P, d)
    subsets = list(itertools.combinations(range(n), k)); index = {T: i for i, T in enumerate(subsets)}
    def mask(r):
        p, S, _, _ = r
        base = set(S) if hole else set(S) | {p}
        rest = [i for i in range(n) if i not in base]; b = 0
        for E in itertools.combinations(rest, k - len(base)):
            T = tuple(sorted(base | set(E)))
            if hole or p in T: b |= 1 << index[T]
        return b
    full = (1 << len(subsets)) - 1; allm = [mask(r) for r in rel]
    tot = sorted(set(sum(r[2]) for r in rel))
    def covers(T):
        u = 0
        for r, b in zip(rel, allm):
            if sum(r[2]) <= T: u |= b
        return u == full
    lo, hi = 0, len(tot) - 1
    if not covers(tot[hi]): return None
    while lo < hi:
        mid = (lo + hi) // 2
        if covers(tot[mid]): hi = mid
        else: lo = mid + 1
    idx = [i for i, r in enumerate(rel) if sum(r[2]) <= tot[lo]]
    rng = random.Random(seed); best = None
    for it in range(rounds):
        left = full; ch = []
        while left:
            sc = [(allm[i] & left).bit_count() for i in idx]; top = max(sc)
            cand = [i for i, s in zip(idx, sc) if s >= top * (1 if it == 0 else rng.uniform(0.8, 1.0))]
            j = rng.choice(cand) if it else max(cand, key=lambda i: -sum(rel[i][2]))
            ch.append(j); left &= ~allm[j]
        rng.shuffle(ch)
        for j in list(ch):  # drop redundant relations
            u = 0
            for i in ch:
                if i != j: u |= allm[i]
            if u == full: ch.remove(j)
        if best is None or len(ch) < len(best): best = ch
    return [rel[i] for i in best], tot[lo]
if __name__ == "__main__":
    a = sys.argv; P = json.load(open(a[1])); d, k, rounds = int(a[2]), int(a[3]), int(a[4])
    hole = "--hole" in a
    ch, T = main(P, d, k, rounds, hole)
    ok = all(check(P, d, p, S, m) for p, S, m, _ in ch)
    cov = all(any(set(S) <= set(T_) and (p in T_ or hole) for p, S, m, _ in ch) for T_ in itertools.combinations(range(len(P)), k))
    print(f"{len(ch)} relations, total weight at most {T}; identities {'OK' if ok else 'FAIL'}; cover {'OK' if cov else 'FAIL'}")
    json.dump([{"p": p, "simplex": list(S), "weights": m, "orient": o} for p, S, m, o in sorted(ch, key=lambda r: (r[0], r[1]))], open(a[5], "w"))
