#!/usr/bin/env python3
"""Unbiased estimate of total solving cost of the ES(k) formula split on all
consecutive triples O(a,a+1,a+2), a = 1..n-3, plus the alternating splits.
Samples base sequences uniformly (DP), one split uniformly, and writes the
cube with its weight N_base * 2^(#alt)."""
import sys, random
import cubes

def count_table(L, tmax, fmax):
    # f[i][(t,f)] = number of valid completions from position i with current runs
    from functools import lru_cache
    @lru_cache(None)
    def f(i, t, fr):
        if i == L:
            return 1
        r = 0
        if t + 1 <= tmax: r += f(i + 1, t + 1, 0)
        if fr + 1 <= fmax: r += f(i + 1, 0, fr + 1)
        return r
    return f

def sample(n, k, rng):
    win = list(range(1, n - 2)); L = len(win)
    tmax, fmax = k - 4, k - 3
    f = count_table(L, tmax, fmax)
    N = f(0, 0, 0)
    bits = []; t = fr = 0
    for i in range(L):
        a = f(i + 1, t + 1, 0) if t + 1 <= tmax else 0
        b = f(i + 1, 0, fr + 1) if fr + 1 <= fmax else 0
        if rng.random() * (a + b) < a:
            bits.append(1); t += 1; fr = 0
        else:
            bits.append(0); fr += 1; t = 0
    lits = [cubes.tvar(n, x, x + 1, x + 2) * (1 if v else -1) for x, v in zip(win, bits)]
    nalt = 0
    for i in range(L - 2):
        if bits[i] and not bits[i + 1] and bits[i + 2]:
            nalt += 1
            v = cubes.tvar(n, win[i], win[i] + 2, win[i] + 4)
            lits.append(v * rng.choice((1, -1)))
    lexok = cubes.lex_ok(n, win, bits)
    return N, nalt, lexok, lits

if __name__ == "__main__":
    n, k, seed = map(int, sys.argv[1:4])
    N, nalt, lexok, lits = sample(n, k, random.Random(seed))
    print(f"c N={N} nalt={nalt} weight={N * 2**nalt} lexok={int(lexok)}")
    for x in lits:
        print(f"{x} 0")
