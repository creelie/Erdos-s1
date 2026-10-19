#!/usr/bin/env python3
"""Solve the cubes of cubes.py incrementally with CaDiCaL (pysat).
Usage: run_cubes.py n k l [shard nshards] [--sample N --seed S --limit SEC]
Prints one line per cube: index, result, seconds, conflicts."""
import sys, time, random, argparse
from pysat.solvers import Solver
import es_encode, cubes

p = argparse.ArgumentParser()
p.add_argument("n", type=int); p.add_argument("k", type=int); p.add_argument("l", type=int)
p.add_argument("--shard", type=int, default=0); p.add_argument("--nshards", type=int, default=1)
p.add_argument("--sample", type=int, default=0); p.add_argument("--seed", type=int, default=1)
p.add_argument("--fresh", action="store_true", help="new solver per cube")
a = p.parse_args()
nv, cl = es_encode.build(a.n, a.k)
cs = cubes.cubes(a.n, a.k, a.l)
idx = list(range(len(cs)))[a.shard::a.nshards]
if a.sample:
    random.seed(a.seed); idx = random.sample(idx, min(a.sample, len(idx)))
print(f"# n={a.n} k={a.k} l={a.l} cubes={len(cs)} here={len(idx)}", flush=True)
s = None
tot = time.time()
for i in idx:
    if s is None or a.fresh:
        if s: s.delete()
        s = Solver(name="cadical153", bootstrap_with=cl)
    t = time.time()
    r = s.solve(assumptions=cs[i])
    st = s.accum_stats()
    print(i, "SAT" if r else "UNSAT", f"{time.time()-t:.2f}", st.get("conflicts", -1), flush=True)
    if r:
        print("# SAT cube found", i, flush=True)
print(f"# total {time.time()-tot:.1f}s", flush=True)
