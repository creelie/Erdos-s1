#!/usr/bin/env python3
"""Estimate total CPU cost of the fully split ES(k) formula on n points.
Usage: estimate.py n k samples seed [limit_s]"""
import sys, time, random, threading
from pysat.solvers import Solver
import es_encode, sample_full

n, k, M, seed = map(int, sys.argv[1:5])
limit = float(sys.argv[5]) if len(sys.argv) > 5 else 600
nv, cl = es_encode.build(n, k)
rng = random.Random(seed)
s = Solver(name="cadical195", bootstrap_with=cl, use_timer=True)
tot = 0.0; done = 0; timeouts = 0
for j in range(M):
    N, nalt, lexok, lits = sample_full.sample(n, k, rng)
    w = N * 2 ** nalt
    if not lexok:
        print(j, "lex", 0, w, flush=True); done += 1; continue
    t = time.time()
    timer = threading.Timer(limit, s.interrupt); timer.start()
    r = s.solve_limited(assumptions=lits, expect_interrupt=True)
    timer.cancel()
    try: s.clear_interrupt()
    except NotImplementedError: pass
    dt = time.time() - t
    tag = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[r]
    timeouts += r is None
    tot += w * dt; done += 1
    print(j, tag, f"{dt:.3f}", w, flush=True)
print(f"# estimate {tot/done:.4g} CPU-s (samples {done}, timeouts {timeouts})", flush=True)
