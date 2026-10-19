#!/usr/bin/env python3
"""Figures of the paper, written to paper/figures/ as compact PNG files.

fig_layers.png    the fifteen points in R^3 and their three convex layers
fig_fourd.png     the eleven points in R^4, seen through a projection to R^3
fig_relation.png  an interior relation and the perturbations it survives
fig_bounds.png    known values and bounds of g^(d)(k) for d = 3, 4, 5
fig_ninethree.png the eighteen points in R^3 without nine in convex position, by layers
fig_ninefour.png  the thirteen points in R^4 without nine in convex position, projected to R^3

Needs Python 3 with numpy and matplotlib; ImageMagick's convert, if present, reduces
each PNG to a 256-colour palette."""
import json, itertools, os, shutil, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 11, "mathtext.fontset": "cm"})
BLUE, ORANGE, GREEN, RED, GREY = "#2b5c8a", "#d9822b", "#3a8a4f", "#b8312f", "#6b6b6b"

def load(name):
    return np.array(json.load(open(os.path.join(ROOT, "highdim", name))), dtype=float)

def orient(P, S):
    M = np.array([[1.0, *P[i]] for i in S])
    return np.linalg.det(M)

def facets(P, H, d):
    out = []
    for S in itertools.combinations(H, d):
        s = {orient(P, list(S) + [q]) > 0 for q in H if q not in S}
        if len(s) == 1: out.append(S)
    return out

def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=290, facecolor="white", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    if shutil.which("convert"):
        subprocess.run(["convert", path, "-background", "white", "-alpha", "remove", "-alpha", "off",
                        "-resize", "1600x1600>", "-colors", "256", "PNG8:" + path], check=True)
    print("wrote", path)

def polytope(ax, P, idx, color, alpha, lw=1.2, faces=True, d=3):
    F = facets(P, idx, d)
    if faces:
        ax.add_collection3d(Poly3DCollection([P[list(f)] for f in F], facecolor=color, alpha=alpha, edgecolor="none"))
    E = {tuple(sorted(e)) for f in F for e in itertools.combinations(f, 2)}
    ax.add_collection3d(Line3DCollection([P[list(e)] for e in E], colors=color, linewidths=lw))

def affine_fit(X, Y):
    """The affine map x -> A x + b that sends the rows of X closest to those of Y."""
    Xh = np.column_stack([X, np.ones(len(X))])
    T = np.linalg.lstsq(Xh, Y, rcond=None)[0]
    return lambda Z: np.column_stack([Z, np.ones(len(Z))]) @ T

def labels(ax, Q, offs, color="black", size=9):
    cc = Q.mean(axis=0)
    for i, q in enumerate(Q):
        v = q - cc; v = v / (np.linalg.norm(v) + 1e-9)
        ax.text(*(q + v * offs[i]), str(i + 1), fontsize=size, color=color, ha="center", va="center")

def screen_labels(fig, ax, Q, dist=11, size=10):
    """Labels 1, 2, ... placed in screen coordinates, pushed away from the centre of the picture."""
    from mpl_toolkits.mplot3d import proj3d
    fig.canvas.draw()
    S = np.array([proj3d.proj_transform(*q, ax.get_proj())[:2] for q in Q])
    tr = ax.transData.transform(S)
    c = tr.mean(axis=0)
    for i, (s, t) in enumerate(zip(S, tr)):
        v = t - c; v = v / (np.linalg.norm(v) + 1e-9)
        ax.annotate(str(i + 1), tuple(s), xycoords="data", xytext=tuple(dist * v), textcoords="offset points",
                    ha="center", va="center", fontsize=size, zorder=30)

def fig_layers():
    """Schematic: each convex layer of the set of Table II drawn as a regular solid of the same
    combinatorial type, nested about a common centre (the labels are those of Table II)."""
    P = load("witness_r3_15.json")
    Q = np.zeros((15, 3))
    for j, i in enumerate([0, 3, 4]):            # the equator triangle 1, 4, 5 of the bipyramid
        t = 2 * np.pi / 3 + 2 * np.pi * j / 3
        Q[i] = (np.cos(t), np.sin(t), 0)
    Q[1], Q[2] = (0, 0, 1.25), (0, 0, -1.25)     # its apexes 2 and 3
    # octahedron: opposite vertices are the pairs that share no facet
    F = facets(P, list(range(5, 11)), 3)
    opp = {i: [j for j in range(5, 11) if j != i and not any(i in f and j in f for f in F)][0] for i in range(5, 11)}
    R = np.linalg.qr(np.array([[1.0, 0.35, 0.2], [-0.3, 1.0, 0.25], [-0.15, -0.3, 1.0]]))[0]
    done, axes = set(), iter(np.eye(3))
    for i in range(5, 11):
        if i in done: continue
        e = R @ next(axes) * 0.56
        Q[i], Q[opp[i]] = e, -e; done |= {i, opp[i]}
    T = np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]) / np.sqrt(3) * 0.27
    Q[11:15] = T @ np.linalg.qr(np.array([[1.0, 0.5, 0.1], [0.2, 1.0, 0.6], [0.3, 0.1, 1.0]]))[0]
    fig = plt.figure(figsize=(8.0, 7.2)); ax = fig.add_subplot(111, projection="3d")
    polytope(ax, Q, range(0, 5), BLUE, 0.05, 1.0)
    polytope(ax, Q, range(5, 11), ORANGE, 0.13, 1.2)
    polytope(ax, Q, range(11, 15), GREEN, 0.40, 1.4)
    for rng, col, lab in [(range(0, 5), BLUE, "bipyramid, points 1 to 5"),
                          (range(5, 11), ORANGE, "octahedron, points 6 to 11"),
                          (range(11, 15), GREEN, "tetrahedron, points 12 to 15")]:
        ax.scatter(*Q[list(rng)].T, s=46, color=col, depthshade=False, label=lab, edgecolors="black", linewidths=0.4)
    for i, q in enumerate(Q):
        v = q / (np.linalg.norm(q) + 1e-9)
        ax.text(*(q + v * (0.15 if i < 5 else 0.12 if i < 11 else 0.10)), str(i + 1), fontsize=10, ha="center", va="center")
    ax.view_init(elev=17, azim=-38)
    ax.set_box_aspect(np.ptp(Q, axis=0)); ax.set_axis_off()
    ax.legend(loc="lower left", frameon=False, fontsize=11, bbox_to_anchor=(0.0, 0.02))
    save(fig, "fig_layers.png")

def fig_fourd():
    P = load("witness_r4_11.json")
    # affine image with the hull close to a regular bipyramid: the tetrahedron 2, 3, 4, 5 in the
    # hyperplane x_4 = 0 and the apexes 1 and 6 on the axis; then projection to R^3 with the fourth
    # coordinate added along a slanted direction
    T = np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], dtype=float)
    Y = np.zeros((6, 4)); Y[1:5, :3] = T; Y[0, 3] = 1.4; Y[5, 3] = -1.4
    Q4 = affine_fit(P[:6], Y)(P)
    c = Q4[6:].mean(axis=0)
    for i in range(6, 11): Q4[i] = c + 1.7 * (Q4[i] - c)
    # an orthonormal projection to R^3 and a view chosen to keep the thirteen points apart on the page
    M = np.array([[-0.8342, 0.1972, -0.2179], [-0.5487, -0.3011, 0.2368],
                  [0.0436, -0.5831, -0.8095], [0.0339, 0.7283, -0.4911]])
    Q = Q4 @ M
    fig = plt.figure(figsize=(8.0, 7.0)); ax = fig.add_subplot(111, projection="3d")
    F = facets(P, list(range(6)), 4)
    E = {tuple(sorted(e)) for f in F for e in itertools.combinations(f, 2)}
    ax.add_collection3d(Line3DCollection([Q[list(e)] for e in E], colors=BLUE, linewidths=1.0, alpha=0.9))
    I = list(range(6, 11))
    ax.add_collection3d(Poly3DCollection([Q[list(t)] for t in itertools.combinations(I, 3)], facecolor=ORANGE, alpha=0.08, edgecolor="none"))
    ax.add_collection3d(Line3DCollection([Q[list(e)] for e in itertools.combinations(I, 2)], colors=ORANGE, linewidths=1.4))
    ax.scatter(*Q[:6].T, s=48, color=BLUE, depthshade=False, edgecolors="black", linewidths=0.4,
               label="hull vertices 1 to 6 (apexes 1 and 6)")
    ax.scatter(*Q[I].T, s=48, color=ORANGE, depthshade=False, edgecolors="black", linewidths=0.4,
               label="inner 4-simplex, points 7 to 11")
    labels(ax, Q, [0.22] * 6 + [0.16] * 5)
    ax.view_init(elev=22, azim=-35); ax.set_box_aspect(np.ptp(Q, axis=0)); ax.set_axis_off()
    ax.legend(loc="upper left", frameon=False, fontsize=11)
    save(fig, "fig_fourd.png")

def fig_relation():
    C = json.load(open(os.path.join(ROOT, "highdim", "cert_r3_15.json")))
    # a relation with small, not too unequal weights
    r = min(sorted(C, key=lambda r: sum(r["weights"]))[:40], key=lambda r: max(r["weights"]) / min(r["weights"]))
    p, S, m = r["p"], r["simplex"], r["weights"]; M = sum(m)
    # drawn after the affine map sending the simplex to a regular tetrahedron, which keeps the
    # barycentric coordinates m_j / M of the inner point
    V = np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], dtype=float)
    x = sum(mj * v for mj, v in zip(m, V)) / M
    fig = plt.figure(figsize=(8.0, 6.6)); ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection([V[list(f)] for f in itertools.combinations(range(4), 3)],
                                         facecolor=BLUE, alpha=0.08, edgecolor=BLUE, linewidths=1.2))
    rad = 0.10
    u, w = np.mgrid[0:2 * np.pi:28j, 0:np.pi:14j]
    for y in list(V) + [x]:
        ax.plot_surface(y[0] + rad * np.cos(u) * np.sin(w), y[1] + rad * np.sin(u) * np.sin(w),
                        y[2] + rad * np.cos(w), color=GREY, alpha=0.18, linewidth=0, shade=False)
    ax.scatter(*V.T, s=50, color=BLUE, depthshade=False, edgecolors="black", linewidths=0.4)
    ax.scatter(*x, s=70, color=RED, depthshade=False, edgecolors="black", linewidths=0.5)
    for j in range(4):
        ax.plot(*np.column_stack([x, V[j]]), color=RED, lw=0.9, ls=(0, (3, 2)))
        mid = x + 0.5 * (V[j] - x) + np.array([0, 0, 0.12])
        ax.text(*mid, f"${m[j]}/{M}$", fontsize=12, color=RED, ha="center")
    for j, i in enumerate(S):
        ax.text(*(V[j] * 1.22), f"$x_{{{i + 1}}}$", fontsize=14, ha="center", va="center")
    ax.text(*(x + np.array([0.0, 0.0, -0.28])), f"$x_{{{p + 1}}}$", fontsize=14, color=RED, ha="center")
    ax.view_init(elev=24, azim=-20); ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
    save(fig, "fig_relation.png")

def fig_bounds():
    # (d, k): (lower, upper, new); upper None means no small upper bound is recorded
    data = {(3, 5): (6, 6, 0), (3, 6): (9, 9, 0), (3, 7): (13, 13, 0), (3, 8): (16, None, 1), (3, 9): (19, None, 1),
            (4, 5): (5, 5, 0), (4, 6): (7, 7, 0), (4, 7): (9, 9, 0), (4, 8): (12, 13, 1), (4, 9): (14, None, 1),
            (5, 5): (5, 5, 0), (5, 6): (6, 6, 0), (5, 7): (8, 8, 0), (5, 8): (10, 10, 0), (5, 9): (12, 13, 0)}
    old = {(3, 8): 14, (3, 9): 15, (4, 8): 11, (4, 9): 13}
    fig = plt.figure(figsize=(8.6, 6.4)); ax = fig.add_subplot(111, projection="3d")
    dx, dy = 0.40, 0.40
    labels = []
    for (d, k), (lo, hi, new) in sorted(data.items(), key=lambda t: (t[0][0], -t[0][1])):
        x, y = k - dx / 2, 1.6 * (5 - d) - dy / 2
        base = old.get((d, k), lo)
        ax.bar3d(x, y, 0, dx, dy, base, color="#9db8d3", alpha=1.0, shade=True, edgecolor="#4a6a8a", linewidth=0.3)
        if new:
            ax.bar3d(x, y, base, dx, dy, lo - base, color=RED, alpha=1.0, shade=True, edgecolor="#6a1a1a", linewidth=0.3)
        if hi is not None and hi > lo:
            ax.bar3d(x, y, lo, dx, dy, hi - lo, color=ORANGE, alpha=0.30, shade=False, edgecolor=ORANGE, linewidth=0.6)
        top = hi if hi is not None else lo
        lab = str(lo) if hi == lo else (f"{lo}\u2013{hi}" if hi else f"$\\geq${lo}")
        labels.append((k, 1.6 * (5 - d), top, lab, new))
    ax.set_xticks(range(5, 10)); ax.set_yticks([0, 1.6, 3.2]); ax.set_yticklabels(["5", "4", "3"])
    ax.set_zticks([0, 5, 10, 15, 20]); ax.set_ylim(-0.5, 3.7); ax.set_xlim(4.5, 9.5)
    ax.set_xlabel("$k$"); ax.set_ylabel("$d$"); ax.set_zlabel("$g^{(d)}(k)$")
    ax.view_init(elev=30, azim=-72); ax.set_box_aspect((1.7, 1.25, 0.9))
    # labels are placed in screen coordinates so that no bar can hide them
    from mpl_toolkits.mplot3d import proj3d
    fig.canvas.draw()
    for k, y, top, lab, new in labels:
        u, v, _ = proj3d.proj_transform(k, y, top + 0.35, ax.get_proj())
        ax.annotate(lab, (u, v), xycoords="data", xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                    fontsize=9.5, color=RED if new else "black", zorder=20,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#9db8d3", label="known before"), Patch(color=RED, label="gain from the new sets"),
                       Patch(color=ORANGE, alpha=0.35, label="range left open")],
              loc="upper left", frameon=False, fontsize=10)
    save(fig, "fig_bounds.png")

def whiten(X):
    """An affine image of the rows of X with centre 0 and identity covariance."""
    c = X.mean(axis=0)
    w, V = np.linalg.eigh(np.cov((X - c).T))
    return (X - c) @ V / np.sqrt(w)

def fig_ninethree():
    """Schematic: each convex layer of the eighteen points of Theorem C in R^3 under its own affine map
    to a round position, nested about a common centre; the combinatorics are those of the set."""
    P = load("witness_r3_18_k9.json")
    lay = [list(range(0, 7)), list(range(7, 14)), list(range(14, 18))]
    Q = np.zeros((18, 3))
    R = np.linalg.qr(np.array([[1.0, 0.4, 0.1], [-0.2, 1.0, 0.3], [0.1, -0.4, 1.0]]))[0]
    for L, sc in zip(lay, (1.0, 0.56, 0.27)):
        W = whiten(P[L])
        Q[L] = W / np.max(np.linalg.norm(W, axis=1)) * sc
    Q[7:14] = Q[7:14] @ R
    fig = plt.figure(figsize=(8.0, 7.2)); ax = fig.add_subplot(111, projection="3d")
    polytope(ax, Q, lay[0], BLUE, 0.05, 1.0)
    polytope(ax, Q, lay[1], ORANGE, 0.13, 1.2)
    polytope(ax, Q, lay[2], GREEN, 0.40, 1.4)
    for L, col, lab in [(lay[0], BLUE, "outer layer, points 1 to 7"),
                        (lay[1], ORANGE, "second layer, points 8 to 14"),
                        (lay[2], GREEN, "tetrahedron, points 15 to 18")]:
        ax.scatter(*Q[L].T, s=46, color=col, depthshade=False, label=lab, edgecolors="black", linewidths=0.4)
    ax.view_init(elev=17, azim=-38)
    ax.set_box_aspect(np.ptp(Q, axis=0)); ax.set_axis_off()
    screen_labels(fig, ax, Q)
    ax.legend(loc="lower left", frameon=False, fontsize=11, bbox_to_anchor=(0.0, 0.0))
    save(fig, "fig_ninethree.png")

def fig_ninefour():
    """The thirteen points of Theorem C in R^4 after an affine map that makes the hull round, projected
    orthogonally to R^3; the inner points are drawn
    enlarged about their centre."""
    P = load("witness_r4_13_k9.json")
    W = whiten(P[:7]); A = affine_fit(P[:7], W)
    Q4 = A(P)
    c = Q4[7:].mean(axis=0)
    for i in range(7, 13): Q4[i] = c + 2.2 * (Q4[i] - c)
    # an orthonormal projection to R^3 and a view chosen to keep the thirteen points apart on the page
    M = np.array([[-0.8342, 0.1972, -0.2179], [-0.5487, -0.3011, 0.2368],
                  [0.0436, -0.5831, -0.8095], [0.0339, 0.7283, -0.4911]])
    Q = Q4 @ M
    fig = plt.figure(figsize=(8.0, 7.0)); ax = fig.add_subplot(111, projection="3d")
    F = facets(P, list(range(7)), 4)
    E = {tuple(sorted(e)) for f in F for e in itertools.combinations(f, 2)}
    ax.add_collection3d(Line3DCollection([Q[list(e)] for e in E], colors=BLUE, linewidths=1.0, alpha=0.9))
    I = list(range(7, 13))
    FI = facets(P, I, 4)
    ax.add_collection3d(Poly3DCollection([Q[list(t)] for f in FI for t in itertools.combinations(f, 3)],
                                         facecolor=ORANGE, alpha=0.05, edgecolor="none"))
    ax.add_collection3d(Line3DCollection([Q[list(e)] for e in itertools.combinations(I, 2)], colors=ORANGE, linewidths=1.3))
    ax.scatter(*Q[:7].T, s=48, color=BLUE, depthshade=False, edgecolors="black", linewidths=0.4,
               label="hull vertices 1 to 7")
    ax.scatter(*Q[I].T, s=48, color=ORANGE, depthshade=False, edgecolors="black", linewidths=0.4,
               label="inner cyclic polytope, points 8 to 13")
    ax.view_init(elev=17.5, azim=-121.5); ax.set_box_aspect(np.ptp(Q, axis=0)); ax.set_axis_off()
    screen_labels(fig, ax, Q)
    ax.legend(loc="upper left", frameon=False, fontsize=11)
    save(fig, "fig_ninefour.png")

if __name__ == "__main__":
    fig_layers(); fig_fourd(); fig_relation(); fig_bounds(); fig_ninethree(); fig_ninefour()
