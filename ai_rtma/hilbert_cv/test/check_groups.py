"""Nearest-neighbour spacing within each Hilbert group vs a random subset of the
same size, and where the ungrouped (group 0) obs are."""
import sys
import numpy as np
from scipy.spatial import cKDTree
d = np.genfromtxt(sys.argv[1] if len(sys.argv) > 1 else "test/obs_groups.txt", dtype=None, encoding=None)
xy = np.c_[d["f4"], d["f5"]]; g = d["f7"]
def nn(p):
    return np.sort(cKDTree(p).query(p, k=2)[0][:, 1])
r = np.random.default_rng(1)
ok = g >= 0
print("grid units (2.5 km each): median / 5th percentile nearest-neighbour distance")
for k in range(1, g.max() + 1):
    a = nn(xy[g == k]); b = nn(xy[ok][r.choice(ok.sum(), (g == k).sum(), replace=False)])
    print(f"group {k}: n={len(a)}  hilbert {np.median(a):6.2f} / {np.percentile(a, 5):5.2f}   random {np.median(b):6.2f} / {np.percentile(b, 5):5.2f}")
# local density (obs within 20 grid units) for grouped vs ungrouped obs
t = cKDTree(xy[ok]); dens = np.array([len(v) for v in t.query_ball_point(xy[ok], 20)])
print(f"local density (obs within 50 km): group 0 median {np.median(dens[g[ok]==0]):.0f}, groups 1+ median {np.median(dens[g[ok]>0]):.0f}")
