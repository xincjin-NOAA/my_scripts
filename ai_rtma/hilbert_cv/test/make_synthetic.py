"""Synthetic test obs on a 2.5 km CONUS-sized grid: 40 dense clusters + uniform
background, plus 300 duplicates (same location, 0.5 h off the analysis time)."""
import numpy as np
r = np.random.default_rng(0); nx, ny = 1799, 1059
c = r.uniform([100, 100], [nx-100, ny-100], (40, 2))
p = np.vstack([c[r.integers(40, size=15000)] + r.normal(0, 60, (15000, 2)),
               r.uniform([1, 1], [nx, ny], (5000, 2))])
p = np.clip(p, 1, [nx, ny])
dup = p[r.integers(len(p), size=300)]
with open("test/obs_in.txt", "w") as f:
    for i, (x, y) in enumerate(np.vstack([p, dup])):
        t = 0.0 if i < len(p) else 0.5
        f.write(f"S{i:07d} 181 {-125+x/nx*60:.4f} {24+y/ny*26:.4f} {x:.4f} {y:.4f} {t}\n")
