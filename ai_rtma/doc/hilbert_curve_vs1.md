

# Findings (2026-11-18)

GSI checkout: `osw_gsi` at 85ff00e16 (2025-12-30). File paths below are under `osw_gsi/src/gsi/`.

## Summary

- `hilbert_curve=.true.` (namelist `&SETUP`) splits the surface obs of each ob type into N spatially spread validation groups. The obs of one group are withheld from the 2DVar analysis, so they can be used to verify it. It works only with `twodvar_regional=.true.` (RTMA/URMA).
- GSI has two Hilbert-curve features. They share R. J. Purser's Hilbert-curve library but do different things:
  1. **Cross-validation groups** (`hilbert_curve` option): `hilbert_curve.f90`, `phil.f90`, `phil1.f90`, and the `hilbertcurve` module in `support_2dvar.f90`.
  2. **Ob-density down-weighting** of aircraft winds (`denest` in `phil2.f90`, which uses `phil0.f90`). This one doesn't use the `hilbert_curve` option.
- **Extracting it into a standalone Fortran program works.** The core needs only about 800 lines, with no MPI, BUFR or GSI state. A working version is in `../hilbert_cv/`.

## What the groups mean

Each group is one set of obs that can be withheld from the analysis and used to check it. This is cross-validation: the analysis is scored against obs it was not allowed to use.

| Group (last column of `test/obs_groups.txt`) | Meaning |
|---|---|
| 1 … N (5 in the test) | A validation set. In a GSI run, the obs of **one** of these groups are withheld: their usage is raised so they aren't assimilated. Afterwards the analysis is compared with them. |
| 0 | Not in any validation set, so always assimilated. These are mostly obs in dense clusters. |
| −1 | A duplicate of another ob at the same grid point, removed before grouping. |

**Why there are several groups.** One group holds only about 1/N of the obs that can be withheld. Run the analysis N times, withholding a different group each time, and every ob in groups 1–N is checked exactly once while 80% of those obs are still assimilated each run. This is N-fold cross-validation.

GSI itself withholds only one group per run:
- normally the group set in convinfo (`ncgroup`);
- with `random_cvgrp`, a group chosen from the date, so it changes from cycle to cycle.

In the test, withholding group 3 would hold out 2,118 obs and assimilate the other 17,882.

**Why the groups look the way they do.** Obs are dealt into the groups in turn along the Hilbert curve, and an ob can join a group only if it isn't too close to that group's previous ob. So:
- **Each group covers the whole domain fairly evenly.** The check then samples everywhere, not just one region.
- **Obs in the same group are rarely close together.** Withholding a group never removes a whole local cluster: nearby obs from the other groups and group 0 are still assimilated. So a withheld ob checks how well the analysis fills a small gap, not a large empty area.
- **Group 0 holds the overflow in dense areas.** Where obs are too close for any group to take them, they stay in group 0. In the test, group-0 obs have about twice as many neighbours within 50 km as grouped obs (27 vs 14).


## How the cross-validation works

The call chain (example: `read_prepbufr.f90`):

1. `lhilbert = twodvar_regional .and. hilbert_curve` (`read_prepbufr.f90:617`). `gsimod.F90:2274` turns it off for `oneobtest`.
2. `init_hilbertcurve(maxobs)` (`support_2dvar.f90:2976`) reads namelist `&parmcardhcurve` from the file `parmcard_input`. If the file is missing, it runs with the defaults:
   - `random_cvgrp=.false.`, `usagecv=3`, `usage_dup=3`, `ngrpsfact0=10`;
   - `ngrps_*ob`, the number of groups per ob type: for example t=5, uv=8, ps=5, q=5, gust/vis/mxtm/mitm/pmsl/howv/tcamt/lcbas/cldch=8;
   - `wspd10m`, `uwnd10m` and `vwnd10m` are always set equal to `ngrps_uvob`.
3. `accum_hilbertcurve` is called for every ob read (`support_2dvar.f90:3106`). It keeps an ob only when:
   - `kx` is in 180–188, 280–288, 192–195 or 292–295 (surface land, marine and mesonet), and
   - `usage < 6`.

   It stores the station id, provider, subprovider, grid x/y (`dlon`, `dlat`), earth lon/lat, time offset and index.
4. After the file has been read, `apply_hilbertcurve(ndata, obstype, cdata_all(usage_row,:))` runs (`support_2dvar.f90:3171`, called at `read_prepbufr.f90:3136`):
   1. **Duplicate removal.** Obs at the same grid point (|Δx|, |Δy| < 1e-3 grid units) are duplicates; only the one nearest the analysis time takes part in the grouping. This is an O(n²) loop.
   2. **Number of groups.** It comes from convinfo `ncnumgrp`, taking the value of the *last* subtype read for that ob type. With `random_cvgrp`, it comes from `ngrps_<type>ob` instead.
   3. **Scaling.** Grid coordinates are scaled to (0, 1] by `nlon` and `nlat`.
   4. **Minimum count.** It runs only when `ncross > ngrpsfact0 * ngroups`. Otherwise it prints "not enough obs" and does nothing.
   5. **`hilbert(ngroups, n, x, y, test_set)`** (`hilbert_curve.f90`), Purser's code:
      - `xy_to_hil` (`phil.f90`) maps each (x, y) to a 15-digit base-4 Hilbert-curve code;
      - `hil_to_r` turns that code into a real number in [0, 1), the ob's position along the curve;
      - `bsort` sorts the obs along the curve into a linked list;
      - `getvalsets` (`phil1.f90`) walks the list and deals the obs into groups 1..N in turn. An ob joins a group only if its position on the curve is at least `xhskip = 3e-4` beyond that group's previous member. Otherwise it goes to the residual set, `test_set = 0`.

      Neighbours along a Hilbert curve are close in space. As a result:
      - each group is spread evenly and avoids close pairs;
      - in dense areas the extra obs fall into group 0, and those are always assimilated.
   6. **Withheld group.** The obs whose `test_set` equals the chosen group (`ncgroup` from convinfo, or `shuffle()` when `random_cvgrp`) get usage set to `ncmiter` (or `usagecv`). Their duplicates get `usage_dup`. The fractional part of usage is kept.
      - `shuffle()` (`support_2dvar.f90:3670`) picks group `mod(sum(iadate), N)`, so the withheld group rotates with the date.
      - In `read_prepbufr`, obs with usage ≥ 100 are marked as not used (`rusage = .false.`).
   7. **Output file.** All groups are written to `<type>obs_allcv_groups`, an unformatted file with this layout:
      - a header: `ncross, ngroups, chosen group`;
      - the counts per group;
      - for each group: its marker, then for each member the station id, provider and subprovider, then `kx, alon, alat, dlon, dlat`.
5. `destroy_hilbertcurve` frees the arrays.

**Restrictions** (from the module notes, `support_2dvar.f90:2885`):
- The number of groups must be the same for every subtype of an ob type.
- One processor must handle exactly one ob type: the module keeps global state.
- `xhskip`, `ngen` and the duplicate tolerance are hard-coded.

## Use cases in GSI

**Cross-validation (`hilbert_curve`).** It is only active for 2DVar (RTMA/URMA). The only regression test that turns it on is `RTMA` (`regression/regression_namelists.sh:268`).

| Reader | Ob types that call it |
|---|---|
| `read_prepbufr.f90` | Each ob-type branch sets `thisobtype_usage` for it (lines 2327–3093). The branches are t, uv, spd, ps, q, pw, sst, gust, vis, td2m, mxtm, mitm, pmsl, howv (3 branches), tcamt, lcbas and cldch: 17 ob types. `init`/`accum`/`apply`/`destroy` are at 922 / 2285 / 3136 / 3382. |
| `read_mitm_mxtm.f90` | mxtm and mitm read from their own files: `init`/`accum`/`apply`/`destroy` at 199 / 323 / 404 / 417. |
| `read_satmar.f90` | Sets `lhilbert` and the usage index (line 453), but never calls `accum` or `apply`, so it does nothing. |

`apply_hilbertcurve` also has entries for wspd10m, uwnd10m and vwnd10m, but no reader in this checkout calls it for those types.

So there are **2 active call sites** (`read_prepbufr`, `read_mitm_mxtm`), covering **up to 17 ob types** in RTMA.

**Density down-weighting (`denest`, phil2).** This is the separate feature:
- **When:** in `read_prepbufr.f90:3194–3360`, for `obstype == 'uv'` with kx 230, 231 or 233 (aircraft winds). Only obs at pressure ≥ 100 hPa are used. It runs when convinfo has `ithin_conv ≥ 5` and `rmesh`, `pmesh` > 0. It doesn't depend on `twodvar_regional` or `hilbert_curve`, so it can run in the global GSI.
- **What:** `denest` estimates the local ob density by averaging over 13 randomly oriented 3-D Hilbert curves (`nrand=13`, `qcmod.f90:470`) and returns a weight.
- **Where the weight goes:** it is stored in `cdata_all(26,:)`. `setupw.f90:1220` scales the ob error by it (`ratio_errors * sqrt(w)`), and it is saved in the diag file (`rdiagbuf(10)`).
- **Other readers:** `read_satwnd`, `read_gsb`, `read_saildrone` and `hdraobmod` just set the weight to 1.

## Standalone program: done, in `../hilbert_cv/`

```
hilbert_cv/
  Makefile                   make            (gfortran; make FC=ifx for Intel)
  src/phil.f90               copied unchanged from GSI
  src/phil1.f90              copied unchanged
  src/hilbert_curve.f90      copied unchanged
  src/kinds.f90              small stand-ins for GSI kinds/constants (r_kind = double)
  src/hilbert_cv_main.f90    driver: does what apply_hilbertcurve does
  test/make_synthetic.py     synthetic CONUS-like obs (40 clusters + uniform, 300 duplicates)
  test/check_groups.py       spacing check
```

**Usage:**

```
hilbert_cv <infile> <outfile> <ngroups> <nlon> <nlat> [ngrpsfact0]
```

- **Input:** one ob per line: `station_id kx lon lat x_grid y_grid time_offset_h`.
- **Output:** the same columns plus `group`: 1..N; 0 = in no group; −1 = duplicate removed.
- **What it does:** the driver repeats steps 4.1, 4.3, 4.4 and 4.5 above. It doesn't pick a withheld group, so the caller can use any group, or loop over all N for full N-fold cross-validation.
- **What it leaves out:** the convinfo and namelist logic, the usage rewriting, and the `allcv_groups` file.

**Test result** (20,000 obs plus 300 duplicates, 1799 × 1059 grid, 5 groups, 0.17 s):

| | Count / spacing |
|---|---|
| Duplicates removed | 300 |
| Groups 1–5 | 2,118–2,119 obs each |
| Group 0 | 9,408 obs |
| Median distance to the nearest ob in the same group | 18 grid units (random subset of the same size: 12) |
| 5th percentile of that distance | 7 grid units (random: 3) |
| Local density of group-0 obs | 27 obs within 50 km (grouped obs: 14) |

So the groups are spread out as intended, and the obs left out of every group sit mostly in the dense clusters.

## Next steps

1. **Check against GSI exactly.** Take one RTMA/URMA cycle's `tobs_allcv_groups` (or another type), convert it to the input format, and confirm the groups match. Needed for this:
   - a small reader for the unformatted file;
   - input in GSI's read order, because duplicate handling and `bsort` ties depend on the order.
2. **Expose hard-coded settings.** If needed, make `xhskip` (3e-4), `ngen` (15) and the duplicate tolerance options. `xhskip` sets how many obs end up in group 0.
3. **Possible additions:**
   - read lon/lat plus a map projection instead of grid x/y;
   - a Python wrapper (f2py) for use with `ai_rtma`;
   - replace the O(n²) duplicate check with a sort, for very large ob counts.
4. **If the density weighting is wanted too:** `denest` (`phil2.f90` plus `phil0.f90`, about 2,500 lines, only `kinds`) can be extracted the same way.
