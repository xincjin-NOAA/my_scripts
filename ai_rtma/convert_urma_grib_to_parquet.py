"""
Converts URMA2p5 hourly analysis/background GRIB2 files (e.g.
urma2p5.t00z.2dvaranl_OK.grb2, urma2p5.t00z.2dvarges_OK.grb2, one pair per
hour/cycle) into ocelot3's Hive-partitioned Parquet layout, using the exact
file_base names obs_config_urma.py already expects for the "state" obs_type:

    <output_dir>/anal_urma_<year>.parquet/date=<YYYY-MM-DD>/cycle=<HH>/part-<source>.parquet
    <output_dir>/ges_urma_<year>.parquet/date=<YYYY-MM-DD>/cycle=<HH>/part-<source>.parquet

so the result is DIRECTLY readable by DARegionalGraphDataset for 4 of 5
features with NO obs_config changes -- see the specific-humidity note below
for the 5th.

Each row = one URMA2p5 grid cell (there is no dropna/subsetting across the 5
feature columns -- NaNs pass through and get masked downstream per-feature by
extract_features_from_df's own valid_mask logic, same as every other
instrument):
    latitude, longitude (degrees), obs_time (float, UNIX epoch seconds, from
    the GRIB's own valid_time -- not inferred from filename/directory, since
    tNNz only encodes the hour, not the date), plus 5 feature columns.

Confirmed via --inspect against a real file on 2026-08-26: this product has
sp/2t/10u/10v/2sh -- NOT 2d (dewpoint). obs_config_urma.py's "ges" feature
list assumes dpt_2maboveground (dewpoint, Kelvin-scale normalization stats);
this file only has 2 metre SPECIFIC HUMIDITY (kg/kg -- a different quantity
and a wildly different numeric scale). Per explicit decision: this converter
writes the raw 2sh value under a NEW column, "spfh_2maboveground", rather than
mislabeling it as dpt_2maboveground or approximating a dewpoint conversion.
*** obs_config_urma.py's "ges" features list and FEATURE_STATS still say
dpt_2maboveground, not spfh_2maboveground -- until that config is updated
(deliberately left as a separate follow-up, not done here), this 5th column
exists in the Parquet but ocelot3 won't read it: extract_features_from_df
only pulls the columns named in features_config["features"]. ***

*** UNVERIFIED ASSUMPTION -- flag before trusting a real training run ***
Grid cells are flattened row-major (numpy default .ravel() order) from
cfgrib's 2D latitude/longitude arrays. DARegionalGraphDataset.setup() asserts
this row order must match urma2p5_terrain.npy / urma2p5_slmask_nolakes.npz
1:1 (STATIC_CONTEXT_DIM context is added elementwise by row position, not by
coordinate lookup) -- if those static arrays were flattened in a different
order (e.g. Fortran/column-major, or from a different grid definition
entirely), this converter's grid points will silently misalign with the
static terrain/land-sea-mask context. Verify row count matches (this script
prints it) and, ideally, spot-check a known station's terrain height at the
corresponding flat index before trusting a real training run.

Usage:
    # sanity-check the shortName mapping on one real file first:
    python convert_urma_grib_to_parquet.py --inspect /path/to/urma2p5.t00z.2dvaranl_OK.grb2

    # then convert (recursively finds anl/ges files under --input_root):
    python convert_urma_grib_to_parquet.py \
        --input_root /path/to/urma_grib_root \
        --output_dir /scratch3/NCEPDEV/da/Xin.C.Jin/my_data/ocelot/data_v.../urma_state
    # -> .../ges_urma_2022.parquet/date=2022-01-01/cycle=00/part-urma2p5.t00z.2dvarges_OK.parquet, etc.

    # try on a handful of files before a full run:
    python convert_urma_grib_to_parquet.py --input_root ... --output_dir ... --limit 4
"""
import argparse
import multiprocessing as mp
import os
from glob import glob

import numpy as np
import pandas as pd
import xarray as xr

# ocelot3 feature name -> GRIB2 shortName (as decoded by ecCodes/cfgrib).
# Confirmed 2026-08-26 via --inspect against a real URMA2p5 OK file: sp/2t/10u/
# 10v/2sh (NOT 2d/dewpoint -- see module docstring's specific-humidity note).
# "spfh_2maboveground" is not yet in obs_config_urma.py's "ges" feature list.
SHORTNAME_MAP = {
    "tmp_2maboveground": "2t",
    "pres_surface": "sp",
    "vgrd_10maboveground": "10v",
    "ugrd_10maboveground": "10u",
    "spfh_2maboveground": "2sh",
}

EPOCH = np.datetime64("1970-01-01T00:00:00")


def inspect_grib(path: str) -> None:
    """Print every message's shortName/name/typeOfLevel/level in a GRIB2 file,
    so SHORTNAME_MAP can be checked/corrected against the real product before
    trusting a converted batch."""
    import eccodes

    print(f"Messages in {path}:")
    with open(path, "rb") as f:
        i = 0
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            i += 1
            try:
                short_name = eccodes.codes_get(gid, "shortName")
                name = eccodes.codes_get(gid, "name")
                type_of_level = eccodes.codes_get(gid, "typeOfLevel")
                level = eccodes.codes_get(gid, "level")
                print(f"  [{i}] shortName={short_name!r:10} name={name!r:35} "
                      f"typeOfLevel={type_of_level!r:12} level={level}")
            finally:
                eccodes.codes_release(gid)
    print(f"\nSHORTNAME_MAP currently expects: {SHORTNAME_MAP}")


def read_field(path: str, short_name: str, spatial_dims: tuple = None) -> xr.DataArray:
    """Returns the field transposed onto a canonical 2D (dim0, dim1) axis order --
    either the grid's native dims (e.g. curvilinear 'y'/'x') or, when `spatial_dims`
    is passed (from a previously-read field in the same file), that exact order, so
    every field in a file ravels in the same order regardless of the order ecCodes
    reports its dims in for that particular message."""
    ds = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={"filter_by_keys": {"shortName": short_name}, "indexpath": ""},
    )
    data_vars = list(ds.data_vars)
    if len(data_vars) != 1:
        raise ValueError(
            f"Expected exactly 1 variable for shortName={short_name!r} in {path}, "
            f"got {data_vars}. Run --inspect on this file and check SHORTNAME_MAP."
        )
    da = ds[data_vars[0]].squeeze()
    if da.ndim != 2:
        raise ValueError(
            f"Field {short_name!r} in {path} has shape {da.shape}, dims {da.dims} after "
            f"squeeze -- expected a 2D spatial grid."
        )
    if spatial_dims is not None and set(da.dims) != set(spatial_dims):
        raise ValueError(
            f"{path}: field shortName={short_name!r} has dims {da.dims}, expected "
            f"{spatial_dims} (from the first field read in this file) -- fields appear to be "
            f"on different grids/dim-names within the same file."
        )
    return da.transpose(*(spatial_dims or da.dims))


def _grid_lat_lon(da: xr.DataArray, spatial_dims: tuple) -> tuple:
    """(lat, lon) as flat arrays matching da.values.ravel()'s order. Handles both a
    curvilinear/projected grid (2D latitude/longitude aux coords, e.g. native URMA2p5
    Lambert grid) and a regular lat-lon grid (1D latitude/longitude dimension coords,
    e.g. this OK-subset product, confirmed 2026-08-26: cfgrib gave 1D lat(171)/lon(331)
    for a 171x331=56601-cell field -- broadcast them to 2D instead of assuming curvilinear."""
    lat_coord, lon_coord = da["latitude"], da["longitude"]
    if lat_coord.ndim == 1 and lon_coord.ndim == 1:
        lat_da, lon_da = xr.broadcast(lat_coord, lon_coord)
    elif lat_coord.ndim == 2 and lon_coord.ndim == 2:
        lat_da, lon_da = lat_coord, lon_coord
    else:
        raise ValueError(
            f"Unexpected latitude/longitude coord ndim ({lat_coord.ndim}/{lon_coord.ndim}) -- "
            f"expected both 1D (regular grid) or both 2D (curvilinear grid)."
        )
    lat_da = lat_da.transpose(*spatial_dims)
    lon_da = lon_da.transpose(*spatial_dims)
    return lat_da.values.ravel(), lon_da.values.ravel()


def convert_file(path: str, output_dir: str, file_base: str) -> None:
    part_name = os.path.splitext(os.path.basename(path))[0]

    fields = {}
    lat = lon = valid_time = spatial_dims = None
    for feature_name, short_name in SHORTNAME_MAP.items():
        da = read_field(path, short_name, spatial_dims)
        if spatial_dims is None:
            spatial_dims = da.dims
        fields[feature_name] = da.values.ravel()
        if lat is None:
            lat, lon = _grid_lat_lon(da, spatial_dims)
            if lat.size != fields[feature_name].size or lon.size != fields[feature_name].size:
                raise ValueError(
                    f"{path}: latitude/longitude ({lat.size}/{lon.size} cells) don't match "
                    f"{feature_name!r}'s grid ({fields[feature_name].size} cells)."
                )
            vt = da["valid_time"].values if "valid_time" in da.coords else da["time"].values
            valid_time = pd.Timestamp(vt)
        elif da.values.size != lat.size:
            raise ValueError(
                f"{path}: field {feature_name!r} ({short_name!r}) has {da.values.size} "
                f"grid cells but {list(SHORTNAME_MAP)[0]!r} had {lat.size} -- fields are "
                f"on different grids within the same file, this converter assumes they match."
            )

    n = lat.size
    valid_pos = ~(np.isnan(lat) | np.isnan(lon))
    df = pd.DataFrame(
        {
            "latitude": lat[valid_pos].astype(np.float64),
            "longitude": lon[valid_pos].astype(np.float64),
            "obs_time": np.full(valid_pos.sum(), (valid_time.to_datetime64() - EPOCH) / np.timedelta64(1, "s")),
            **{k: v[valid_pos].astype(np.float64) for k, v in fields.items()},
        }
    )

    date_part = valid_time.strftime("%Y-%m-%d")
    cycle = valid_time.strftime("%H")
    year = date_part.split("-")[0]
    part_dir = os.path.join(output_dir, f"{file_base}_{year}.parquet", f"date={date_part}", f"cycle={cycle}")
    os.makedirs(part_dir, exist_ok=True)
    df.to_parquet(os.path.join(part_dir, f"part-{part_name}.parquet"), engine="pyarrow", index=False)

    print(f"  {os.path.basename(path)}: valid_time={valid_time} n_grid={n} n_written={len(df)} "
          f"-> {file_base}_{year}.parquet/date={date_part}/cycle={cycle}/")


def discover_files(root: str, glob_pattern: str) -> list:
    return sorted(glob(os.path.join(root, "**", glob_pattern), recursive=True))


def _convert_file_safe(args: tuple) -> bool:
    """multiprocessing.Pool-friendly wrapper: one bad/corrupt file shouldn't kill
    a ~35,000-file batch, so catch and report per-file instead of letting it raise."""
    path, output_dir, file_base = args
    try:
        convert_file(path, output_dir, file_base)
        return True
    except Exception as e:
        print(f"  [ERROR] {path}: {e}")
        return False


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inspect", metavar="FILE", default=None,
                    help="Print GRIB message shortNames/levels for one file and exit (no conversion).")
    p.add_argument("--input", nargs="+", default=[],
                    help="Explicit GRIB2 file paths and/or glob patterns. Combinable with --input_root.")
    p.add_argument("--input_root", default=None,
                    help="Directory to search recursively for --anl_glob/--ges_glob files.")
    p.add_argument("--anl_glob", default="*2dvaranl*.grb2", help="Glob (used with --input_root) for analysis files.")
    p.add_argument("--ges_glob", default="*2dvarges*.grb2", help="Glob (used with --input_root) for background files.")
    p.add_argument("--output_dir", default=None, help="Root directory to write anal_urma_<year>.parquet/ and ges_urma_<year>.parquet/ into.")
    p.add_argument("--anl_file_base", default="anal_urma", help="Output dataset name for analysis files (default matches obs_config_urma.py's anal_zarr_name).")
    p.add_argument("--ges_file_base", default="ges_urma", help="Output dataset name for background files (default matches obs_config_urma.py's zarr_name).")
    p.add_argument("--limit", type=int, default=None, help="Only convert the first N files found (for a quick trial run).")
    p.add_argument("--num_cores", type=int, default=8,
                    help="Parallel worker processes (files are independent -- embarrassingly parallel). "
                         "Matches converter_zarr_multi_parall_vs2.py's --num_cores convention. Default 8.")
    args = p.parse_args()

    if args.inspect:
        inspect_grib(args.inspect)
        return

    if not args.output_dir:
        p.error("--output_dir is required (unless using --inspect).")
    if not args.input and not args.input_root:
        p.error("Pass --input and/or --input_root.")

    anl_files, ges_files = [], []
    for pattern in args.input:
        matches = sorted(glob(pattern))
        matches = matches if matches else [pattern]
        for f in matches:
            (anl_files if "anl" in os.path.basename(f) else ges_files).append(f)
    if args.input_root:
        anl_files.extend(discover_files(args.input_root, args.anl_glob))
        ges_files.extend(discover_files(args.input_root, args.ges_glob))

    anl_files = sorted(set(anl_files))
    ges_files = sorted(set(ges_files))
    if args.limit:
        anl_files = anl_files[: args.limit]
        ges_files = ges_files[: args.limit]

    missing = [f for f in anl_files + ges_files if not os.path.isfile(f)]
    if missing:
        raise FileNotFoundError(f"Input file(s) not found: {missing}")
    if not anl_files and not ges_files:
        raise ValueError(f"No input files matched --input={args.input} / --input_root={args.input_root}")

    n_cores = min(mp.cpu_count(), args.num_cores)
    print(f"Found {len(anl_files)} analysis file(s), {len(ges_files)} background file(s) -> {args.output_dir}")
    print(f"Using {n_cores} parallel worker process(es)")
    os.makedirs(args.output_dir, exist_ok=True)

    with mp.Pool(processes=n_cores) as pool:
        print(f"\n{'='*60}\nAnalysis files -> {args.anl_file_base}\n{'='*60}")
        anl_results = pool.map(_convert_file_safe, [(f, args.output_dir, args.anl_file_base) for f in anl_files])

        print(f"\n{'='*60}\nBackground files -> {args.ges_file_base}\n{'='*60}")
        ges_results = pool.map(_convert_file_safe, [(f, args.output_dir, args.ges_file_base) for f in ges_files])

    n_anl_ok, n_ges_ok = sum(anl_results), sum(ges_results)
    print(f"\nDone. Analysis: {n_anl_ok}/{len(anl_files)} succeeded. "
          f"Background: {n_ges_ok}/{len(ges_files)} succeeded.")
    if n_anl_ok < len(anl_files) or n_ges_ok < len(ges_files):
        print("Some files failed -- see [ERROR] lines above.")


if __name__ == "__main__":
    main()
