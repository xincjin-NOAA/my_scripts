"""
Converts monthly, per-variable-station-indexed surface-obs NetCDF files (e.g.
surface_oklahoma_processed.nc) into ocelot3's Hive-partitioned Parquet layout,
one dataset per physical variable:

    <output_dir>/<file_base_prefix>_<var>_<year>.parquet/date=<YYYY-MM-DD>/cycle=<HH>/part-0.parquet

matching the "new" two-level structure ParquetDataManager.get_data_for_bin
tries first (dataset_timeseries.py) -- so the result is directly readable by
ParquetDataManager once an obs_config entry points at file_base
"<file_base_prefix>_<var>".

Input format assumed (auto-detected per file, not hardcoded to t/q/u/v/ps):
for each physical variable `v`, a data_var `v(station_index_v, time)` plus two
per-station coordinate arrays `v_lon(station_index_v)` / `v_lat(station_index_v)`.
Every variable can have its own station count/index -- there is no assumption
that different variables share stations. Missing obs are NaN (`_FillValue`)
and are dropped, not written as rows.

Each row written = one (station, hour) observation:
    latitude, longitude (degrees), obs_time (float, UNIX epoch seconds -- the
    column name/units ParquetDataManager's extract_features_from_df expects),
    <var> (the raw physical value, native units from the NetCDF).

Data is laid out on disk as one directory per month, named YYYYMM (e.g.
202201, 202202, ..., 202312), each with a "surface_obs_netcdf" subdirectory
holding one or more region files (e.g. surface_oklahoma_processed.nc), i.e.
<input_root>/202201/surface_obs_netcdf/surface_oklahoma_processed.nc. Pass the
directory containing the YYYYMM subdirs via --input_root and every file
matching --file_glob (default 'surface_obs_netcdf/*.nc') under every YYYYMM
dir is found and converted automatically. --input (explicit files/globs) still
works too, and can be combined with --input_root.

Because multiple region files can fall in the same month -- and therefore the
same date=/cycle= partition -- each source file writes its own part file
(named after the source file, not a fixed "part-0.parquet") so files from
different regions accumulate in a partition instead of overwriting each other.
Re-running the same source file again overwrites only that file's own part,
so it's safe to re-run after fixing an error.

Usage:
    python convert_surface_netcdf_to_parquet.py \
        --input_root /path/to/monthly_data_root \
        --output_dir /scratch3/NCEPDEV/da/Xin.C.Jin/my_data/ocelot/data_v.../surface_netcdf \
        --file_base_prefix surface_obs
    # -> .../surface_obs_t_2022.parquet/date=2022-01-01/cycle=00/part-surface_oklahoma_processed.parquet, etc.
"""
import argparse
import os
import re
from glob import glob

import numpy as np
import pandas as pd
import xarray as xr

EPOCH = np.datetime64("1970-01-01T00:00:00")
MONTH_DIR_RE = re.compile(r"^\d{6}$")


def discover_variables(ds: xr.Dataset) -> list:
    """Find every physical variable `v` that has sibling `v_lon`/`v_lat` coordinate
    arrays and a `time` dimension -- this is what makes it a convertible variable,
    not a fixed list of names, so it adapts to whatever variables a given file has."""
    found = []
    for name, da in ds.data_vars.items():
        if name.endswith("_lon") or name.endswith("_lat"):
            continue
        if "time" not in da.dims:
            continue
        lon_name, lat_name = f"{name}_lon", f"{name}_lat"
        if lon_name in ds.data_vars and lat_name in ds.data_vars:
            found.append(name)
    return found


def extract_obs_rows(ds: xr.Dataset, var: str) -> pd.DataFrame:
    """One row per valid (station, hour) observation for `var`."""
    station_dim = [d for d in ds[var].dims if d != "time"][0]
    da = ds[var].transpose(station_dim, "time")
    values = da.values  # (n_station, n_time)

    lon = ds[f"{var}_lon"].values  # (n_station,)
    lat = ds[f"{var}_lat"].values

    time_vals = ds["time"].values  # datetime64[*], via xarray's CF time decoding
    if not np.issubdtype(time_vals.dtype, np.datetime64):
        raise ValueError(
            f"ds['time'] did not decode to datetime64 (dtype={time_vals.dtype}); "
            "check the NetCDF's time units/calendar attributes."
        )

    station_idx, time_idx = np.where(~np.isnan(values))
    if station_idx.size == 0:
        return pd.DataFrame(columns=["latitude", "longitude", "obs_time", var])

    obs_val = values[station_idx, time_idx]
    obs_lon = lon[station_idx]
    obs_lat = lat[station_idx]
    obs_time_dt = time_vals[time_idx]

    valid_pos = ~(np.isnan(obs_lon) | np.isnan(obs_lat))
    obs_val, obs_lon, obs_lat, obs_time_dt = (
        obs_val[valid_pos],
        obs_lon[valid_pos],
        obs_lat[valid_pos],
        obs_time_dt[valid_pos],
    )

    obs_time_sec = (obs_time_dt - EPOCH) / np.timedelta64(1, "s")

    return pd.DataFrame(
        {
            "latitude": obs_lat.astype(np.float64),
            "longitude": obs_lon.astype(np.float64),
            "obs_time": obs_time_sec.astype(np.float64),
            var: obs_val.astype(np.float64),
            "_date": pd.DatetimeIndex(obs_time_dt).strftime("%Y-%m-%d"),
            "_cycle": pd.DatetimeIndex(obs_time_dt).strftime("%H"),
        }
    )


def write_partitions(df: pd.DataFrame, output_dir: str, file_base: str, part_name: str) -> int:
    n_written = 0
    for (date_part, cycle), group in df.groupby(["_date", "_cycle"]):
        year = date_part.split("-")[0]
        part_dir = os.path.join(
            output_dir, f"{file_base}_{year}.parquet", f"date={date_part}", f"cycle={cycle}"
        )
        os.makedirs(part_dir, exist_ok=True)
        group.drop(columns=["_date", "_cycle"]).to_parquet(
            os.path.join(part_dir, f"part-{part_name}.parquet"), engine="pyarrow", index=False
        )
        n_written += 1
    return n_written


def convert_file(nc_path: str, output_dir: str, file_base_prefix: str) -> None:
    # Part filename is derived from the source file so that multiple region
    # files landing in the same date=/cycle= partition (common, since several
    # regions' monthly files can share the same month) accumulate as separate
    # files in that partition dir instead of overwriting each other -- Parquet
    # readers (including ParquetDataManager) read every file in the directory.
    part_name = os.path.splitext(os.path.basename(nc_path))[0]

    print(f"\n{'='*60}")
    print(f"Opening {nc_path}")
    ds = xr.open_dataset(nc_path)
    try:
        variables = discover_variables(ds)
        if not variables:
            print(f"  No convertible variables found in {nc_path} (need `<var>(..., time)` "
                  f"plus `<var>_lon`/`<var>_lat`) -- skipping.")
            return
        print(f"  Found variables: {variables}")

        for var in variables:
            df = extract_obs_rows(ds, var)
            if df.empty:
                print(f"  [{var}] no valid (non-NaN) observations -- skipping.")
                continue
            file_base = f"{file_base_prefix}_{var}"
            n_partitions = write_partitions(df, output_dir, file_base, part_name)
            print(f"  [{var}] {len(df)} obs across {n_partitions} date/cycle partitions "
                  f"-> {os.path.join(output_dir, file_base + '_<year>.parquet')}")
    finally:
        ds.close()
    print(f"{'='*60}")


def discover_monthly_files(root: str, file_glob: str) -> list:
    """Every file matching `file_glob` under every YYYYMM-named subdirectory of `root`
    -- `file_glob` defaults to 'surface_obs_netcdf/*.nc' (the one confirmed
    subdirectory files live in under each month dir), not a recursive search."""
    if not os.path.isdir(root):
        raise NotADirectoryError(f"--input_root {root!r} is not a directory.")
    month_dirs = sorted(
        d for d in os.listdir(root)
        if MONTH_DIR_RE.match(d) and os.path.isdir(os.path.join(root, d))
    )
    if not month_dirs:
        raise ValueError(f"No YYYYMM-named subdirectories found under {root!r}.")
    print(f"Found {len(month_dirs)} month directories under {root}: {month_dirs[0]}..{month_dirs[-1]}")

    files = []
    for md in month_dirs:
        matched = sorted(glob(os.path.join(root, md, file_glob)))
        if not matched:
            print(f"  [warn] no files matching {file_glob!r} under {md}/")
        files.extend(matched)
    return files


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--input", nargs="+", default=[],
        help="One or more NetCDF file paths and/or glob patterns (e.g. 'surface_*_processed.nc'). "
             "Can be combined with --input_root.",
    )
    p.add_argument(
        "--input_root", default=None,
        help="Directory containing YYYYMM-named month subdirectories (e.g. 202201, 202202, ...); "
             "every file matching --file_glob under every such subdirectory is converted.",
    )
    p.add_argument(
        "--file_glob", default="surface_obs_netcdf/*.nc",
        help="Glob (relative to each YYYYMM dir) used with --input_root to find NetCDF files "
             "(default: 'surface_obs_netcdf/*.nc' -- confirmed subdirectory files live in).",
    )
    p.add_argument("--output_dir", required=True, help="Root directory to write <file_base>_<year>.parquet/ datasets into.")
    p.add_argument(
        "--file_base_prefix", default="surface_obs",
        help="Output dataset name per variable is '<prefix>_<var>' (default: 'surface_obs', "
             "e.g. 'surface_obs_t_2022.parquet').",
    )
    args = p.parse_args()

    if not args.input and not args.input_root:
        p.error("Pass --input and/or --input_root.")

    nc_files = []
    for pattern in args.input:
        matches = sorted(glob(pattern))
        nc_files.extend(matches if matches else [pattern])
    if args.input_root:
        nc_files.extend(discover_monthly_files(args.input_root, args.file_glob))
    nc_files = sorted(set(nc_files))

    missing = [f for f in nc_files if not os.path.isfile(f)]
    if missing:
        raise FileNotFoundError(f"Input file(s) not found: {missing}")
    if not nc_files:
        raise ValueError(f"No input files matched --input={args.input} / --input_root={args.input_root}")

    print(f"Converting {len(nc_files)} file(s) -> {args.output_dir}")
    os.makedirs(args.output_dir, exist_ok=True)

    for nc_path in nc_files:
        convert_file(nc_path, args.output_dir, args.file_base_prefix)

    print("\nDone.")


if __name__ == "__main__":
    main()
