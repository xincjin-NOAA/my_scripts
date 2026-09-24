# read_rtma_diag

This code was refactored from Matthew Morris’s original work. Many thanks to Matthew for sharing the code and for his help and support during its development.


Reads GSI binary conventional-data diagnostic files (big-endian unformatted Fortran) and writes the observation metadata to a NetCDF4 file via the ncdiag library.

Supported observation types match those written by GSI setup routines:
`setupps`, `setupt`, `setupq`, `setuppw`, `setupuv`, `setupsst`, `setupgps`.

## Output

The program writes `diag_results.nc4` containing one record per observation with fields:

| Field | Description |
|---|---|
| `Variable` | Observation variable (t, q, uv, ps, …) |
| `Station_ID` | Station identifier |
| `Provider_Name` / `Subprovider_Name` | Data provider |
| `Observation_Type` | Prepbufr observation type |
| `Latitude` / `Longitude` | Degrees |
| `Pressure` | hPa |
| `Height` | Meters |
| `Time` | Hours relative to analysis time |
| `Analysis_Use_Flag` | 1 = used, -1 = monitored |
| `Observation` | Observed value |
| `Obs_Minus_Forecast_adjusted` | O−B after bias correction |
| `Observation_Error` | Final observation error |

## Dependencies

- [NCEPLIBS-ncdiag](https://github.com/NOAA-EMC/NCEPLIBS-ncdiag) — provides the `ncdiag::ncdiag` CMake target
- Fortran compiler: Intel (`ifort`/`ifx`) or GFortran

## Build

### Using `build.sh` (recommended)

```sh
CMAKE_PREFIX_PATH=/path/to/ncdiag/install bash build.sh
```

The script configures, builds, and installs to `./install` by default.
Available environment overrides:

| Variable | Default | Description |
|---|---|---|
| `CMAKE_PREFIX_PATH` | _(unset)_ | Path to ncdiag (and NetCDF) install |
| `BUILD_TYPE` | `Release` | CMake build type (`Debug`, `Release`, …) |
| `BUILD_DIR` | `./build` | Directory for CMake build files |
| `INSTALL_PREFIX` | `./install` | Installation destination |
| `BUILD_JOBS` | `8` | Parallel make jobs |
| `BUILD_CLEAN` | `YES` | Set to `NO` to reuse an existing build dir |

### Using CMake directly

```sh
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/ncdiag/install
cmake --build build
cmake --install build
```

The `-convert big_endian` flag (Intel) or `-fconvert=big-endian` (GFortran) is applied automatically to read the big-endian GSI diagnostic files.

## Usage

```sh
./read_diag_conv_rtma.exe diag_conv_ges.2024010100
```

The input file is a GSI diagnostic file. The output `diag_results.nc4` is written to the current directory.
