#!/bin/bash

set -eux

#GSI_ROOT=/lfs/h2/emc/da/save/xin.c.jin/git/gsi/osw_gsi
GSI_ROOT=/scratch3/NCEPDEV/da/Xin.C.Jin/git/osw_gsi
# Detect machine (sets MACHINE_ID)
source $GSI_ROOT/ush/detect_machine.sh

COMPILER=${COMPILER:-"intel"}

# Load modules
set +x
source $GSI_ROOT/ush/module-setup.sh
module use $GSI_ROOT/modulefiles
module load "gsi_${MACHINE_ID}.${COMPILER}"
module list
set -x

# Directory containing this script (and CMakeLists.txt)
readonly DIR_ROOT=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}")")" && pwd -P)

# User options
BUILD_TYPE=${BUILD_TYPE:-"Release"}
CMAKE_OPTS=${CMAKE_OPTS:-}
BUILD_DIR=${BUILD_DIR:-"${DIR_ROOT}/build"}
INSTALL_PREFIX=${INSTALL_PREFIX:-"${DIR_ROOT}/install"}

# Path to the dependencies; ensure module-defined roots are included
CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-""}
for root_var in "netCDF_ROOT" "NetCDF_ROOT" "netcdf_ROOT" "ncdiag_ROOT" "ncdiag_DIR"; do
    if [[ -n "${!root_var:-}" ]]; then
        CMAKE_PREFIX_PATH="${!root_var}:${CMAKE_PREFIX_PATH}"
    fi
done

#==============================================================================#

CMAKE_OPTS+=" -DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
CMAKE_OPTS+=" -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}"
[[ -n "${CMAKE_PREFIX_PATH}" ]] && CMAKE_OPTS+=" -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"

# Re-use or create a new BUILD_DIR (default: create new)
[[ ${BUILD_CLEAN:-"YES"} =~ [yYtT] ]] && rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"

cmake ${CMAKE_OPTS} "${DIR_ROOT}"
make -j "${BUILD_JOBS:-8}"
make install

exit
