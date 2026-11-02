#!/bin/bash
#SBATCH -J cdf2p
#SBATCH -A da-cpu
#SBATCH -q debug 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -t 00:30:00
#SBATCH --mail-user=$LOGNAME@noaa.gov
#SBATCH --mem=0
#SBATCH -o jobs/cdf2p%J.out
#SBATCH -e jobs/cdf2p%J.err
set -x 
set -e


# load environment needed to run python scripts
source  /home/Xin.C.Jin/modules/my_eva.sh

python convert_surface_netcdf_to_parquet.py \
    --input_root /scratch3/NCEPDEV/gpu-emc-ai/Annette.Gibbs/data_prep/aardvark_OK/OK_data_final \
    --file_glob 'obs_netcdf/*.nc' \
    --output_dir /scratch5/purged/Xin.C.Jin/data/v2/diag_urma_ok \
    --file_base_prefix surface_obs

