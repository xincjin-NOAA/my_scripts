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

python ocelot/aux/convert_surface_netcdf_to_parquet.py \
    --input_root /scratch3/NCEPDEV/gpu-emc-ai/Annette.Gibbs/aardvark_OK/OK_data \
    --output_dir /scratch3/NCEPDEV/stmp/Xin.C.Jin/data/ocelot/data_v6/urma_ok \
    --file_base_prefix surface_obs

