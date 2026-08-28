#!/bin/bash
#SBATCH -J grib2p
#SBATCH -A da-cpu
#SBATCH -q debug 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -t 00:30:00
#SBATCH --mail-user=$LOGNAME@noaa.gov
#SBATCH --mem=0
#SBATCH -o jobs/grib2p%J.out
#SBATCH -e jobs/grib2p%J.err
set -x 
set -e


# load environment needed to run python scripts
source  /home/Xin.C.Jin/modules/my_eva.sh

python ocelot/aux/convert_urma_grib_to_parquet.py \
    --input_root /scratch3/NCEPDEV/gpu-emc-ai/Annette.Gibbs/aardvark_OK/OK_data/urma2p5_OK_grb2 \
    --output_dir /scratch3/NCEPDEV/stmp/Xin.C.Jin/data/ocelot/data_v6/urma_ok \
    --num_cores 8
