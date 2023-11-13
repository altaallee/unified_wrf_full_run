#!/bin/bash

##PBS -N wrf_postprocess_soundings
##PBS -o pbs_out_wrf_postprocess_soundings
##PBS -S /bin/bash
##PBS -l select=1:ncpus=24:model=has
##PBS -q devel
##PBS -l walltime=1:00:00
##PBS -j oe
##PBS -W group_list=s2395
##PBS -m abe


ulimit -c unlimited
cd $PBS_O_WORKDIR

module use -a /swbuild/analytix/tools/modulefiles
module load miniconda3/v4
export CONDA_ENVS_PATH=/home1/alee31/mambaforge/envs
source activate wrf_op

python3 -W ignore plot_soundings.py --date $DATE$RUN --hours $HOURS --ens $ENS
