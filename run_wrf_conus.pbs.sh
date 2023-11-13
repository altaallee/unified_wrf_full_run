#!/bin/bash

#PBS -N wrf 
#PBS -o pbs_out_wrf_conus
#PBS -S /bin/bash
#PBS -l select=30:ncpus=28:model=bro
##PBS -q devel
#PBS -l walltime=6:00:00
#PBS -j oe
#PBS -W group_list=s2395
#PBS -m abe


module purge
module load comp-intel/2020.4.304
module load mpi-hpe/mpt
module load hdf5/1.8.18_mpt
module load hdf4/4.2.12
module load netcdf/4.4.1.1_mpt

ulimit -s unlimited
cd $PBS_O_WORKDIR

cd ../WRF/run_conus
ln -sf ../../WPS/met_em* .
mpiexec -np 840 ./real.exe > out_real
mpiexec -np 840 ./wrf.exe > out_wrf
