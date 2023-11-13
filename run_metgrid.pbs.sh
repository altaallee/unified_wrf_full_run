#!/bin/bash

#PBS -N metgrid 
#PBS -o pbs_out_metgrid
#PBS -S /bin/bash
#PBS -l select=1:ncpus=40:model=sky_ele
##PBS -q devel
#PBS -l walltime=0:15:00
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

cd ../WPS
#mpiexec -np 40 ./geogrid.exe > out_geogrid
mpiexec -np 40 ./metgrid.exe > out_metgrid
