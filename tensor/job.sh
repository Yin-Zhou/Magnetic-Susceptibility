#!/bin/bash
#
#SBATCH --verbose
#SBATCH --account=hc76
#SBATCH --job-name="MATRIX"  
#SBATCH --output="%j.%N.out"  
#SBATCH --partition=chen,parallel,share,aquila
#SBATCH --nodes=1 
#SBATCH -n 1
#SBATCH --export=ALL  
#SBATCH -t 10-00:00:00
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=yz3433@nyu.edu

#module load icc/18.0.5

module purge
module load intel2018/psx  mkl/19.0.0  gsl/2.5  fftw/intel/3.3.4  python/gnu/2.7.10

touch RUNNING

#mpirun -np 16 /gpfsnyu/scratch/jm7305/tmp/testforjm/z-vasp.5.4.1new/bin/vasp_std 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfsnyu/packages/gsl/2.5/lib

# export PYTHONPATH=$PYTHONPATH:/gpfsnyu/scratch/yz3433/LOOP-4.1.1

ulimit -s unlimited
ulimit -l unlimited
mpirun -np 1 /gpfsnyu/scratch/yz3433/eigenvalue/power > errors.out 2 >> errors.out


touch done

\rm RUNNING
