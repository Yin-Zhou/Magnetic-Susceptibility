#!/bin/bash
#
#SBATCH --verbose
#SBATCH --account=hc76
#SBATCH --job-name="do_EP"  
#SBATCH --output="%j.%N.out"  
#SBATCH --partition=chen,parallel,share,aquila
#SBATCH --nodes=1 
#SBATCH -n 40
#SBATCH --export=ALL  
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=yz3433@nyu.edu

#module load icc/18.0.5

module purge
module load intel2018/psx  mkl/19.0.0  gsl/2.5  fftw/intel/3.3.4  python/gnu/2.7.10

touch RUNNING

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfsnyu/packages/gsl/2.5/lib

export PYTHONPATH=$PYTHONPATH:/gpfsnyu/scratch/yz3433/do_EP

ulimit -s unlimited
ulimit -l unlimited
mpirun -np 40 /gpfsnyu/home/yz3433/bin/do_EP_spectral > errors_spectral.out 2 >> errors_spectral.out
mpirun -np 40 /gpfsnyu/home/yz3433/bin/do_EP_self > errors_spectral.out 2 >> errors_spectral.out


touch done

\rm RUNNING
