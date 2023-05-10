#!/bin/bash
#SBATCH -A tra23_units
#SBATCH -p m100_usr_prod
#SBATCH --time 02:00:00       # format: HH:MM:SS
#SBATCH -N 8                  # nodes
#SBATCH --gres=gpu:4          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-node=32  # 8 tasks out of 128
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=wanda_parallel_jacobi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=walter.nadalin@studenti.units.it

module purge
module load hpc-sdk
module load spectrum_mpi

make clean
make mpi
make openacc

rm data/times.dat
echo -e "mode\tsize\titrs\tprcs\ttoto\t\tdata\t\tcomp\t\tcomm" >> data/times.dat

for dim in 15000 20000 25000 30000
do
        for iters in 1000
        do
                prc=32

                for value in {1..4}
                do
                        make mpirun prc=$prc pernode=32 persocket=16 dim=$dim iters=$iters
                        ((prc*=2))
                done

                prc=4

                for value in {1..4}
                do
                        make openaccrun prc=$prc pernode=4 persocket=2 dim=$dim iters=$iters
                        ((prc*=2))
                done

        done
done

make clean
