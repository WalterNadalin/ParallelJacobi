#SBATCH -A tra23_units
#SBATCH -p m100_usr_prod
#SBATCH --time 02:00:00       # format: HH:MM:SS
#SBATCH -N 8                  # nodes
#SBATCH --ntasks-per-node=32 # tasks out of 128
#SBATCH --gres=gpu:4          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=wanda_parallel_multiplication
#SBATCH --mail-type=ALL
#SBATCH --mail-user=walter.nadalin@studenti.units.it

module load autoload hpc-sdk
module load autoload spectrum_mpi

make clean
make benchmark mode=mpi
make benchmark mode=openacc

rm data/times.dat
echo -e "sign\tsize\titrs\tprcs\tcomm\t\tcomp" >> data/times.dat

for iters in {5000..10000..5000}
do
        for dim in {5000..10000..5000}
        do
                prc=32

                for value in {1..4}
                do
                        mpirun -np $prc -npernode 32 --map-by socket --bind-to core ./mpi_jacobi.x $dim $iters
                        ((prc*=2))
                done

                prc=4

                for value in {1..4}
                do
                        mpirun -np $prc -npernode 4 -npersocket 2 --bind-to core ./cuda_jacobi.x $dim $iters
                        ((prc*=2))
                done

        done
done

make clean
             
