#cat $PBS_NODEFILE > nodes.txt
#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error
cd phon_dens1/
source ../mpi_openmp_set.sh
$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error
