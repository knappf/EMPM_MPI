#cat $PBS_NODEFILE > nodes.txt
#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error
$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error

