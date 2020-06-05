echo "Calculation of AD matrices"
#cat $PBS_NODEFILE > nodes.txt
/opt/intel/compilers_and_libraries_2020/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./eqm3_admat > log_ad 2>error
#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./eqm3_admat > log_ad 2>error_ad
