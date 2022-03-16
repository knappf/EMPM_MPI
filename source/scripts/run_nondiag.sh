source ../mpi_openmp_set.sh
#export NUM_MPI_PROCS=30
#/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun -n $NUM_MPI_PROCS ./Nondiag12_MPI > log_nondiag12 2>error_nondiag12 
#/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun -n $NUM_MPI_PROCS ./Nondiag23_MPI > log_nondiag23 2>error_nondiag23
mpirun -n $NUM_MPI_PROCS ./Nondiag12_MPI > log_nondiag12 2>error_nondiag12
mpirun -n $NUM_MPI_PROCS ./Nondiag23_MPI > log_nondiag23 2>error_nondiag23


