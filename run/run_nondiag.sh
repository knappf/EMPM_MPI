source ../mpi_openmp_set.sh
#export NUM_MPI_PROCS=30
#/afs/ics.muni.cz/software/openmpi-1.8.2/intel/bin/mpirun -n $NUM_MPI_PROCS ./Nondiag12_MPI > log_nondiag12 2>error_nondiag12 
#/afs/ics.muni.cz/software/openmpi-1.8.2/intel/bin/mpirun -n $NUM_MPI_PROCS ./Nondiag23_MPI > log_nondiag23 2>error_nondiag23

mpirun -np $NUM_MPI_PROCS ./Nondiag12_MPI > log_nondiag12 2>error_nondiag12
#/opt/intel/oneapi/mpi/2021.7.0/bin/mpirun -n $NUM_MPI_PROCS ./Nondiag12_MPI > log_nondiag12 2>error_nondiag12
#mpirun -n $NUM_MPI_PROCS ./Nondiag23_MPI > log_nondiag23 2>error_nondiag23


