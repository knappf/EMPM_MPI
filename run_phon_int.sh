
#  ******** IMPORTANT **************************************************************************
#  MPI version crashes for specific number of MPI porcesses, try to use smaller number if possible 
cd phon_int/
source ../mpi_openmp_set.sh
#$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./phon_int_MPI > log_phon_int 2>error
$MPIRUN_PATH/mpirun -np 4 ./phon_int_MPI > log_phon_int 2>error

