#  ******** IMPORTANT **************************************************************************  
#  MPI version crashes for specific number of MPI porcesses, try to use smaller number if possible 
cd eqm2_MPI/
source ../mpi_openmp_set.sh
#$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./phon_int_MPI > log_phon_int 2>error_phon_int    
/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./phon_int_MPI > log_phon_int 2>error_phon_int
#/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./phon_int_MPI > log_phon_int 2>error_phon_int
