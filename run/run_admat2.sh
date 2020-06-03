echo "Calculation of 2-phonon AD matrices"
#cat $PBS_NODEFILE > nodes.txt
#/opt/intel/compilers_and_libraries_2020/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./eqm_admat > log_ad 2>error_ad
#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./eqm_admat > log_ad 2>error_ad
#source ../mpi_openmp_set.sh
$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./eqm_admat > log_ad 2>error_ad


