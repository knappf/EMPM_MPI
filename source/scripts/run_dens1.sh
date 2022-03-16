#cat $PBS_NODEFILE > nodes.txt
cd eqm2_MPI/
source ../mpi_openmp_set.sh
/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error_phon_dens1

#/opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error_phon_dens1
