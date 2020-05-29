echo "Calculation of 2-phonon densities"
#cat $PBS_NODEFILE > nodes.txt
/opt/intel/compilers_and_libraries_2020/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens2_MPI > log_dens2 2>error_dens2
#/software/openmpi/3.1.2/intel/bin/mpirun -np $NUM_MPI_PROCS phon_dens2_MPI > log_dens2 2>error_dens2


