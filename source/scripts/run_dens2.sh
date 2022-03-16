echo "Calculation of 2-phonon densities"
#cat $PBS_NODEFILE > nodes.txt
mpirun -np $NUM_MPI_PROCS ./phon_dens2_MPI > log_dens2 2>error_dens2     
