echo "Calculation of AD matrices"
#cat $PBS_NODEFILE > nodes.txt
mpirun -np $NUM_MPI_PROCS ./eqm3_admat_real > log_ad 2>error_ad
