#export OMP_NUM_THREADS=24
#echo "# of MPI procs"
#export NUM_MPI_PROCS=24
#echo $NUM_MPI_PROCS
cd eqm2_MPI/
source ../mpi_openmp_set.sh
echo "Energy threshold"
cat en_trun.dat
./eqm_svd > log_eqm2 2>error_eqm2
