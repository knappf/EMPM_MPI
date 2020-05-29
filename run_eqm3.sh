#export OMP_NUM_THREADS=24
#echo "# of MPI procs"
#export NUM_MPI_PROCS=24
#echo $NUM_MPI_PROCS
cd eqm3_MPI/
source ../mpi_openmp_set.sh
echo "Energy threshold"
cat en_trun.dat
./eqm3_svd > log_eqm3 2>error_eqm3
