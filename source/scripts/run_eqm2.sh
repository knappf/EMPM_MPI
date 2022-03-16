cd eqm2_MPI/
source ../mpi_openmp_set.sh
echo "Energy threshold"
cat en_trun.dat
./eqm_svd < ip_j_int_2ph > log_eqm2 2>error_eqm2
