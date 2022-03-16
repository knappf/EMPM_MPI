cd eqm3_MPI/
source ../mpi_openmp_set.sh
echo "Energy threshold"
cat en_trun.dat
./eqm3_svd <ip_j_int_3ph  > log_eqm3 2>error_eqm3
