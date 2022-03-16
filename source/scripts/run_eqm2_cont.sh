cd eqm2_cont/
source ../mpi_openmp_set.sh
echo "Energy threshold"
cat en_trun.dat
./eqm_svd > log_eqm2 2>error_eqm2
