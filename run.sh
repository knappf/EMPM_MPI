
export MPIRUN_PATH="/software/openmpi-1.8.2/intel/bin/"

export NUM_MPI_PROCS=24
export OMP_NUM_THREADS=24

# HFB code 
echo "HF calculation"
cd input 
#./HFB_DD > log_HFB
#
echo "F matrix calculation"
cd ../fmat
./fmat < input_space > log_fmat
echo "TDA calculation"
cd ../tda
./tda > log_TDA
#
echo "1-phonon densities calculation"
cd ../phon_dens1
./run_dens1.sh


#

