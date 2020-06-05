# compiler 
#XLF= ifort -O3 -parallel #-O3 -parallel
#XLF=/software/openmpi-1.8.2/intel/bin/mpifort -O3
#XLF=/software/openmpi/3.1.2/gcc/bin/mpifort 
#XLF= mpiifort -O3 #-g -debug full
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= phon_dens2_MPI
EXE_PATH=../../${EMPM_BIN_DIR}
#
all: anglib.o ro7_MPI.o input_sp.o phon_dens_MPI.o
	${XLFMPI} anglib.o ro7_MPI.o input_sp.o phon_dens_MPI.o ${OPT} -o ${EXE_PATH}${EXF}

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

ro7_MPI.o: ro7_MPI.f90
	${XLFMPI} -c ro7_MPI.f90

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

phon_dens_MPI.o: phon_dens_MPI.f90
	${XLFMPI} -c phon_dens_MPI.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}     

