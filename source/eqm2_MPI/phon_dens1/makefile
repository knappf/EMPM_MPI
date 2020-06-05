# compiler 
#XLFMPI= mpiifort -O3 #g  
#XLFMPI=/afs/ics.muni.cz/software/openmpi-1.8.2/intel/bin/mpifort -O3
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= phon_dens1_MPI
EXE_PATH=../../${EMPM_BIN_DIR}
#
all: eof_mod.o input_sp.o anglib.o dens1_MPI.o phon_dens1_part_cut_MPI.o
	${XLFMPI} eof_mod.o input_sp.o anglib.o dens1_MPI.o phon_dens1_part_cut_MPI.o ${OPT} -o ${EXE_PATH}${EXF}

eof_mod.o: eof_mod.f90
	${XLFMPI} -c eof_mod.f90

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

dens1_MPI.o: dens1_MPI.f90
	${XLFMPI} -c dens1_MPI.f90

phon_dens1_part_cut_MPI.o: phon_dens1_part_cut_MPI.f90
	${XLFMPI} -c phon_dens1_part_cut_MPI.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}     

