# compiler 
#XLF= mpiifort -O3 -parallel
#XLF= ifort -g -traceback 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= Nondiag12_MPI
EXE_PATH=../../${EMPM_BIN_DIR}
#
all: anglib.o input_sp.o interacndg.o nondiag12_MPI.o
	$(XLFMPI) anglib.o input_sp.o interacndg.o nondiag12_MPI.o ${OPT} -o ${EXE_PATH}${EXF}

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

#interacndg_new_fast_cr_new.o: interacndg_new_fast_cr_new.f
#	${XLF} -c interacndg_new_fast_cr_new.f

interacndg.o: interacndg.f90
	${XLFMPI} -c interacndg.f90

nondiag12_MPI.o: nondiag12_MPI.f90
	${XLFMPI} -c nondiag12_MPI.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}     

