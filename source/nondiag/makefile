# compiler 
#XLF= ifort -O3 -parallel
#XLF= ifort -g -traceback 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= Nondiag
EXE_PATH=../${EMPM_BIN_DIR}
#
all: anglib.o input_sp.o interacndg.o nondiag.o
	$(XLF) anglib.o input_sp.o interacndg.o nondiag.o ${OPT} -o ${EXE_PATH}${EXF}

anglib.o: anglib.f90
	${XLF} -c anglib.f90

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

#interacndg_new_fast_cr_new.o: interacndg_new_fast_cr_new.f
#	${XLF} -c interacndg_new_fast_cr_new.f

interacndg.o: interacndg.f90
	${XLF} -c interacndg.f90

nondiag.o: nondiag.f
	${XLF} -c nondiag.f

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}     

