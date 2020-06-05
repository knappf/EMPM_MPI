# compiler 
#XLF= ifort -O3 -parallel 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= Fmat
EXE_PATH=../${EMPM_BIN_DIR}
#
all: anglib.o fmat_unfi.o 
	$(XLF) anglib.o fmat_unfi.o ${OPT} -o ${EXE_PATH}${EXF}

anglib.o: anglib.f90
	${XLF} -c anglib.f90

fmat_unfi.o: fmat_unfi.f90
	${XLF} -c fmat_unfi.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}




	
