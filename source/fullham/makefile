# compiler 
#XLF= ifort -O3 -parallel 
# option
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= Fullham
EXE_PATH=../${EMPM_BIN_DIR}
#
all: fullmat.o fullham.o
	$(XLF) fullmat.o fullham.o ${OPT} -o ${EXE_PATH}${EXF}

fullmat.o: fullmat.f90
	${XLF} -c fullmat.f90

fullham.o: fullham.f90
	${XLF} -c fullham.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}     

