# compiler 
#
#XLF= ifort -O3 -parallel 
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= Tda
#EXE_PATH=../../run/tda/
EXE_PATH=../${EMPM_BIN_DIR}
#
all: input_sp.o cm.o read_inter.o tda_pn_coupl_ortog.o
	$(XLF) input_sp.o cm.o read_inter.o tda_pn_coupl_ortog.o ${OPT} -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

cm.o: cm.f90
	${XLF} -c cm.f90

read_inter.o: read_inter.f90
	${XLF} -c read_inter.f90

tda_pn_coupl_ortog.o: tda_pn_coupl_ortog.f90
	${XLF} -c tda_pn_coupl_ortog.f90


clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}     




	
