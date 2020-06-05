# compiler 
#XLF= ifort -O3 -parallel #-g -traceback #-O3 -parallel -heap-arrays #-qopenmp #-O3 -parallel
#XLF= ifort -g -traceback 
#XLF= mpiifort -O3 -parallel#-g -debug full #-O3 -parallel
#XLF=/software/openmpi-1.8.2/intel/bin/mpifort -O3#/software/openmpi/3.1.2/intel/bin/mpifort
# options
#OPT= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= eqm3_admat
EXE_PATH=../../${EMPM_BIN_DIR}
#
all: types_eqm.o input_sp.o anglib.o base.o admat_cp_par.o eqm.o
	$(XLFMPI) types_eqm.o input_sp.o anglib.o base.o admat_cp_par.o eqm.o -o ${EXE_PATH}${EXF}
	
types_eqm.o: types_eqm.f90
	${XLFMPI} -c types_eqm.f90

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

base.o: base.f90
	${XLFMPI} -c base.f90

admat_cp_par.o: admat_cp_par.f90
	${XLFMPI} -c admat_cp_par.f90

eqm.o: eqm.f90
	${XLFMPI} -c eqm.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}




	
