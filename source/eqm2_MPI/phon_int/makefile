# compiler 
#XLF= mpif90 -O3 #-parallel 
#XLFMPI= mpiifort -O3 #-parallel
#XLF=/afs/ics.muni.cz/software/openmpi-1.8.2/intel/bin/mpifort -O3
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= phon_int_MPI
EXE_PATH=../../${EMPM_BIN_DIR}
#
all: eof_mod.o input_sp.o anglib.o interac_part_fast5_cut.o phon_int.o
	$(XLFMPI) eof_mod.o input_sp.o anglib.o interac_part_fast5_cut.o phon_int.o ${OPT} -o ${EXE_PATH}${EXF}

eof_mod.o: eof_mod.f90
	${XLFMPI} -c eof_mod.f90

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

interac_part_fast5_cut.o: interac_part_fast5_cut.f90 
	${XLFMPI} -c interac_part_fast5_cut.f90

phon_int.o: phon_int.f90
	${XLFMPI} -c phon_int.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}     



