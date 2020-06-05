# compiler 
#XLF= ifort -O3 -parallel -qopenmp #-g -traceback #-O3 -parallel -heap-arrays #-qopenmp #-O3 -parallel
#XLF= ifort -g -traceback 
# options
OPT= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= eqm3_svd
EXE_PATH=../../${EMPM_BIN_DIR}
#
all: types_eqm.o input_sp.o anglib.o base.o dens_list.o read_admat.o choles.o cm_ort_svd.o ham.o dmat_test.o eqm.o 
	$(XLF) types_eqm.o input_sp.o anglib.o base.o dens_list.o read_admat.o choles.o cm_ort_svd.o ham.o dmat_test.o eqm.o ${OPT} -o ${EXE_PATH}${EXF}

types_eqm.o: types_eqm.f90
	${XLF} -c types_eqm.f90

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

base.o: base.f90
	${XLF} -c base.f90

dens_list.o: dens_list.f90
	${XLF} -c dens_list.f90

read_admat.o: read_admat.f90
	${XLF} -c read_admat.f90	

choles.o: choles.f90
	${XLF} -c choles.f90

cm_ort_svd.o: cm_ort_svd.f90
	${XLF} -c cm_ort_svd.f90

#ham_geev_compl.o: ham_geev_compl.f90
#	${XLF} -c ham_geev_compl.f90
ham.o : ham.f90
	${XLF} -c ham.f90

dmat_test.o: dmat_test.f
	${XLF} -c dmat_test.f

eqm.o: eqm.f90
	${XLF} -c eqm.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXE_PATH}${EXF}




	
