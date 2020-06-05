
all:
	@$(MAKE) -C HFB_code  
	@$(MAKE) -C fmat
	@$(MAKE) -C tda
	@$(MAKE) -C eqm2_MPI/phon_dens1
	@$(MAKE) -C eqm2_MPI/phon_int
	@$(MAKE) -C eqm2_MPI/eqm2_svd
	@$(MAKE) -C eqm2_MPI/eqm2_admat
	@$(MAKE) -C eqm3_MPI/eqm3_svd
	@$(MAKE) -C eqm3_MPI/eqm3_admat
	@$(MAKE) -C eqm3_MPI/phon_dens2_MPI	
	@$(MAKE) -C nondiag
	@$(MAKE) -C fullham

HFB:
	@$(MAKE) clean -C HFB_code
	@$(MAKE) -C HFB_code

fmat:
	@$(MAKE) clean -C fmat
	@$(MAKE) -C fmat

tda:
	@$(MAKE) clean -C tda
	@$(MAKE) -C tda

phon_dens1:
	@$(MAKE) clean -C eqm2_MPI/phon_dens1
	@$(MAKE) -C eqm2_MPI/phon_dens1

phon_int:
	@$(MAKE) clean -C eqm2_MPI/phon_int
	@$(MAKE) -C eqm2_MPI/phon_int

eqm2_svd:
	@$(MAKE) clean -C eqm2_MPI/eqm2_svd
	@$(MAKE) -C eqm2_MPI/eqm2_svd

eqm2_admat:
	@$(MAKE) clean -C eqm2_MPI/eqm2_admat
	@$(MAKE) -C eqm2_MPI/eqm2_admat

eqm3_svd:
	@$(MAKE) clean -C eqm3_MPI/eqm3_svd
	@$(MAKE) -C eqm3_MPI/eqm3_svd

eqm3_admat:
	@$(MAKE) clean -C eqm3_MPI/eqm3_admat
	@$(MAKE) -C eqm3_MPI/eqm3_admat

phon_dens2:
	@$(MAKE) clean -C eqm3_MPI/phon_dens2_MPI
	@$(MAKE) -C eqm3_MPI/phon_dens2_MPI

clean:
	@$(MAKE) clean -C HFB_code
	@$(MAKE) clean -C fmat
	@$(MAKE) clean -C tda
	@$(MAKE) clean -C eqm2_MPI/phon_dens1
	@$(MAKE) clean -C eqm2_MPI/phon_int
	@$(MAKE) clean -C eqm2_MPI/eqm2_svd
	@$(MAKE) clean -C eqm2_MPI/eqm2_admat
	@$(MAKE) clean -C eqm3_MPI/eqm3_svd
	@$(MAKE) clean -C eqm3_MPI/eqm3_admat
	@$(MAKE) clean -C eqm3_MPI/phon_dens2_MPI	
	@$(MAKE) clean -C nondiag
	@$(MAKE) clean -C fullham		


