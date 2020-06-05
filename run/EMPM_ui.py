#!/usr/bin/env python2
# 
# (kshell_home_dir)/bin/kshell_ui.py 
#    user interface to generate shell script.
# 

import sys, os, os.path, shutil
import readline


bindir = os.path.dirname( __file__ )
print 'directory where is the source code',bindir


print "\n" \
    + "----------------------------- \n" \
    + "  EMPM user interface \n" \
    + "     to generate job inputs and folders. \n" \
    + "-----------------------------\n "


os.system('mkdir inputs')   
outputfileMPI= open("inputs/mpi_openmp_set.sh", "w")

n_nodes = input('Set # of MPI nodes ')
n_th = input('Set # of OMP threads ')


lineMPI='export NUM_MPI_PROCS='+ str(n_nodes)
lineOMP='export OMP_NUM_THREADS='+str(n_th)
outputfileMPI.write('export MPIRUN_PATH="/software/openmpi-1.8.2/intel/bin/"      # Metacentrum urga, ursa, uruk\n')
outputfileMPI.write('#export MPIRUN_PATH="/opt/intel/compilers_and_libraries_2020.1.217/linux/mpi/intel64/bin/"    # ipnp16\n')
outputfileMPI.write('# number of MPI procs\n')
outputfileMPI.write(lineMPI)
outputfileMPI.write('\n')
outputfileMPI.write('echo "# MPI procs"\n')
outputfileMPI.write('echo $NUM_MPI_PROCS \n')
outputfileMPI.write('# number of OMP threads \n')
outputfileMPI.write(lineOMP)
outputfileMPI.write('\n')
outputfileMPI.write('#echo "# OMP threads"\n')
outputfileMPI.write('echo $OMP_NUM_THREADS\n')

outputfileMPI.close()

print ('***********')
print ('***********')
#sys.exit()

def print_var_dict( var_dict, skip=() ):
    ret = ""
    keys = var_dict.keys()
    keys.sort()
    for key in keys:
        if key in skip: continue
        v = var_dict[key]
        if isinstance(v, list): 
            vv = ""
            for i in v: vv += str(i)+", "
        elif isinstance(v, int) or isinstance(v, float):
            vv = str(v)
        else: vv = v
        ret += "  " + key + " = " + vv + "\n"
    return ret

class SimpleCompleter(object):
    def __init__(self, options):
        self.options = sorted(options)
        return

    def complete(self, text, state):
        response = None
        if state == 0:
            if text:
                self.matches = [s for s in self.options
                                if s and s.startswith(text)]
            else:
                self.matches = self.options[:]
        try:
            response = self.matches[state]
        except IndexError:
            response = None
        return response
        
output_ans = ""        
def raw_input_save(c=None):
    if c is None: r = raw_input()
    else: r = raw_input(c)
    global output_ans
    output_ans += r + '\n'
    return r        


#outputfileMPI= open(bindir+"../inputs/input.dat", "w")
outputfileHF= open("inputs/input.dat", "w")
print 'We count the oscillator shells as N = 0,1,2,3,4,... So the total number of shells is equal to (noscmax + 1)'

HF_inputs={
"A":16,
"Z":8,
"noscmax":4,
"noscmax12":8,
"noscmax123":12,
"if_2B":0,
"if_QTDA":0,
"tol":'1.0d-7',
"qpp":1.0,
"qnn":1.0,
"qpn":1.0,
"min_p":1,
"min_n":1,
"min_Y":1,
"max_p":10,
"max_n":10,
"max_Y":10,
"hw":16.0,
"s3":0.0,
"s3Y":1.0,
"if_ort":1,
"if_self":0,
"n_lam_occup":1,
"if_NAT":1,
}

while True: 
    print "\n --- input HF parameters --- "
    print print_var_dict(HF_inputs,skip=())
    
    ask = "modify parameter? \n" \
        + " (e.g.  A = 40 for parameter change\n" \
        + "        <CR>          for no more modification ) :\n"

    list_param = [ k +" = " for k in HF_inputs.keys() ]
    #print list_param
    readline.set_completer( SimpleCompleter(list_param).complete )
    
    if 'libedit' in readline.__doc__: # for Mac
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
            
    ans = raw_input_save(ask)
    readline.parse_and_bind("tab: None")
        
    ans = ans.strip()
    if not ans:
        break
    elif '=' in ans:
        arr = ans.split('=')
        arr = [ a.strip() for a in arr ] 
        if len(arr[1]) != 0:
            HF_inputs[ arr[0] ] = arr[1]
        else: 
            del HF_inputs[ arr[0] ]
    else: 
            print "ILLEGAL INPUT"
            
                    

#print  HF_inputs["A"]         
line1=str(HF_inputs["A"]) + ',' + str(HF_inputs["Z"] )+ '   ! mass and proton number\n'
line2=str(HF_inputs["noscmax"]) + ',' +str(HF_inputs["noscmax12"]) +','+ str(HF_inputs["noscmax123"]) + '   ! noscmax (maximal shell number) for HFB, noscmax12, noscmax123\n'
line3=str(HF_inputs["if_2B"]) + '   ! if_2B (= 0 noscmax12 is not applied on V^{2B}, = 1 noscmax12 is applied also on V^{2B})\n'
line4=str(HF_inputs["min_p"])+','+ str(HF_inputs["max_p"])+ ',' + str(HF_inputs["min_n"])+ ',' + str(HF_inputs["max_n"]) +','+str(HF_inputs["min_Y"])+ ',' + str(HF_inputs["max_Y"]) + '   ! min_p,max_p,min_n,max_n,min_Y,max_Y - config. space for TDA\n'
line5=str(HF_inputs["if_QTDA"])+'   ! if_QTDA (= 0 for TDA, = 1 for QTDA)\n'
line6=str(HF_inputs["tol"])+'   ! precision parameter\n'
line7=str(HF_inputs["hw"])+'   ! hbar*omega [MeV]\n'
line8=str(HF_inputs["qpp"])+','+str(HF_inputs["qnn"])+','+str(HF_inputs["qpn"]) + '   ! quenching of proton-proton, neutron-neutron, proton-neutron and F interaction\n'
line9=str(HF_inputs["s3"])+ '   ! strength of the NNN interaction\n'
line10=str(HF_inputs["s3Y"])+'    ! strength of the LambdaNN interaction\n'
line11=str(HF_inputs["if_ort"])+'     ! if_ort (= 0 without orthog., = 1 with orthog)\n'
line12=str(HF_inputs["if_self"])+'     ! if_self (= 0 Lambda non treated self-consistently, = 1 Lambda treated self-consistently)\n'
line13=str(HF_inputs["n_lam_occup"])+'     ! n_lam_occup  - which s.p. level is ocuppied by Lambda\n'
line14=str(HF_inputs["if_NAT"])+'    ! if_NAT (=1 definition of NO due to Tichai, =2 definition due to normalized MBPT funciton,##=3 definition due to EMPM correlated wf.)\n'



outputfileHF.write(line1)
outputfileHF.write(line2)
outputfileHF.write(line3)
outputfileHF.write(line4)
outputfileHF.write(line5)
outputfileHF.write(line6)
outputfileHF.write(line7)
outputfileHF.write(line8)
outputfileHF.write(line9)
outputfileHF.write(line10)
outputfileHF.write(line11)
outputfileHF.write(line12)
outputfileHF.write(line13)
outputfileHF.write(line14)
outputfileHF.close()        

#******************
print ('***********')
print ('***********')

#outputfileF= open(bindir+"../inputs/input_space", "w")
outputfileF= open("inputs/input_space", "w")
Fmat_inputs={
"min_p":1,
"min_n":1,
"max_p":10,
"max_n":10,
}

Fmat_inputs["min_p"]=HF_inputs["min_p"]
Fmat_inputs["min_n"]=HF_inputs["min_n"]
Fmat_inputs["max_p"]=HF_inputs["max_p"]
Fmat_inputs["max_n"]=HF_inputs["max_n"]

while True: 
    print "\n --- input FMAT parameters --- "
    print print_var_dict(Fmat_inputs,skip=())

    ask = "modify parameter? \n" \
        + " (e.g.  min_n = 2 for parameter change\n" \
        + "        <CR>          for no more modification ) :\n"

    list_param = [ k +" = " for k in Fmat_inputs.keys() ]
#print list_param
    readline.set_completer( SimpleCompleter(list_param).complete )

    if 'libedit' in readline.__doc__: # for Mac
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
        
    ans = raw_input_save(ask)
    readline.parse_and_bind("tab: None")
    
    ans = ans.strip()
    if not ans:
        break
    elif '=' in ans:
        arr = ans.split('=')
        arr = [ a.strip() for a in arr ] 
        if len(arr[1]) != 0:
            Fmat_inputs[ arr[0] ] = arr[1]
        else: 
            del Fmat_inputs[ arr[0] ]
    else: 
        print "ILLEGAL INPUT"

line1=str(Fmat_inputs["min_p"]) + ',' + str(Fmat_inputs["max_p"] )+ ' ! proton space\n'
line2=str(Fmat_inputs["min_n"]) + ',' + str(Fmat_inputs["max_n"] )+ ' ! neutron space\n'


outputfileF.write(line1)
outputfileF.write(line2)
outputfileF.close()

print ('***********')
print ('***********')
#******
#******
#outputfileF= open(bindir+"../inputs/input_space", "w")
outputfileTDA= open("inputs/input_tda_coup.dat", "w")
TDA_inputs={
#"A":16,
#"Z":8,
"minh_p":1,
"minh_n":1,
"maxh_p":3,
"maxh_n":3,
"maxp_p":10,
"maxp_n":10,
"minp_p":4,
"minp_n":4,
}


TDA_inputs["minh_p"]=Fmat_inputs["min_p"]
TDA_inputs["minh_n"]=Fmat_inputs["min_n"]
TDA_inputs["maxp_p"]=Fmat_inputs["max_p"]
TDA_inputs["maxp_n"]=Fmat_inputs["max_n"]


while True: 
    print "\n --- input TDA parameters --- "
    print print_var_dict(TDA_inputs,skip=())

    ask = "modify parameter? \n" \
        + " (e.g.  minh_n = 2 for parameter change\n" \
        + "        <CR>          for no more modification ) :\n"

    list_param = [ k +" = " for k in TDA_inputs.keys() ]
#print list_param
    readline.set_completer( SimpleCompleter(list_param).complete )

    if 'libedit' in readline.__doc__: # for Mac
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
        
    ans = raw_input_save(ask)
    readline.parse_and_bind("tab: None")
    
    ans = ans.strip()
    if not ans:
        break
    elif '=' in ans:
        arr = ans.split('=')
        arr = [ a.strip() for a in arr ] 
        if len(arr[1]) != 0:
            TDA_inputs[ arr[0] ] = arr[1]
        else: 
            del TDA_inputs[ arr[0] ]
    else: 
        print "ILLEGAL INPUT"

cm=raw_input("Apply Gramm-Scmidt CM orthogonalzation?: y/n ") 

line1='     '+str(HF_inputs["A"]) + '     ' + str(HF_inputs["Z"] )+ ' ! A,Z\n'
line2='     '+str(TDA_inputs["minh_n"]) + '      ' + str(TDA_inputs["maxh_n"] )+ '  !hole space for neutrons\n'
line3='     '+str(TDA_inputs["minh_p"]) + '      ' + str(TDA_inputs["maxh_p"] )+ '  !hole space for protns\n'
line4='     '+str(TDA_inputs["minp_n"]) + '      ' + str(TDA_inputs["maxp_n"] )+ '   !particle space for neutrons\n'
line5='     '+str(TDA_inputs["minp_p"]) + '      ' + str(TDA_inputs["maxp_p"] )+ '   !particle space for protons\n'
line6=cm + '    !Apply Gramm-Scmidt CM orthogonalzation? '

outputfileTDA.write(line1)
outputfileTDA.write(line2)
outputfileTDA.write(line3)
outputfileTDA.write(line4)
outputfileTDA.write(line5)
outputfileTDA.write(line6)

outputfileTDA.close()


print ('***********')
print ('***********')

eth2=input("2 phonon energy threshold=     ") 
outputfilee2= open("inputs/en_th2.dat", "w")
outputfilee2.write("%8.5f" % eth2)


eth3=input("3 phonon energy threshold=   ") 

outputfilee3= open("inputs/en_th3.dat", "w")
outputfilee3.write("%8.5f" % eth3 )

print ('Creation of the bash file for running the job')
print ('***********')
outputfile4= open("run.sh", "w")
########################                                                                                                                        
outputfile4.write('source mpi_openmp_set.sh\n')
outputfile4.write('# HFB code\n')
outputfile4.write('\n')
outputfile4.write('echo "HF calculation" \n')
outputfile4.write('cd hf\n')
outputfile4.write('./Hf > log_HF\n')
outputfile4.write('\n')                                                    
outputfile4.write('echo "F matrix calculation"\n')
outputfile4.write('cd ../fmat\n')
outputfile4.write('./Fmat < input_space > log_fmat\n')
outputfile4.write('echo "TDA calculation"\n')
outputfile4.write('cd ../tda\n')
outputfile4.write('./Tda > log_TDA\n')
outputfile4.write('cd ..\n')
outputfile4.write('\n')                                                         
outputfile4.write('echo "1-phonon densities calculation"\n')
outputfile4.write('./run_dens1.sh \n')                                                      
outputfile4.write('\n') 
outputfile4.write('echo "Phonon interaction calculation"\n')
outputfile4.write('./run_phon_int.sh    \n')
outputfile4.write('\n')                                                    
outputfile4.write('echo "EMPM 2-phonon calculation"\n')
outputfile4.write('./run_eqm2.sh\n')
outputfile4.write('\n')                                                        
outputfile4.write('echo "EMPM 3-phonon calculation"\n')
outputfile4.write('./run_eqm3.sh\n')
outputfile4.close()

outputfile5= open("run_admat2.sh", "w")
outputfile5.write('echo "Calculation of 2-phonon AD matrices"\n')
outputfile5.write('#cat $PBS_NODEFILE > nodes.txt\n')
#outputfile5.write('#/opt/intel/compilers_and_libraries_2020/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./eqm_admat > log_ad 2>error_ad\n')
#outputfile5.write('#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./eqm_admat > log_ad 2>error_ad\n')
#outputfile5.write('#source ../mpi_openmp_set.sh\n')
outputfile5.write('$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./eqm_admat > log_ad 2>error_ad\n')
outputfile5.close()

outputfile6= open("run_admat3.sh", "w")
outputfile6.write('echo "Calculation of AD matrices"\n')
outputfile6.write('#cat $PBS_NODEFILE > nodes.txt\n')
#outputfile6.write('/opt/intel/compilers_and_libraries_2020/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./eqm3_admat > log_ad 2>error\n')
#outputfile6.write('#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./eqm3_admat > log_ad 2>error_ad\n')
outputfile6.write('$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./eqm3_admat > log_ad 2>error_ad\n')
outputfile6.close()

outputfile7= open("run_dens1.sh", "w")
outputfile7.write('#cat $PBS_NODEFILE > nodes.txt\n')
#outputfile7.write('#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error\n')
#outputfile7.write('cd phon_dens1/\n')
outputfile7.write('cd eqm2_MPI/\n')
outputfile7.write('source ../mpi_openmp_set.sh\n')
outputfile7.write('$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_dens1 2>error\n')
outputfile7.close()

outputfile8= open("run_dens2.sh", "w")
outputfile8.write('echo "Calculation of 2-phonon densities"\n')
outputfile8.write('#cat $PBS_NODEFILE > nodes.txt\n')
#outputfile8.write('/opt/intel/compilers_and_libraries_2020/linux/mpi/intel64/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens2_MPI > log_dens2 2>error_dens2\n')
#outputfile8.write('#/software/openmpi/3.1.2/intel/bin/mpirun -np $NUM_MPI_PROCS phon_dens2_MPI > log_dens2 2>error_dens2     \n')
outputfile8.write('$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS phon_dens2_MPI > log_dens2 2>error_dens2     \n')
outputfile8.close()

outputfile9= open("run_eqm2.sh", "w")
#outputfile9.write('#export OMP_NUM_THREADS=24\n')
#outputfile9.write('#echo "# of MPI procs"\n')
#outputfile9.write('#export NUM_MPI_PROCS=24\n')
#outputfile9.write('#echo $NUM_MPI_PROCS\n')
outputfile9.write('cd eqm2_MPI/\n')
outputfile9.write('source ../mpi_openmp_set.sh\n')
outputfile9.write('echo "Energy threshold"\n')
outputfile9.write('cat en_trun.dat\n')
outputfile9.write('./eqm_svd > log_eqm2 2>error_eqm2\n')
outputfile9.close()


outputfile10= open("run_eqm3.sh", "w")
#outputfile10.write('#export OMP_NUM_THREADS=24\n')
#outputfile10.write('#echo "# of MPI procs"\n')
#outputfile10.write('#export NUM_MPI_PROCS=24\n')
#outputfile10.write('#echo $NUM_MPI_PROCS\n')
outputfile10.write('cd eqm3_MPI/\n')
outputfile10.write('source ../mpi_openmp_set.sh\n')
outputfile10.write('echo "Energy threshold"\n')
outputfile10.write('cat en_trun.dat\n')
outputfile10.write('./eqm3_svd > log_eqm3 2>error_eqm3\n')
outputfile10.close()

outputfile11= open("run_phon_int.sh","w")
outputfile11.write('#  ******** IMPORTANT **************************************************************************  \n')
outputfile11.write('#  MPI version crashes for specific number of MPI porcesses, try to use smaller number if possible \n')
outputfile11.write('cd eqm2_MPI/\n')
outputfile11.write('source ../mpi_openmp_set.sh\n')
outputfile11.write('#$MPIRUN_PATH/mpirun -np $NUM_MPI_PROCS ./phon_int_MPI > log_phon_int 2>error    \n')
outputfile11.write('$MPIRUN_PATH/mpirun -np 4 ./phon_int_MPI > log_phon_int 2>error\n')
outputfile11.close()

print ('Creations of the folders and copy of the inputs files')
print ('***********')
#Creations of the folders and copy of the inputs files

os.system('mkdir eqm2_MPI')
os.system('mkdir eqm3_MPI')
os.system('mkdir fmat')
os.system('mkdir fullham')
os.system('mkdir nondiag')
#os.system('mkdir phon_dens1')
#os.system('mkdir phon_int')
os.system('mkdir scratch')
os.system('mkdir tda')
os.system('mkdir hf')


os.system('ln -s inputs/mpi_openmp_set.sh')
###############
###input#####
###############
os.system('cp inputs/input.dat hf/')
os.system('ln -s ../vlk.dat hf/')


###############
#####fmat#########
###############

os.system('cp inputs/input_tda_coup.dat tda/')

os.system('ln -s ../fmat/fmat_n.dat tda/')
os.system('ln -s ../fmat/fmat_p.dat tda/')
os.system('ln -s ../fmat/fmat_pn.dat tda/')
os.system('ln -s ../fmat/fmat_np.dat tda/')
os.system('ln -s ../fmat/singpart_coup.dat tda/')
os.system('ln -s ../hf/r1Y1_n.dat tda/')
os.system('ln -s ../hf/r1Y1_p.dat tda/')


outputfile2= open("fmat/Input_files", "w")

outputfile2.write('! Name of interaction file\n')
outputfile2.write('        vlk_hfb.dat\n')
outputfile2.write('! Name of output proton interaction file\n')
outputfile2.write('        fmat_p.dat\n ')
outputfile2.write('! Name of output neutron interaction file\n')
outputfile2.write('        fmat_n.dat\n' )
outputfile2.write('! Name of output proton-neutron interaction file\n')
outputfile2.write('        fmat_pn.dat\n')
outputfile2.write('! Name of output neutron-proton interaction file\n')
outputfile2.write('        fmat_np.dat\n')

outputfile2.close()
os.system('cp inputs/input_space fmat/')
os.system('ln -s ../hf/HF_n.dat fmat/')
os.system('ln -s ../hf/HF_p.dat fmat/')
os.system('ln -s ../hf/vlk_hfb.dat fmat/')


###############
####tda#######
###############
os.system('mkdir tda/1phonon')
outputfile1= open("tda/Input_files_tda", "w")


outputfile1.write('! Name of s.p. file\n')
outputfile1.write('\n')
outputfile1.write('! Name of interaction file\n')
outputfile1.write('\n')
outputfile1.write('! Name of output F proton interaction file\n')
outputfile1.write('        fmat_p.dat\n')
outputfile1.write('! Name of output F neutron interaction file\n')
outputfile1.write('        fmat_n.dat\n')
outputfile1.write('! Name of output F proton-neutron interaction file\n')
outputfile1.write('        fmat_pn.dat\n')
outputfile1.write('! Name of output V interaction file\n')
outputfile1.write('        vlk_nat_orb.dat\n')

outputfile1.close()

###############
######phon_dens1#####
###############
#os.system('ln -s ../inputs/input_tda_coup.dat phon_dens1/')
#os.system('ln -s ../fmat/signpart_coup.dat phon_dens1/')
#os.system('ln -s ../scratch phon_dens1/')
#os.system('ln -s ../tda/1phonon phon_dens1/')

###############
######phon_int#####
###############
#os.system('ln -s ../fmat/fmat_n.dat phon_int/')
#os.system('ln -s ../fmat/fmat_p.dat phon_int/')
#os.system('ln -s ../fmat/fmat_pn.dat phon_int/')
#os.system('ln -s ../fmat/fmat_np.dat phon_int/')
#os.system('ln -s ../fmat/singpart_coup.dat phon_int/')
#os.system('ln -s ../tda/1phonon phon_int/')
#os.system('ln -s ../scratch phon_int/')
#os.system('ln -s ../inputs/input_tda_coup.dat phon_int/')


###############
######eqm2_MPI#
###############

os.system('ln -s ../tda/1phonon eqm2_MPI/')
os.system('ln -s ../scratch eqm2_MPI/')
os.system('ln -s ../inputs/input_tda_coup.dat eqm2_MPI/')
os.system('ln -s ../fmat/singpart_coup.dat eqm2_MPI/')
os.system('ln -s ../tda/tda_r_overl.dat eqm2_MPI/')
os.system('ln -s ../inputs/mpi_openmp_set.sh eqm2_MPI/')
os.system('ln -s ../inputs/en_th2.dat eqm2_MPI/en_trun.dat')
os.system('mkdir eqm2_MPI/2phonon')

###############
######eqm3_MPI#
###############

os.system('ln -s ../tda/1phonon eqm3_MPI/')
os.system('ln -s ../eqm2_MPI/2phonon eqm3_MPI/')
os.system('mkdir eqm3_MPI/3phonon')
os.system('ln -s ../inputs/en_th3.dat eqm3_MPI/en_trun.dat')
os.system('ln -s ../eqm2_MPI/h_corr.dat eqm3_MPI/')
os.system('ln -s ../tda/input_tda_coup.dat eqm3_MPI/')
os.system('ln -s ../run_admat3.sh eqm3_MPI/')
os.system('ln -s ../run_dens2.sh eqm3_MPI/')
os.system('ln -s ../fmat/singpart_coup.dat eqm3_MPI/')
os.system('ln -s ../tda/tda_r_overl.dat eqm3_MPI/')
os.system('ln -s ../scratch eqm3_MPI/')

###############
######nondiag##
###############
os.system('ln -s ../tda/1phonon nondiag/')
os.system('ln -s ../eqm2_MPI/2phonon nondiag/')
os.system('ln -s ../fmat/fmat_n.dat nondiag/')
os.system('ln -s ../fmat/fmat_p.dat nondiag/')
os.system('ln -s ../fmat/fmat_pn.dat nondiag/')
os.system('ln -s ../fmat/fmat_np.dat nondiag/')
os.system('ln -s ../fmat/singpart_coup.dat nondiag/')
os.system('ln -s ../inputs/input_tda_coup.dat nondiag/')
os.system('ln -s ../scratch nondiag/')

outputfile1= open("nondiag/input", "w")
outputfile1.write('-1\n')
outputfile1.write('1\n')
outputfile1.close()


###############
######fullham##
###############
os.system('ln -s ../tda/1phonon fullham/')
os.system('ln -s ../eqm2_MPI/2phonon fullham/')
os.system('ln -s ../eqm3_MPI/3phonon fullham/')
os.system('ln -s ../nondiag/Vint_phon12.dat fullham/')
os.system('ln -s ../nondiag/Vint_phon23.dat fullham/')
os.system('ln -s ../nondiag/Vint_phon13.dat fullham/')




#############################################
###### Move executable files from run directory##
#############################################
btda='mv '+bindir+'/Tda tda/' 
bfmat='mv '+bindir+'/Fmat fmat/'
binput='mv '+bindir+'/Hf hf/'
bdens1='mv '+bindir+'/phon_dens1_MPI eqm2_MPI/'
bint='mv '+bindir+'/phon_int_MPI eqm2_MPI/'
beq2a='mv '+bindir+'/eqm_admat eqm2_MPI/'
beq2s='mv '+bindir+'/eqm_svd eqm2_MPI/'
beq3a='mv '+bindir+'/eqm3_admat eqm3_MPI/'
beq3s='mv '+bindir+'/eqm3_svd eqm3_MPI/'
beq3p='mv '+bindir+'/phon_dens2_MPI eqm3_MPI/'
bnond='mv '+bindir+'/Nondiag nondiag/'
bfullham='mv '+bindir+'/Fullham fullham/'
#print btda
os.system(btda)
os.system(bfmat)
os.system(binput)
os.system(binput)
os.system(bdens1)
os.system(bint)
os.system(beq2a)
os.system(beq2s)
os.system(beq3a)
os.system(beq3s)
os.system(beq3p)
os.system(bnond)
os.system(bfullham)


os.system('chmod 777 run*')
sys.exit()    
