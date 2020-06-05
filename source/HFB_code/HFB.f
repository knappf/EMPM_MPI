      program Phonon_EDF

       USE technical
       USE math

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

        ihf=0
        if_dd=0
        
        call read_input 
            
            
        if (vdd_s1.ne.0.0d0.or.vdd_s0.ne.0.d0) if_dd=1

        call dimm(id,noscmax)  !calculates the dimension 'id' for no. of shells 'noscmax'

        call initialize_basis
        call occupations

        call kinetic
        call interaction
        

        do i=1,id                ! The spin-orbit term
         jjip=levp(i)%j2
         jjin=levn(i)%j2
         llip=levp(i)%l
         llin=levn(i)%l
         kin_p(i,i)=kin_p(i,i)+vso*0.5d0
     &          *(dble(jjip*(jjip+2))/4.d0-dble(llip*(llip+1))-0.75d0)
         kin_n(i,i)=kin_n(i,i)+vso*0.5d0
     &          *(dble(jjin*(jjin+2))/4.d0-dble(llin*(llin+1))-0.75d0)
        enddo
 
        call transition 

        call hfb_iteration

    
!        stop

        call transf_interaction

        call radial_wf(1)

        call TDA

!        if(if_QTDA.eq.1) call phonon_density_calc


        if(if_NAT.eq.1) then 
          call NatOrb
!        if(if_NAT.eq.2) call NatOrb2
          call NatOrb_transf
           ihf=2
           call HFB_energy_NAT
           call radial_wf(2)
           call TDA_NO
        endif    

!        call transit_E0_01
!        call transit_E1_01

!        call deallocate_all

       stop
      end

