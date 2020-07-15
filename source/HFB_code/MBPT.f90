      subroutine MBPT

       USE technical
       USE math
       USE geom
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
       
       
       
       call read_inter(1,hbarom,az,an,imax,int(jmax,1) &
     &,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp &
     &,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn &
     &,v_pn,irowc_pn,irowe_pn)
     
       
       write(*,*) 'MBPT(2) calculation'

!****************************************************************************
!      The many-body perturbation theory E(2) energy is calculated here     *
       ener2pt=0.d0
!       if((.not.ifp_hfb).and.(.not.ifn_hfb)) then
        call MBPT_energy(ener2pt)

        open(11,file='Energy_2.out',status='unknown',form='formatted')
         write(11,*) 'E^(2)=',ener2pt,'MeV'
        close(11)
!       endif
!****************************************************************************
       return
      end
