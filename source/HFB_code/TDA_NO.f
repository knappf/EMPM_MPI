      subroutine TDA_NO

       USE technical
       USE math
       USE geom
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'


       type(twoquas_type), allocatable, save :: ph(:)
       double precision, allocatable, save :: TDA_matrix(:,:),W(:)
       double precision, allocatable, save :: Hamtot(:,:),wtot(:)
       double precision, allocatable, save :: sixj1(:,:,:,:,:,:)
       double precision, allocatable, save :: spur_CM(:)
       integer :: p1,h1,tz1,p2,h2,tz2

       if(if_QTDA.eq.1) then
        write(*,*) 'TDA_NO runs only for if_QTDA = 0  !!'
        stop
       endif


       call read_inter(4,hbarom,az,an,imax,int(jmax,1)
     &,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp
     &,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn
     &,v_pn,irowc_pn,irowe_pn)

       allocate(Fpp(id,id,id,id))
       allocate(Fnn(id,id,id,id))
       allocate(Fpn(id,id,id,id))
       Fpp=0.d0
       Fnn=0.d0
       Fpn=0.d0

       allocate(sixj1(jmax,jmax,0:jmax,jmax,jmax,0:jmax))
       sixj1=0.d0

       do j1=1,jmax,2
        do j2=1,jmax,2
         do j3=0,jmax
          do j4=1,jmax,2
           do j5=1,jmax,2
            do j6=0,jmax
             sixj1(j1,j2,j3,j4,j5,j6)=sixj_int(j1,j2,2*j3,j4,j5,2*j6)
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
 
       write(*,*) 'TDA Natural Orbitals calculation starts.' 
       
       open(1,file='TDA_NO_energies.out',status='unknown',
     &              form='formatted')
        write(1,*)

       open(2,file='TDA_NO_dimens.out',status='unknown',
     &              form='formatted')
        write(2,*)

       open(3,file='NO_phon_struct.out',status='unknown',
     &              form='formatted')
        write(3,*)

       open(4,file='NO_ph_store.out',status='unknown',
     &              form='unformatted')

       open(5,file='NO_ph_store2.out',status='unknown',
     &              form='unformatted')

       jmax=0
       do i=min_p,max_p
        if(lnop(i)%j2.gt.jmax) jmax=lnop(i)%j2
       enddo
       do i=min_n,max_n
        if(lnon(i)%j2.gt.jmax) jmax=lnon(i)%j2
       enddo

       write(4) min_p,max_p,min_n,max_n

!********** Loop over all multipolarities ******************

       do ipar=-1,1,2
        do Jp=0,4!,jmax

        write(*,*) 'TDA NO calculation for pi=',ipar,'J=',Jp

!********** Selection rules *********************************

       i1=0
        do i=min_p,max_p
         if(dabs(lnop(i)%vi**2.d0-1.d0).lt.1.d-2) then
          do j=min_p,max_p
           if(dabs(lnop(j)%ui**2.d0-1.d0).lt.1.d-2) then 
            if(lnop(i)%ipar*lnop(j)%ipar.eq.ipar) then
             if((lnop(i)%j2+lnop(j)%j2)/2.ge.Jp) then
              if(abs(lnop(i)%j2-lnop(j)%j2)/2.le.Jp) then
               i1=i1+1
              endif  
             endif
            endif
           endif
          enddo
         endif
        enddo 


        do i=min_n,max_n
         if(dabs(lnon(i)%vi**2.d0-1.d0).lt.1.d-2) then
          do j=min_n,max_n
           if(dabs(lnon(j)%ui**2.d0-1.d0).lt.1.d-2) then
            if(lnon(i)%ipar*lnon(j)%ipar.eq.ipar) then
             if((lnon(i)%j2+lnon(j)%j2)/2.ge.Jp) then
              if(abs(lnon(i)%j2-lnon(j)%j2)/2.le.Jp) then
               i1=i1+1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo

       if(i1.gt.0) allocate(ph(i1))

!********** Phonons *****************************************       

       i1=0
        do i=min_p,max_p
         if(dabs(lnop(i)%vi**2.d0-1.d0).lt.1.d-2) then
          do j=min_p,max_p
           if(dabs(lnop(j)%ui**2.d0-1.d0).lt.1.d-2) then  
            if(lnop(i)%ipar*lnop(j)%ipar.eq.ipar) then
             if((lnop(i)%j2+lnop(j)%j2)/2.ge.Jp) then
              if(abs(lnop(i)%j2-lnop(j)%j2)/2.le.Jp) then
               i1=i1+1
               ph(i1)%q1=i
               ph(i1)%q2=j
               ph(i1)%tz=-1  !tz=-1 => protons
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo

        do i=min_n,max_n
         if(dabs(lnon(i)%vi**2.d0-1.d0).lt.1.d-2) then
          do j=min_n,max_n
           if(dabs(lnon(j)%ui**2.d0-1.d0).lt.1.d-2) then
            if(lnon(i)%ipar*lnon(j)%ipar.eq.ipar) then
             if((lnon(i)%j2+lnon(j)%j2)/2.ge.Jp) then
              if(abs(lnon(i)%j2-lnon(j)%j2)/2.le.Jp) then
               i1=i1+1
               ph(i1)%q1=i
               ph(i1)%q2=j
               ph(i1)%tz=+1  !tz=+1 => neutrons
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo

!**********  missing orthogonalization of spurious state ****

        i1_spur=0
        if(if_ort.eq.1.and.(ipar.eq.-1.and.Jp.eq.1)) i1_spur=1

        write(2,*) 'pi=',ipar,'J=',Jp
        write(2,*) 'dimension=',i1-i1_spur

        if(i1.gt.0) allocate(TDA_matrix(i1,i1))
        if(i1.gt.0) allocate(W(i1))
        if(i1.gt.0) TDA_matrix=0.d0

!        if(if_ort.eq.1) then
!        if(ipar.eq.-1.and.Jp.eq.1) allocate(spur_CM(i1))
!        if(ipar.eq.-1.and.Jp.eq.1)
!     &   call spur_vec_NO(if_QTDA,spur_CM,i1,Jp,ph)
!        endif

!********* Transformation from V to F **********************

       if(i1.gt.0) then
        Fpp=0.d0
        Fnn=0.d0
        Fpn=0.d0
        do i=1,id
         do j=1,id
          do k=1,id
           do l=1,id
            valp=0.d0
            valn=0.d0
            valpn=0.d0
            do Jpp=0,jmax
             phasep=(-1)**((lnop(j)%j2+lnop(k)%j2)/2-Jp-Jpp)
             phasen=(-1)**((lnon(j)%j2+lnon(k)%j2)/2-Jp-Jpp)
             phasepn=(-1)**((lnop(j)%j2+lnon(k)%j2)/2-Jp-Jpp)
             trip=dble(2*Jpp+1)*sixj1(lnop(i)%j2,lnop(j)%j2,Jp,
     &            lnop(l)%j2,lnop(k)%j2,Jpp)
             trin=dble(2*Jpp+1)*sixj1(lnon(i)%j2,lnon(j)%j2,Jp,
     &            lnon(l)%j2,lnon(k)%j2,Jpp)
             tripn=dble(2*Jpp+1)*sixj1(lnop(i)%j2,lnop(j)%j2,Jp,
     &            lnon(l)%j2,lnon(k)%j2,Jpp)
             valp=valp+phasep*trip*Vpp(i,k,j,l,Jpp)
             valn=valn+phasen*trin*Vnn(i,k,j,l,Jpp)
             valpn=valpn+phasepn*tripn*Vpn(i,k,j,l,Jpp)
            enddo
            Fpp(i,j,k,l)=valp
            Fnn(i,j,k,l)=valn
            Fpn(i,j,k,l)=valpn
            if(ipar.eq.+1.and.Jp.eq.0) then
             Fpp(i,j,k,l)=F_0*valp
             Fnn(i,j,k,l)=F_0*valn
             Fpn(i,j,k,l)=F_0*valpn
            endif
           enddo
          enddo
         enddo
        enddo

!********* Construction of the TDA matrix ******************
                                                 
        do i=1,i1                      
         p1=ph(i)%q2                  
         h1=ph(i)%q1                       
         tz1=ph(i)%tz                
        do j=1,i1                               
         p2=ph(j)%q2                                 
         h2=ph(j)%q1                                
         tz2=ph(j)%tz                           
                                                    
        if(tz1.eq.1.and.tz2.eq.-1) then                       
        TDA_matrix(i,j)=TDA_matrix(i,j)                  
     &                  +(-1)**(Jp+(lnop(p2)%j2-lnop(h2)%j2)/2)
     &                  *Fpn(h2,p2,p1,h1)  
        endif !neutron-proton subblock                        
                                                           
        if(tz1.eq.-1.and.tz2.eq.1) then                 
        TDA_matrix(i,j)=TDA_matrix(i,j)                  
     &                  +(-1)**(Jp-(lnon(p2)%j2-lnon(h2)%j2)/2)
     &                  *Fpn(p1,h1,h2,p2)      
        endif !proton-neutron subblock                            
         


        if(tz1.eq.1.and.tz2.eq.1) then                             

        TDA_matrix(i,j)=TDA_matrix(i,j)
     &                 +(-1)**(Jp+(lnon(h1)%j2-lnon(p1)%j2)/2)
     &                 *Fnn(p2,h2,h1,p1)
 
        if(h1.eq.h2) then
        TDA_matrix(i,j)=TDA_matrix(i,j)+kin_n(p2,p1)
        endif

        if(p1.eq.p2) then
        TDA_matrix(i,j)=TDA_matrix(i,j)-kin_n(h1,h2)
        endif

        if(h1.eq.h2) then
         if(lnon(p1)%j2.eq.lnon(p2)%j2) then
          do jang=0,jmax!Jp
           do ihole=min_n,max_n  
            if(dabs(lnon(ihole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                      +dble(2*jang+1)/dble(lnon(p1)%j2+1)
     &                      *Vnn(p2,ihole,p1,ihole,jang)
            endif
           enddo
          enddo 
         endif
        endif

        if(h1.eq.h2) then
         if(lnon(p1)%j2.eq.lnon(p2)%j2) then
          do jjang=0,jmax!Jp
           do iihole=min_p,max_p
            if(dabs(lnop(iihole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                      +dble(2*jjang+1)/dble(lnon(p1)%j2+1)
     &                      *Vpn(iihole,p2,iihole,p1,jjang)
            endif
           enddo
          enddo    
         endif
        endif
         
        if(p1.eq.p2) then
         if(lnon(h1)%j2.eq.lnon(h2)%j2) then
          do jjjang=0,jmax!Jp
           do iiihole=min_n,max_n
            if(dabs(lnon(iiihole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                      -dble(2*jjjang+1)/dble(lnon(h1)%j2+1)
     &                      *Vnn(h1,iiihole,h2,iiihole,jjjang) 
            endif 
           enddo
          enddo
         endif
        endif

        if(p1.eq.p2) then
         if(lnon(h1)%j2.eq.lnon(h2)%j2) then
          do jjjjang=0,jmax!Jp
           do iiiihole=min_p,max_p
            if(dabs(lnop(iiiihole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                      -dble(2*jjjjang+1)/dble(lnon(h1)%j2+1)
     &                      *Vpn(iiiihole,h1,iiiihole,h2,jjjjang)

            endif 
           enddo
          enddo
         endif        
        endif

        endif !neutron-neutron subblock                     
                                                               
        if(tz1.eq.-1.and.tz2.eq.-1) then                
        
        TDA_matrix(i,j)=TDA_matrix(i,j)
     &                 +(-1)**(Jp+(lnop(h1)%j2-lnop(p1)%j2)/2)
     &                 *Fpp(p2,h2,h1,p1)

        if(h1.eq.h2) then
        TDA_matrix(i,j)=TDA_matrix(i,j)+kin_p(p2,p1)
        endif

        if(p1.eq.p2) then
        TDA_matrix(i,j)=TDA_matrix(i,j)-kin_p(h1,h2)
        endif
       
        if(h1.eq.h2) then
         if(lnop(p1)%j2.eq.lnop(p2)%j2) then
          do jpang=0,jmax!Jp
           do iphole=min_p,max_p
            if(dabs(lnop(iphole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                      +dble(2*jpang+1)/dble(lnop(p1)%j2+1)
     &                      *Vpp(p2,iphole,p1,iphole,jpang)
            endif
           enddo
          enddo
         endif
        endif

        if(h1.eq.h2) then
         if(lnop(p1)%j2.eq.lnop(p2)%j2) then
          do jjpang=0,jmax!Jp
           do iiphole=min_n,max_n
            if(dabs(lnon(iiphole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                      +dble(2*jjpang+1)/dble(lnop(p1)%j2+1)
     &                      *Vpn(p2,iiphole,p1,iiphole,jjpang)
            endif
           enddo
          enddo
         endif                                          
        endif
       
        if(p1.eq.p2) then
         if(lnop(h1)%j2.eq.lnop(h2)%j2) then
          do jjjpang=0,jmax!Jp
           do iiiphole=min_p,max_p
            if(dabs(lnop(iiiphole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                       -dble(2*jjjpang+1)/dble(lnop(h1)%j2+1)
     &                       *Vpp(h1,iiiphole,h2,iiiphole,jjjpang)
            endif
           enddo
          enddo
         endif
        endif
 
        if(p1.eq.p2) then
         if(lnop(h1)%j2.eq.lnop(h2)%j2) then
          do jjjjpang=0,jmax!Jp
           do iiiiphole=min_n,max_n
            if(dabs(lnon(iiiiphole)%vi**2.d0-1.d0).lt.1.d-2) then
             TDA_matrix(i,j)=TDA_matrix(i,j)
     &                      -dble(2*jjjjpang+1)/dble(lnop(h1)%j2+1)
     &                      *Vpn(h1,iiiiphole,h2,iiiiphole,jjjjpang)    
            endif 
           enddo 
          enddo
         endif
        endif

        endif !proton-proton subblock                          
                                                                     
        enddo !ends i-loop                          
        enddo !ends j-loop                       
 
!********* Writing output files ****************************
                                              
        if((ipar.eq.1.and.Jp.eq.0).or.(ipar.eq.-1.and.Jp.eq.1)) then
         write(5) i1
         write(5) TDA_matrix
        endif

        if(.not.(if_ort.eq.1.and.(ipar.eq.-1.and.Jp.eq.1))) then

        call diagonalization(TDA_matrix,W,i1)

         write(1,*) 'pi=',ipar,'J=',Jp
         write(1,*) W(1:i1)
        
        call transit_calc_NO(ipar,Jp,ph,TDA_matrix,W,i1)
!        call phonon_density_calc(ipar,Jp,ph,TDA_matrix,W,i1)        

        write(3,*) 'pi=',ipar,'J=',Jp
        write(3,*)
        do i=1,i1
         write(3,*) 'phonon E_i',W(i)
         write(3,*)
         do j=1,i1
          m1=ph(j)%q1
          m2=ph(j)%q2
          m3=ph(j)%tz
          if(TDA_matrix(j,i)**2.d0.gt.1.d-2) then
           if(m3.eq.-1) then
            write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
     &             'p',lnop(m1)%l,lnop(m1)%j2,lnop(m1)%ei,lnop(m2)%l,
     &                 lnop(m2)%j2,lnop(m2)%ei,TDA_matrix(j,i)**2.d0
           endif
           if(m3.eq.+1) then
            write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
     &             'n',lnon(m1)%l,lnon(m1)%j2,lnon(m1)%ei,lnon(m2)%l,
     &                 lnon(m2)%j2,lnon(m2)%ei,TDA_matrix(j,i)**2.d0
           endif
          endif
         enddo     ! ends j=1,i1 loop
         write(3,*) 
        enddo    ! ends i=1,i1 loop

        elseif(ipar.eq.-1.and.Jp.eq.1) then

        call diag_subspace(TDA_matrix,W,spur_CM,i1)
        write(1,*) ' pi=',ipar,'J=',Jp
        write(1,*) w(2:i1)
        call transit_calc_spur_NO(ipar,Jp,ph,TDA_matrix,W,i1)
        write(3,*) ' pi=',ipar,'J=',Jp
        write(3,*)
         do i=2,i1
          write(3,*) 'phonon  E_i=',w(i)
          write(3,*)
          do j=1,i1
           m1=ph(j)%q1
           m2=ph(j)%q2
           m3=ph(j)%tz
           if(TDA_matrix(j,i)**2.d0.gt.1.d-2) then
            if(m3.eq.-1) then
             write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
     &             'p',lnop(m1)%l,lnop(m1)%j2,lnop(m1)%ei,lnop(m2)%l,
     &                 lnop(m2)%j2,lnop(m2)%ei,TDA_matrix(j,i)**2.d0
            endif
            if(m3.eq.+1) then
             write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
     &             'n',lnon(m1)%l,lnon(m1)%j2,lnon(m1)%ei,lnon(m2)%l,
     &                 lnon(m2)%j2,lnon(m2)%ei,TDA_matrix(j,i)**2.d0
            endif
           endif
          enddo
          write(3,*)
         enddo
        endif

        write(4) ipar,Jp,i1
        write(4) (ph(i),i=1,i1)
        write(4) (W(i),i=1,i1)

        do i=1,i1
          write(4) (TDA_matrix(j,i),j=1,i1)
        enddo

!****************************************************************
!       Diagonalization of the Hamiltonian in the 0+1 space
        if(Jp.eq.0.and.ipar.eq.1) then
        allocate(Hamtot(i1+1,i1+1),wtot(i1+1))
        Hamtot=0.d0 
        wtot=0.d0
        do i=1,i1
         Hamtot(i+1,i+1)=W(i)
        enddo
        do i=1,i1
         val=0.d0
         do j=1,i1
          it=ph(j)%tz
          ip=ph(j)%q2
          ih=ph(j)%q1
          if(it.eq.-1) then
           jhp=lnop(ih)%j2
           val=val+TDA_matrix(j,i)*kin_p(ih,ip)*dsqrt(dble(jhp+1))
           do k=1,id
            if(dabs(lnop(k)%vi**2.d0-1.d0).lt.1.d-2) then
             val=val+0.5d0*TDA_matrix(j,i)*Fpp(k,k,ih,ip)
     &                                  *dsqrt(dble(lnop(k)%j2+1))
            endif
            if(dabs(lnon(k)%vi**2.d0-1.d0).lt.1.d-2) then
             val=val+TDA_matrix(j,i)*Fpn(ih,ip,k,k)
     &                                  *dsqrt(dble(lnon(k)%j2+1))
            endif
           enddo
          endif
          if(it.eq.1) then
           jhn=lnon(ih)%j2
           val=val+TDA_matrix(j,i)*kin_n(ih,ip)*dsqrt(dble(jhn+1))
           do k=1,id
            if(dabs(lnop(k)%vi**2.d0-1.d0).lt.1.d-2) then
             val=val+TDA_matrix(j,i)*Fpn(k,k,ih,ip)
     &                                  *dsqrt(dble(lnop(k)%j2+1))
            endif
            if(dabs(lnon(k)%vi**2.d0-1.d0).lt.1.d-2) then
             val=val+0.5d0*TDA_matrix(j,i)*Fnn(k,k,ih,ip)
     &                                  *dsqrt(dble(lnon(k)%j2+1))
            endif
           enddo
          endif
         enddo
         Hamtot(1,1+i)=val
         Hamtot(1+i,1)=val
        enddo
        do i=1,i1+1
         do j=1,i1+1 
          if(dabs(Hamtot(i,j)).gt.1.d-7) write(501,*) i,j,Hamtot(i,j)
         enddo
        enddo
        call diagonalization(Hamtot,wtot,i1+1)
         open(11,file='Hamtot_en.out',status='unknown',form='formatted')
         write(11,*)
         write(11,*) 'pi=',ipar,'J=',Jp
         write(11,*) wtot(1:i1+1)
         close(11)
         open(14,file='Hamtot_store.out',status='unknown',
     &              form='unformatted')
         do i=1,i1+1
          write(14) (Hamtot(j,i),j=1,i1+1)
         enddo
         close(14)
         open(61,file='Ampl.out',status='unknown',form='formatted')
          do ii1=1,i1+1
           da_0=0.d0
           da_1=0.d0
           da_0=Hamtot(1,ii1)**2.d0
            do ii2=2,i1+1
             da_1=da_1+Hamtot(ii2,ii1)**2.d0
            enddo
           write(61,'(1x,4(f12.5,1x))') wtot(ii1),da_0,da_1
          enddo
         close(61)
        deallocate(Hamtot,wtot)
        endif
!****************************************************************

        deallocate(ph,TDA_matrix,W)
         if(if_ort.eq.1) then
         if(ipar.eq.-1.and.Jp.eq.1) deallocate(spur_CM)
         endif
        endif ! ends if condition if(i1.gt.0) then
                                               
        enddo !ends Jp loop                
       enddo !ends ipar loop                 
                                                
       deallocate(sixj1)
       deallocate(Fpp,Fnn,Fpn)                                          
                              
       close(1)
       close(2)
       close(3)                    
       close(4)
       close(5)                           

       write(*,*) 'TDA Natural Orbitals calculation has fininshed.'

        return
       end
