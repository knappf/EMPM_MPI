       subroutine interaction

       USE technical
       USE geom
       use int_arr
       USE inter_stor
       

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
!       include 'inter_stor.inc'
       
       

       double precision :: unitar_matrix(id,id)
  
       integer (kind=1) :: iip,iij,iit
       integer (kind=2) :: ii1,ii2,ii3,ii4
       real(kind=8) :: v_elem2,v_elem
!       real (kind=8), dimension (:,:,:), allocatable:: cgg1_int,cgg2_int
       
       character*2 :: ch

       id2=2*id

       jmax=0
       do i=1,id
        if(levp(i)%j2.gt.jmax) jmax=levp(i)%j2
       enddo

       allocate(tran_p(id,id),tran_n(id,id))
       tran_p=0.d0
       tran_n=0.d0
       
       imax=id
       
       
       call read_inter(0,hbarom,az,an,imax,int(jmax,1)
     &,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp
     &,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn
     &,v_pn,irowc_pn,irowe_pn)

!       allocate(Vpp(id,id,id,id,0:jmax))
!       allocate(Vnn(id,id,id,id,0:jmax))
!       allocate(Vpn(id,id,id,id,0:jmax))
!       Vpp=0.d0
!       Vnn=0.d0
!       Vpn=0.d0

       allocate(Fpp(id,id,id,id))
       allocate(Fnn(id,id,id,id))
       allocate(Fpn(id,id,id,id))
       Fpp=0.d0
       Fnn=0.d0
       Fpn=0.d0

!       allocate(Vpp_DD(id,id,id,id,0:jmax))
!       allocate(Vnn_DD(id,id,id,id,0:jmax))
!       allocate(Vpn_DD(id,id,id,id,0:jmax))
!       Vpp_DD=0.d0
!       Vnn_DD=0.d0
!       Vpn_DD=0.d0

!       if (if_dd.eq.1) then 

       allocate(Fpp_DD(id,id,id,id))
       allocate(Fnn_DD(id,id,id,id))
       allocate(Fpn_DD(id,id,id,id))
       Fpp_DD=0.d0
       Fnn_DD=0.d0
       Fpn_DD=0.d0
        
!       endif
       


 
       do i=1,id
        do j=1,id
         do k=1,id
          do l=1,id
   
           valp=0.d0
           valn=0.d0
           valpn=0.d0
       if(levp(i)%j2.eq.levp(j)%j2.and.levp(k)%j2.eq.levp(l)%j2)then
            do Jp=0,jmax

             phase=(-1)**(1+(levp(k)%j2+levp(l)%j2)/2)
             anorm=1.d0/dsqrt(dble(levp(j)%j2+1)*dble(levp(l)%j2+1))
             tri=1.d0
             if(abs((levp(j)%j2-levp(l)%j2)/2).gt.Jp) tri=0.d0
            if((levp(j)%j2+levp(l)%j2)/2.lt.Jp) tri=0.d0
            valp=valp+dble(2*Jp+1)*phase*anorm*tri*Vpp(i,k,j,l,Jp)
            valn=valn+dble(2*Jp+1)*phase*anorm*tri*Vnn(i,k,j,l,Jp)
            valpn=valpn+dble(2*Jp+1)*phase*anorm*tri*Vpn(i,k,j,l,Jp)
            enddo
           endif 
           Fpp(i,j,k,l)=F_0*valp
           Fnn(i,j,k,l)=F_0*valn
           Fpn(i,j,k,l)=F_0*valpn
 
          enddo
         enddo
        enddo
       enddo

!***********************************************************************
!    Here the DD interaction elements are calculated                   *
      
!       igrid=20      ! number of the node points
       igrid2=250     ! number of the grid points
       sizebox=3.0d0*dble(AZ+AN)**0.33333333d0
       bos1=dsqrt(0.5d0*(mp+mn)*hbarom/(hbarc**2.d0))
       dx=sizebox/dble(igrid2)

      
!       unitar_matrix=0.d0
!       do i1=1,id
!        unitar_matrix(i1,i1)=1.d0
!       enddo


       if (if_dd.eq.1) then 
       
       unitar_matrix=0.d0
       do i1=1,id
        unitar_matrix(i1,i1)=1.d0
       enddo

       igrid=20      ! number of the node points
       igrid2=250     ! number of the grid points

!       sizebox=10.d0    ! the interval in which the radial integral will be numerically summed
       sizebox=2.0d0*dble(AZ+AN)**0.33333333d0

       bos1=dsqrt(0.5d0*(mp+mn)*hbarom/(hbarc**2.d0))

       dx=sizebox/dble(igrid2)

       allocate(zcross(igrid))
       zcross(1)=0.070539889692d0 
       zcross(2)=0.372126818002d0
       zcross(3)=0.916582102483d0
       zcross(4)=1.70730653103d0
       zcross(5)=2.74919925531d0
       zcross(6)=4.04892531384d0
       zcross(7)=5.61517497087d0
       zcross(8)=7.45901745389d0
       zcross(9)=9.59439286749d0
       zcross(10)=12.0388025566d0
       zcross(11)=14.8142934155d0
       zcross(12)=17.9488955686d0
       zcross(13)=21.4787881904d0
       zcross(14)=25.4517028094d0
       zcross(15)=29.9325546634d0
       zcross(16)=35.0134341868d0
       zcross(17)=40.8330570974d0
       zcross(18)=47.6199940299d0
       zcross(19)=55.8107957541d0
       zcross(20)=66.5244165252d0
      
       do i=1,igrid
        zcross(i)=dsqrt(zcross(i))/bos1
       enddo

       allocate(weight(igrid))
       weight(1)=0.181080062419d0
       weight(2)=0.422556767879d0
       weight(3)=0.666909546702d0
       weight(4)=0.9153523727d0
       weight(5)=1.1695397071d0
       weight(6)=1.43135498624d0
       weight(7)=1.7029811359d0
       weight(8)=1.98701589585d0
       weight(9)=2.28663576323d0
       weight(10)=2.60583465152d0
       weight(11)=2.94978381794d0
       weight(12)=3.32539569477d0
       weight(13)=3.74225636246d0
       weight(14)=4.21424053477d0
       weight(15)=4.76252016007d0
       weight(16)=5.42172779036d0
       weight(17)=6.25401146407d0
       weight(18)=7.38731523837d0
       weight(19)=9.15132879607d0
       weight(20)=12.8933886244d0
  
       do i=1,igrid
        weight(i)=weight(i)/(2.d0*zcross(i)*bos1**2.d0)
       enddo

       allocate(R_field(0:noscmax/2,0:noscmax,jmax,igrid))
       do nnn1=0,noscmax/2
        do nnn2=0,noscmax
         do nnn3=1,jmax
          do nnn4=1,igrid
           radi=zcross(nnn4)   !dble(nnn4)*sizebox/dble(igrid)
           R_field(nnn1,nnn2,nnn3,nnn4)=R_val(nnn1,nnn2,nnn3,bos1,radi)
          enddo
         enddo
        enddo
       enddo

       allocate(rad_den1(igrid))
       rad_den1=0.d0

       do i=1,igrid
        val=0.d0
        radi=zcross(i)   !dble(i)*sizebox/dble(igrid)
        do j=1,id
         val_p=0.d0
         val_n=0.d0
         do k=1,id
          val_p=val_p+unitar_matrix(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
          val_n=val_n+unitar_matrix(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
         enddo
         val=val+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
         val=val+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
        enddo
        rad_den1(i)=val/(4.d0*pi)
       enddo
       
       call precal_C1_C2(jmax,cgg1_int,cgg2_int)
       
       open(1,file='vlk_dd.dat',status='unknown',form='unformatted')
       
       do i=1,id
        do j=1,id
         do k=1,id
          do l=1,id
           if (100000*i+j.le.100000*k+l.and.
     &     levp(i)%ipar*levp(j)%ipar.eq.levp(k)%ipar*levp(l)%ipar) then 
           RR=R_integral(i,j,k,l,1.d0,rad_den1)
           do Jp=0,jmax
           
            if((levp(i)%j2+levp(j)%j2.ge.2*Jp.and.
     &       abs(levp(i)%j2-levp(j)%j2).le.2*Jp).and.
     &       (levp(k)%j2+levp(l)%j2.ge.2*Jp.and.
     &         abs(levp(k)%j2-levp(l)%j2).le.2*Jp)) then            

           
            jphs=(-1)**((levp(i)%j2+levp(j)%j2+levp(k)%j2+levp(l)%j2)/2)  ! (-1)^(j_i+j_j+j_k+j_l)
            lphs=(-1)**(levp(i)%l+levp(j)%l+levp(k)%l+levp(l)%l)          ! (-1)^(l_i+l_j+l_k+l_l)
            ljphs=(-1)**(levp(i)%l+levp(k)%l+(levp(i)%j2-levp(k)%j2)/2)   ! (-1)^(l_i+l_k+j_i-j_k)
            iphS0=(-1)**(levp(j)%l+levp(l)%l+(levp(j)%j2+levp(l)%j2)/2)   ! (-1)^(l_j+l_l+j_j+j_l)
            ijphs=(-1)**(Jp+levp(i)%l+levp(j)%l)                          ! (-1)^(l_i+l_j+J)
            klphs=(-1)**(Jp+levp(k)%l+levp(l)%l)                          ! (-1)^(l_k+l_l+J)
            C1=0.d0
            C2=0.d0
!            if(Jp.ne.0) C1=cg_int(levp(i)%j2,1,levp(j)%j2,1,2*Jp)
!            if(Jp.ne.0) C2=cg_int(levp(k)%j2,1,levp(l)%j2,1,2*Jp)
             if(Jp.ne.0) C1=cgg1_int(levp(i)%j2,levp(j)%j2,Jp)
             if(Jp.ne.0) C2=cgg1_int(levp(k)%j2,levp(l)%j2,Jp)

!            C3=cg_int(levp(i)%j2,1,levp(j)%j2,-1,2*Jp)
!            C4=cg_int(levp(k)%j2,1,levp(l)%j2,-1,2*Jp)
             C3=cgg2_int(levp(i)%j2,levp(j)%j2,Jp)
             C4=cgg2_int(levp(k)%j2,levp(l)%j2,Jp)


            anorm=1.d0
            if(i.eq.j) anorm=anorm/dsqrt(2.d0)
            if(k.eq.l) anorm=anorm/dsqrt(2.d0)
            dterm=1.d0
            dterm=dterm*dsqrt(dble((levp(i)%j2+1)*(levp(j)%j2+1)
     &                                 *(levp(k)%j2+1)*(levp(l)%j2+1)))
            dterm=dterm/dble(96*pi*(2*Jp+1))
            geometric1=dterm*(2.d0*dble(jphs*(1+lphs))*C1*C2
     &                    +dble(ljphs*(1-ijphs)*(1-klphs))*C3*C4)
            geometric0=dterm*dble(iphS0*(1+ijphs)*(1+klphs))*C3*C4
            
            VpnDD=VDD_S1*RR*geometric1+VDD_S0*RR*geometric0
            VppDD=2.d0*VDD_S0*RR*geometric0*anorm
            VnnDD=2.d0*VDD_S0*RR*geometric0*anorm
            
           if (dabs(VpnDD).gt.1.d-14) write(1)int(0,1),int(2*Jp,1),
     &      int(2*i-1,2),int(2*j,2),int(2*k-1,2),int(2*l,2),VpnDD
    
           if (i.le.j.and.k.le.l) then 
 
           if (dabs(VppDD).gt.1.d-14) write(1)int(-1,1),int(2*Jp,1),
     &      int(2*i-1,2),int(2*j-1,2),int(2*k-1,2),int(2*l-1,2),VppDD

            if (dabs(VnnDD).gt.1.d-14) write(1)int(1,1),int(2*Jp,1),
     &      int(2*i,2),int(2*j,2),int(2*k,2),int(2*l,2),VnnDD

           endif

            
            
!            Vpn_DD(i,j,k,l,Jp)=VDD_S1*RR*geometric1
!     &                                    +VDD_S0*RR*geometric0
!            Vpp_DD(i,j,k,l,Jp)=2.d0*VDD_S0*RR*geometric0*anorm
!            Vnn_DD(i,j,k,l,Jp)=2.d0*VDD_S0*RR*geometric0*anorm
            endif
           enddo
           
           endif
          enddo
         enddo
        enddo
       enddo
       
       close(1)

!***********************************************************************
      call read_inter(2,hbarom,az,an,imax,int(jmax,1)
     &,n_nn_dd,icol_nn_dd,v_nn_dd,irowc_nn_dd,irowe_nn_dd,n_pp_dd
     &,icol_pp_dd,v_pp_dd,irowc_pp_dd,irowe_pp_dd,n_pn_dd,icol_pn_dd
     &,v_pn_dd,irowc_pn_dd,irowe_pn_dd)

       do i=1,id
        do j=1,id
         do k=1,id
          do l=1,id
           valp=0.d0
           valn=0.d0
           valpn=0.d0
           if(levp(i)%j2.eq.levp(j)%j2.and.levp(k)%j2.eq.levp(l)%j2)then
            do Jp=0,jmax
             phase=(-1)**(1+(levp(k)%j2+levp(l)%j2)/2)
             anorm=1.d0/dsqrt(dble(levp(j)%j2+1)*dble(levp(l)%j2+1))
             tri=1.d0
             if(abs((levp(j)%j2-levp(l)%j2)/2).gt.Jp) tri=0.d0
             if((levp(j)%j2+levp(l)%j2)/2.lt.Jp) tri=0.d0
             valp=valp+dble(2*Jp+1)*phase*anorm*tri*Vpp_DD(i,k,j,l,Jp)
             valn=valn+dble(2*Jp+1)*phase*anorm*tri*Vnn_DD(i,k,j,l,Jp)
             valpn=valpn+dble(2*Jp+1)*phase*anorm*tri*Vpn_DD(i,k,j,l,Jp)
            enddo
           endif
           Fpp_DD(i,j,k,l)=F_0*valp
           Fnn_DD(i,j,k,l)=F_0*valn
           Fpn_DD(i,j,k,l)=F_0*valpn
          enddo
         enddo
        enddo
       enddo
       
       endif

       return
      end
