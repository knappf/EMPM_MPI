      subroutine hfb_iteration

       USE technical
       USE geom
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision, allocatable, save :: argu2(:)
       double precision, allocatable, save :: argup(:)
       double precision, allocatable, save :: argun(:)
       double precision, allocatable, save :: rad_denp(:)
       double precision, allocatable, save :: rad_denn(:)
       double precision, allocatable, save :: rad_den2(:)

       double precision :: hp(id,id),Dp(id,id)
       double precision :: hn(id,id),Dn(id,id)
       double precision :: wp(id),wn(id)

       double precision :: rp(id,id),rn(id,id)
       double precision :: kp(id,id),kn(id,id)

!       double precision :: h2p(id,id),h2n(id,id)
       double precision :: w2p(id),w2n(id)

!       double precision :: D2p(id,id),D2n(id,id)

       double precision :: E_diff, E_orig
       integer :: i_count1,i_count2,count_max

       alpha=0.90d0

       count_max=2000
       E_diff=1.d9

       call HFB_energy

       call hp_field(hp)
       call hn_field(hn)
       call Dp_field(Dp)
       call Dn_field(Dn)

       allocate(H11p(id,id),H11n(id,id))
       H11p=Dp
       H11n=Dn

       write(*,*) 'E_HFB,       lambda_p,         lambda_n'
       write(*,*) E_HFB,ferp,fern

       i_count1=0
       E_diff1=1.d9

       do while(((ifp_hfb.or.ifn_hfb).and.(dabs(E_diff1).gt.precis
     &         .and.i_count1.le.2)).or.((.not.(ifp_hfb.or.ifn_hfb)).and.
     &         (dabs(E_diff1).gt.precis.and.i_count1.le.count_max)))
        i_count1=i_count1+1
        E_orig1=E_HFB

        call diagonalization(hp,wp,id)
        call diagonalization(hn,wn,id)

        call make_sp_levels(hp,hn,wp,wn,i_count1)
        tran_p=hp
        tran_n=hn

        call make_densities
        call check_densities

        bos1=dsqrt(0.5d0*(mp+mn)*hbarom/(hbarc**2.d0))
        rad_den1=0.d0
        do i=1,igrid
         val=0.d0
         radi=zcross(i)  !dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          val=val+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
          val=val+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
         enddo
         rad_den1(i)=val/(4.d0*pi)
        enddo
!********************************************************************************
       r4tot=0.d0
       r2tot=0.d0
       r4Z=0.d0
       r2Z=0.d0
       r4N=0.d0
       r2N=0.d0
       allocate(argu2(0:igrid2),argup(0:igrid2),argun(0:igrid2))
       argu2=0.d0
       argup=0.d0
       argun=0.d0
       allocate(rad_denp(0:igrid2),rad_denn(0:igrid2),
     &                                    rad_den2(0:igrid2))
       rad_denp=0.d0 
       rad_denn=0.d0
       rad_den2=0.d0
        do i=0,igrid2
         val1=0.d0
         val2=0.d0
         radi=dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          val1=val1+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
          val2=val2+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
         enddo
         rad_denp(i)=val1/(4.d0*pi)
         rad_denn(i)=val2/(4.d0*pi)
         rad_den2(i)=(val1+val2)/(4.d0*pi)
        enddo
       do i=0,igrid2
        radi=dble(i)*sizebox/dble(igrid2)
        argu2(i)=radi**2.d0*rad_den2(i)
        argup(i)=radi**2.d0*rad_denp(i)
        argun(i)=radi**2.d0*rad_denn(i)
       enddo
       dx=sizebox/dble(igrid2)
       do i=0,igrid2-1
        radi=(dble(i)+0.5d0)*sizebox/dble(igrid2)
        r4tot=r4tot+dx*0.5d0*(argu2(i)+argu2(i+1))*radi**2.d0
        r2tot=r2tot+dx*0.5d0*(argu2(i)+argu2(i+1))
        r4Z=r4Z+dx*0.5d0*(argup(i)+argup(i+1))*radi**2.d0
        r2Z=r2Z+dx*0.5d0*(argup(i)+argup(i+1))
        r4N=r4N+dx*0.5d0*(argun(i)+argun(i+1))*radi**2.d0
        r2N=r2N+dx*0.5d0*(argun(i)+argun(i+1))
       enddo

       open(1,file='radial.out',status='unknown',form='formatted')
        write(1,*) 'sqrt(<r^2>_tot)=',dsqrt(r4tot/r2tot),'fm'
        write(1,*) 'sqrt(<r^2>_p)=',dsqrt(r4Z/r2Z),'fm'
        write(1,*) 'sqrt(<r^2>_n)=',dsqrt(r4N/r2N),'fm'
       close(1)
       deallocate(argu2,argup,argun)
       deallocate(rad_denp,rad_denn,rad_den2)
!********************************************************************************
!       call precal_C1_C2(jmax,cgg1_int,cgg2_int)

       if (if_dd.eq.1) then 

       
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
        
        
!   nacitanie novych V_DD

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
             if (tri.ne.0.0d0) then 
             valp=valp+dble(2*Jp+1)*phase*anorm*tri*Vpp_DD(i,k,j,l,Jp)
             valn=valn+dble(2*Jp+1)*phase*anorm*tri*Vnn_DD(i,k,j,l,Jp)
             valpn=valpn+dble(2*Jp+1)*phase*anorm*tri*Vpn_DD(i,k,j,l,Jp)
             endif 
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

        call HFB_energy

        E_diff1=E_HFB-E_orig1

        write(*,*) E_HFB,ferp,fern

        call hp_field(hp)
        call hn_field(hn)
        call Dp_field(Dp)
        call Dn_field(Dn)

        H11p=Dp
        H11n=Dn
       enddo

       if(ifp_hfb.or.ifn_hfb) then
        i_count2=0
        E_diff2=1.d9

        do while(dabs(E_diff2).gt.precis.and.i_count2.le.count_max)
         i_count2=i_count2+1
         E_orig2=E_HFB

         rp=rhop_HFB
         rn=rhon_HFB
         kp=kapp_HFB
         kn=kapn_HFB

         if(.not.ifp_hfb) then
          call diagonalization(hp,wp,id)

          call make_sp_levels(hp,hn,wp,wn,i_count2+2)
          tran_p=hp
         endif

         if(.not.ifn_hfb) then
          call diagonalization(hn,wn,id)

          call make_sp_levels(hp,hn,wp,wn,i_count2+2)
          tran_n=hn
         endif

         if(ifp_hfb) then
          if(i_count2.eq.1) then
           do i=1,id
            do j=1,id
             if(levp(i)%j2.eq.levp(j)%j2.and.levp(i)%l.eq.levp(j)%l)
     &                                                            then
              Dp(i,j)=Dp(i,j)+0.1d0
             endif
            enddo
           enddo
          endif

          do i=1,id
           hp(i,i)=hp(i,i)-ferp
          enddo

!          call dgemm('N','N',id,id,id,1.d0,hp,max(1,id),hp,
!     &                        max(1,id),0.d0,h2p,max(1,id))
!          call dgemm('N','T',id,id,id,1.d0,Dp,max(1,id),Dp,
!     &                        max(1,id),0.d0,D2p,max(1,id))
!          do i=1,id
!           do j=1,id
!            h2p(i,j)=h2p(i,j)+D2p(i,j)
!           enddo
!          enddo

!          call diagonalization(h2p,w2p,id)
          call diagonalization(hp,wp,id)

          do i=1,id
           w2p(i)=Dp(i,i)
           w2n(i)=Dn(i,i)
          enddo

          call make_qsp_levels(wp,wn,hp,hn,w2p,w2n)
          tran_p=hp
         endif

         if(ifn_hfb) then
          if(i_count2.eq.1) then
           do i=1,id
            do j=1,id
             if(levn(i)%j2.eq.levn(j)%j2.and.levn(i)%l.eq.levn(j)%l)
     &                                                            then
              Dn(i,j)=Dn(i,j)+0.1d0
             endif
            enddo
           enddo
          endif

          do i=1,id
           hn(i,i)=hn(i,i)-fern
          enddo

!          call dgemm('N','N',id,id,id,1.d0,hn,max(1,id),hn,
!     &                        max(1,id),0.d0,h2n,max(1,id))
!          call dgemm('N','T',id,id,id,1.d0,Dn,max(1,id),Dn,
!     &                        max(1,id),0.d0,D2n,max(1,id))
!          do i=1,id
!           do j=1,id
!            h2n(i,j)=h2n(i,j)+D2n(i,j)
!           enddo
!          enddo

!          call diagonalization(h2n,w2n,id)
          call diagonalization(hn,wn,id)

          do i=1,id
           w2p(i)=Dp(i,i)
           w2n(i)=Dn(i,i)
          enddo

          call make_qsp_levels(wp,wn,hp,hn,w2p,w2n)
          tran_n=hn
         endif

         call make_densities
         call check_densities

         bos1=dsqrt(0.5d0*(mp+mn)*hbarom/(hbarc**2.d0))
         rad_den1=0.d0
         do i=1,igrid
          val=0.d0
          radi=zcross(i)    !dble(i)*sizebox/dble(igrid)
          do j=1,id
           val_p=0.d0
           val_n=0.d0
           do k=1,id
            val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
            val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
           enddo
           val=val+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
           val=val+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
          enddo
          rad_den1(i)=val/(4.d0*pi)
         enddo
!********************************************************************************
       r4tot=0.d0
       r2tot=0.d0
       r4Z=0.d0
       r2Z=0.d0
       r4N=0.d0
       r2N=0.d0
       allocate(argu2(0:igrid2),argup(0:igrid2),argun(0:igrid2))
       argu2=0.d0
       argup=0.d0
       argun=0.d0
       allocate(rad_denp(0:igrid2),rad_denn(0:igrid2),
     &                                    rad_den2(0:igrid2))
       rad_denp=0.d0
       rad_denn=0.d0
       rad_den2=0.d0
        do i=0,igrid2
         val1=0.d0
         val2=0.d0
         radi=dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_p(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_n(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          val1=val1+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
          val2=val2+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
         enddo
         rad_denp(i)=val1/(4.d0*pi)
         rad_denn(i)=val2/(4.d0*pi)
         rad_den2(i)=(val1+val2)/(4.d0*pi)
        enddo
       do i=0,igrid2
        radi=dble(i)*sizebox/dble(igrid2)
        argu2(i)=radi**2.d0*rad_den2(i)
        argup(i)=radi**2.d0*rad_denp(i)
        argun(i)=radi**2.d0*rad_denn(i)
       enddo
       dx=sizebox/dble(igrid2)
       do i=0,igrid2-1
        radi=(dble(i)+0.5d0)*sizebox/dble(igrid2)
        r4tot=r4tot+dx*0.5d0*(argu2(i)+argu2(i+1))*radi**2.d0
        r2tot=r2tot+dx*0.5d0*(argu2(i)+argu2(i+1))
        r4Z=r4Z+dx*0.5d0*(argup(i)+argup(i+1))*radi**2.d0
        r2Z=r2Z+dx*0.5d0*(argup(i)+argup(i+1))
        r4N=r4N+dx*0.5d0*(argun(i)+argun(i+1))*radi**2.d0
        r2N=r2N+dx*0.5d0*(argun(i)+argun(i+1))
       enddo

       open(1,file='radial.out',status='unknown',form='formatted')
        write(1,*) 'sqrt(<r^2>_tot)=',dsqrt(r4tot/r2tot),'fm'
        write(1,*) 'sqrt(<r^2>_p)=',dsqrt(r4Z/r2Z),'fm'
        write(1,*) 'sqrt(<r^2>_n)=',dsqrt(r4N/r2N),'fm'
       close(1)
       deallocate(argu2,argup,argun)
       deallocate(rad_denp,rad_denn,rad_den2)
!********************************************************************************


       if (if_dd.eq.1) then 
       
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
!            C3=cg_int(levp(i)%j2,1,levp(j)%j2,-1,2*Jp)
!            C4=cg_int(levp(k)%j2,1,levp(l)%j2,-1,2*Jp)
            
            if(Jp.ne.0) C1=cgg1_int(levp(i)%j2,levp(j)%j2,Jp)
            if(Jp.ne.0) C2=cgg1_int(levp(k)%j2,levp(l)%j2,Jp)

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

            if (dabs(VpnDD).gt.1.d-10) write(1)int(0,1),int(2*Jp,1),
     &      int(2*i,2),int(2*j,2),int(2*k,2),int(2*l,2),VpnDD
     
           if (dabs(VppDD).gt.1.d-10) write(1)int(-1,1),int(2*Jp,1),
     &      int(2*i,2),int(2*j,2),int(2*k,2),int(2*l,2),VppDD

            if (dabs(VnnDD).gt.1.d-10) write(1)int(1,1),int(2*Jp,1),
     &      int(2*i,2),int(2*j,2),int(2*k,2),int(2*l,2),VnnDD



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
         
         
!  nacitanie novych V_DD

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

         rhop_HFB=alpha*rp+(1.d0-alpha)*rhop_HFB
         rhon_HFB=alpha*rn+(1.d0-alpha)*rhon_HFB
         kapp_HFB=alpha*kp+(1.d0-alpha)*kapp_HFB
         kapn_HFB=alpha*kn+(1.d0-alpha)*kapn_HFB

         call HFB_energy

         E_diff2=E_HFB-E_orig2

         write(*,*) E_HFB,ferp,fern

         call hp_field(hp)
         call hn_field(hn)
         call Dp_field(Dp)
         call Dn_field(Dn)

         H11p=Dp
         H11n=Dn
        enddo
       endif

       open(1,file='Energy.out',status='unknown',form='formatted')
        write(1,*) E_HFB,E_pair
       close(1)

       open(1,file='HF_basis_p.out',status='unknown',form='formatted')
       open(2,file='HF_basis_n.out',status='unknown',form='formatted')
        write(1,*) tran_p 
        write(2,*) tran_n
       close(1)
       close(2)

       do i=1,id
!        lhfp(i)%vi=lhfp(i)%vi*(-1)**((lhfp(i)%j2-1)/2)
!        lhfn(i)%vi=lhfn(i)%vi*(-1)**((lhfn(i)%j2-1)/2)
       enddo
       open(1,file='HF_p.out',status='unknown',form='formatted')
        write(1,*)'i,     l,    2*j,      ei,      qei,     vi'
        do ii = 1, id
         write(1,*) lhfp(ii)%index,lhfp(ii)%l,lhfp(ii)%j2,
     &    lhfp(ii)%ei,lhfp(ii)%qei,lhfp(ii)%vi
        end do
       close(1)

       open(1,file='HF_n.out',status='unknown',form='formatted')
        write(1,*)'i,     l,    2*j,      ei,      qei,     vi'
        do ii = 1, id
         write(1,*) lhfn(ii)%index,lhfn(ii)%l,lhfn(ii)%j2,
     &    lhfn(ii)%ei,lhfn(ii)%qei,lhfn(ii)%vi
        end do
       close(1)

       s1=0.d0
       s2=0.d0
       do i=1,id
        s1=s1+rhop_HFB(i,i)*dble(levp(i)%j2+1)
        s2=s2+rhon_HFB(i,i)*dble(levn(i)%j2+1)
       enddo

       if(dabs(s1-dble(AZ)).gt.dsqrt(precis)) then
        write(*,*) 'Error: wrong proton number AZ!!!'
        write(*,*) s1,'=/=',AZ
        stop
       endif
       if(dabs(s2-dble(AN)).gt.dsqrt(precis)) then
        write(*,*) 'Error: wrong proton number AN!!!'
        write(*,*) s2,'=/=',AN
        stop
       endif

       call check_densities

!       allocate(rad_denp(0:igrid2),rad_denn(0:igrid2),
!     &                                    rad_den2(0:igrid2))
!       rad_denp=0.d0
!       rad_denn=0.d0
!       rad_den2=0.d0
!        do i=0,igrid2
!         val1=0.d0
!         val2=0.d0
!         radi=dble(i)*sizebox/dble(igrid2)
!         do j=1,id
!          val_p=0.d0
!          val_n=0.d0
!          do k=1,id
!           val_p=val_p+tran_p(k,j)
!     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
!           val_n=val_n+tran_n(k,j)
!     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
!          enddo
!          val1=val1+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
!          val2=val2+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
!         enddo
!         rad_denp(i)=val1/(4.d0*pi)
!         rad_denn(i)=val2/(4.d0*pi)
!         rad_den2(i)=(val1+val2)/(4.d0*pi)
!        enddo

!       open(47,file='rad_density.out',status='unknown',form='formatted')
!       rad_sum=0.d0
!       dx=sizebox/dble(igrid2)
!       do i=0,igrid2
!        radi=dble(i)*sizebox/dble(igrid2)
!        rad_sum=rad_sum+rad_den2(i)*dx*radi**2.d0
!        write(47,*) radi,rad_den2(i)
!       enddo
!       write(47,*) 'Int rho = ',rad_sum*4.d0*pi
!       close(47)

!       open(47,file='rad_dens_n.out',status='unknown',form='formatted')
!       rad_sum=0.d0
!       dx=sizebox/dble(igrid2)
!       do i=0,igrid2
!        radi=dble(i)*sizebox/dble(igrid2)
!        rad_sum=rad_sum+rad_denn(i)*dx*radi**2.d0
!        write(47,*) radi,rad_denn(i)
!       enddo
!       write(47,*) 'Int rho = ',rad_sum*4.d0*pi
!       close(47)

!       open(47,file='rad_dens_p.out',status='unknown',form='formatted')
!       rad_sum=0.d0
!       dx=sizebox/dble(igrid2)
!       do i=0,igrid2
!        radi=dble(i)*sizebox/dble(igrid2)
!        rad_sum=rad_sum+rad_denp(i)*dx*radi**2.d0
!        write(47,*) radi,rad_denp(i)
!       enddo
!       write(47,*) 'Int rho = ',rad_sum*4.d0*pi
!       close(47)

!       deallocate(rad_denp,rad_denn,rad_den2)

       return
      end
