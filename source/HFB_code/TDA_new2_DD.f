      subroutine TDA

       USE technical
       USE math
       USE geom
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       type(twoquas_type), allocatable, save :: ph(:)
!       type(phonon_type), allocatable, save :: list_phon(:)
       double precision, allocatable, save :: TDA_matrix(:,:),w(:)
       double precision, allocatable, save :: spur_CM(:),spur_N(:)
       
       
       
       
       call read_inter(1,hbarom,az,an,imax,int(jmax,1)
     &,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp
     &,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn
     &,v_pn,irowc_pn,irowe_pn)
     
        
       
       if(if_QTDA.eq.0) write(*,*) 'TDA calculation starts'
       if(if_QTDA.eq.1) write(*,*) 'QTDA calculation starts'

       open(1,file='TDA_energies.out',status='unknown',form='formatted')
        write(1,*)
       open(2,file='TDA_dimens.out',status='unknown',form='formatted')
        write(2,*)
       open(3,file='phon_struct.out',status='unknown',form='formatted')
        write(3,*)
       open(4,file='ph_store.out',status='unknown',form='unformatted')
 
!       call dimm(id,nosc_TDA)
 
       jmax=0
       do i=min_p,max_p !1,id
        if(lhfp(i)%j2.gt.jmax) jmax=lhfp(i)%j2
       enddo
       do i=min_n,max_n !1,id
        if(lhfn(i)%j2.gt.jmax) jmax=lhfn(i)%j2
       enddo

       write(4) min_p,max_p,min_n,max_n
       
       ihf=1
       do ipar=-1,1,2
        do Jp=0,3!jmax

        write(*,*) 'TDA calculation for pi=',ipar,'J=',Jp

       if(if_QTDA.eq.0) then
        i1=0
        do i=min_p,max_p !1,id
         if(dabs(lhfp(i)%vi**2.d0-1.d0).lt.1.d-2) then 
          do j=min_p,max_p !1,id
           if(dabs(lhfp(j)%ui**2.d0-1.d0).lt.1.d-2) then
            if(lhfp(i)%ipar*lhfp(j)%ipar.eq.ipar) then
             if((lhfp(i)%j2+lhfp(j)%j2)/2.ge.Jp) then
              if(abs(lhfp(i)%j2-lhfp(j)%j2)/2.le.Jp) then
               i1=i1+1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo
          
        do i=min_n,max_n !1,id
         if(dabs(lhfn(i)%vi**2.d0-1.d0).lt.1.d-2) then
          do j=min_n,max_n !1,id
           if(dabs(lhfn(j)%ui**2.d0-1.d0).lt.1.d-2) then
            if(lhfn(i)%ipar*lhfn(j)%ipar.eq.ipar) then
             if((lhfn(i)%j2+lhfn(j)%j2)/2.ge.Jp) then
              if(abs(lhfn(i)%j2-lhfn(j)%j2)/2.le.Jp) then
               i1=i1+1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo
        
        if(i1.gt.0) allocate(ph(i1))
!        if(i1.gt.0) allocate(list_phon(i1))

        i1=0
        do i=min_p,max_p !1,id
         if(dabs(lhfp(i)%vi**2.d0-1.d0).lt.1.d-2) then
          do j=min_p,max_p !1,id
           if(dabs(lhfp(j)%ui**2.d0-1.d0).lt.1.d-2) then
            if(lhfp(i)%ipar*lhfp(j)%ipar.eq.ipar) then
             if((lhfp(i)%j2+lhfp(j)%j2)/2.ge.Jp) then
              if(abs(lhfp(i)%j2-lhfp(j)%j2)/2.le.Jp) then
               i1=i1+1
               ph(i1)%q1=i
               ph(i1)%q2=j
               ph(i1)%tz=-1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo

        do i=min_n,max_n !1,id
         if(dabs(lhfn(i)%vi**2.d0-1.d0).lt.1.d-2) then
          do j=min_n,max_n !1,id
           if(dabs(lhfn(j)%ui**2.d0-1.d0).lt.1.d-2) then
            if(lhfn(i)%ipar*lhfn(j)%ipar.eq.ipar) then
             if((lhfn(i)%j2+lhfn(j)%j2)/2.ge.Jp) then
              if(abs(lhfn(i)%j2-lhfn(j)%j2)/2.le.Jp) then
               i1=i1+1
               ph(i1)%q1=i
               ph(i1)%q2=j
               ph(i1)%tz=+1
              endif
             endif
            endif
           endif
          enddo
         endif
        enddo       
       endif

       if(if_QTDA.eq.1) then
        i1=0
        do i=min_p,max_p !1,id
         do j=min_p,max_p !1,id
          if(j.le.i) then
           if(lhfp(i)%ipar*lhfp(j)%ipar.eq.ipar) then
            if((lhfp(i)%j2+lhfp(j)%j2)/2.ge.Jp) then
             if(abs(lhfp(i)%j2-lhfp(j)%j2)/2.le.Jp) then
              if(.not.(i.eq.j.and.mod(Jp,2).eq.1)) then
               i1=i1+1
              endif
             endif
            endif
           endif
          endif
         enddo
        enddo

        do i=min_n,max_n !1,id
         do j=min_n,max_n !1,id
          if(j.le.i) then
           if(lhfn(i)%ipar*lhfn(j)%ipar.eq.ipar) then
            if((lhfn(i)%j2+lhfn(j)%j2)/2.ge.Jp) then
             if(abs(lhfn(i)%j2-lhfn(j)%j2)/2.le.Jp) then
              if(.not.(i.eq.j.and.mod(Jp,2).eq.1)) then
               i1=i1+1
              endif
             endif
            endif
           endif
          endif
         enddo
        enddo

        if(i1.gt.0) allocate(ph(i1))
!        if(i1.gt.0) allocate(list_phon(i1))

        i1=0
        do i=min_p,max_p !1,id
         do j=min_p,max_p !1,id
          if(j.le.i) then
           if(lhfp(i)%ipar*lhfp(j)%ipar.eq.ipar) then
            if((lhfp(i)%j2+lhfp(j)%j2)/2.ge.Jp) then
             if(abs(lhfp(i)%j2-lhfp(j)%j2)/2.le.Jp) then
              if(.not.(i.eq.j.and.mod(Jp,2).eq.1)) then
               i1=i1+1
               ph(i1)%q1=i
               ph(i1)%q2=j
               ph(i1)%tz=-1
              endif
             endif
            endif
           endif
          endif
         enddo
        enddo

        do i=min_n,max_n !1,id
         do j=min_n,max_n !1,id
          if(j.le.i) then
           if(lhfn(i)%ipar*lhfn(j)%ipar.eq.ipar) then
            if((lhfn(i)%j2+lhfn(j)%j2)/2.ge.Jp) then
             if(abs(lhfn(i)%j2-lhfn(j)%j2)/2.le.Jp) then
              if(.not.(i.eq.j.and.mod(Jp,2).eq.1)) then
               i1=i1+1
               ph(i1)%q1=i
               ph(i1)%q2=j
               ph(i1)%tz=+1
              endif
             endif
            endif
           endif
          endif
         enddo
        enddo
       endif

       i1_spur=0
       if(if_ort.eq.1.and.(ipar.eq.-1.and.Jp.eq.1)) i1_spur=1
       if(if_ort.eq.1.and.(if_QTDA.eq.1.and.(ipar.eq.1.and.Jp.eq.0))) 
     &                                                      i1_spur=1

       write(2,*) ' pi=',ipar,'J=',Jp
       write(2,*) 'dimension=',i1-i1_spur

       if(i1.gt.0) allocate(TDA_matrix(i1,i1))
       if(i1.gt.0) allocate(w(i1))
       if(i1.gt.0) TDA_matrix=0.d0

       if(if_ort.eq.1) then
       if(ipar.eq.-1.and.Jp.eq.1) allocate(spur_CM(i1))
       if(ipar.eq.-1.and.Jp.eq.1) 
     &  call spur_vec(if_QTDA,spur_CM,i1,Jp,ph)
       if((ipar.eq.1.and.Jp.eq.0).and.if_QTDA.eq.1)allocate(spur_N(i1))
       if((ipar.eq.1.and.Jp.eq.0).and.if_QTDA.eq.1) 
     &  call spur_vec(if_QTDA,spur_N,i1,Jp,ph) 
       endif

       
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
            Jpp_min=
     &min(iabs(lhfp(i)%j2-lhfp(k)%j2)/2,
     &iabs(lhfn(i)%j2-lhfn(k)%j2)/2,
     &iabs(lhfp(i)%j2-lhfn(k)%j2)/2)
     
            Jpp_max=
     &max((lhfp(i)%j2+lhfp(k)%j2)/2,
     &(lhfn(i)%j2+lhfn(k)%j2)/2,
     &(lhfp(i)%j2+lhfn(k)%j2)/2)

!            do Jpp=0,jmax
             do Jpp=Jpp_min,Jpp_max
             phasep=(-1)**((lhfp(j)%j2+lhfp(k)%j2)/2-Jp-Jpp)
             phasen=(-1)**((lhfn(j)%j2+lhfn(k)%j2)/2-Jp-Jpp)
             phasepn=(-1)**((lhfp(j)%j2+lhfn(k)%j2)/2-Jp-Jpp)
             trip=dble(2*Jpp+1)*sixj_int(lhfp(i)%j2,lhfp(j)%j2,2*Jp,
     &            lhfp(l)%j2,lhfp(k)%j2,2*Jpp)
             trin=dble(2*Jpp+1)*sixj_int(lhfn(i)%j2,lhfn(j)%j2,2*Jp,
     &            lhfn(l)%j2,lhfn(k)%j2,2*Jpp)
             tripn=dble(2*Jpp+1)*sixj_int(lhfp(i)%j2,lhfp(j)%j2,2*Jp,
     &            lhfn(l)%j2,lhfn(k)%j2,2*Jpp)
             if (trip.ne.0.0d0) valp=valp+phasep*trip*Vpp(i,k,j,l,Jpp)
             if (trin.ne.0.0d0) valn=valn+phasen*trin*Vnn(i,k,j,l,Jpp)
          if (tripn.ne.0.0d0) valpn=valpn+phasepn*tripn*Vpn(i,k,j,l,Jpp)
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

        if(if_QTDA.eq.0) then
         do i=1,i1
          p1=ph(i)%q2
          h1=ph(i)%q1
          tz1=ph(i)%tz
          do j=1,i1
           p2=ph(j)%q2
           h2=ph(j)%q1
           tz2=ph(j)%tz
           if(tz1.eq.-1.and.tz2.eq.-1) then
            phas=dble((-1)**(Jp+(lhfp(p1)%j2-lhfp(h1)%j2)/2))
            if(i.eq.j) TDA_matrix(i,j)=lhfp(p1)%ei-lhfp(h1)%ei
            TDA_matrix(i,j)=TDA_matrix(i,j)+Fpp(h1,p1,p2,h2)*phas
           endif
           if(tz1.eq.1.and.tz2.eq.1) then
            phas=dble((-1)**(Jp+(lhfn(p1)%j2-lhfn(h1)%j2)/2))
            if(i.eq.j) TDA_matrix(i,j)=lhfn(p1)%ei-lhfn(h1)%ei
            TDA_matrix(i,j)=TDA_matrix(i,j)+Fnn(h1,p1,p2,h2)*phas
           endif
           if(tz1.eq.-1.and.tz2.eq.1) then
            phas=dble((-1)**(Jp+(lhfp(p1)%j2-lhfp(h1)%j2)/2))
            TDA_matrix(i,j)=TDA_matrix(i,j)+Fpn(h1,p1,p2,h2)*phas
           endif
           if(tz1.eq.1.and.tz2.eq.-1) then
            phas=dble((-1)**(Jp+(lhfn(p1)%j2-lhfn(h1)%j2)/2))
            TDA_matrix(i,j)=TDA_matrix(i,j)+Fpn(p2,h2,h1,p1)*phas
           endif
          enddo
         enddo
        endif

        if(if_QTDA.eq.1) then
         do i=1,i1
          ia=ph(i)%q2
          ib=ph(i)%q1
          tz1=ph(i)%tz
          do j=1,i1
           ic=ph(j)%q2
           idd=ph(j)%q1
           tz2=ph(j)%tz
           if(tz1.eq.-1.and.tz2.eq.-1) then
            phas=dble((-1)**(-Jp+(lhfp(ic)%j2+lhfp(idd)%j2)/2))
            dnorm=1.d0
            if(ia.eq.ib) dnorm=dnorm/dsqrt(2.d0)
            if(ic.eq.idd) dnorm=dnorm/dsqrt(2.d0)
!            if(i.eq.j) TDA_matrix(i,j)=lhfp(ia)%qei+lhfp(ib)%qei
            if(ib.eq.idd) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                              +dnorm*H11p(ia,ic)
            if(ia.eq.ic) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                              +dnorm*H11p(ib,idd)
            if(ia.eq.idd) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                         -phas*dnorm*H11p(ib,ic)
            if(ib.eq.ic) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                         -phas*dnorm*H11p(ia,idd)
          TDA_matrix(i,j)=TDA_matrix(i,j)+Vpp(ia,ib,ic,idd,Jp)*dnorm
     &        *(lhfp(ia)%ui*lhfp(ib)%ui*lhfp(ic)%ui*lhfp(idd)%ui
     &            +lhfp(ia)%vi*lhfp(ib)%vi*lhfp(ic)%vi*lhfp(idd)%vi)
     &        +Fpp(ia,ib,ic,idd)*dnorm
     &        *(lhfp(ia)%ui*lhfp(ib)%vi*lhfp(ic)%vi*lhfp(idd)%ui
     &            +lhfp(ia)%vi*lhfp(ib)%ui*lhfp(ic)%ui*lhfp(idd)%vi)
     &        -phas*Fpp(ia,ib,idd,ic)*dnorm
     &        *(lhfp(ia)%ui*lhfp(ib)%vi*lhfp(ic)%ui*lhfp(idd)%vi
     &            +lhfp(ia)%vi*lhfp(ib)%ui*lhfp(ic)%vi*lhfp(idd)%ui)
           endif
           if(tz1.eq.1.and.tz2.eq.1) then
            phas=dble((-1)**(-Jp+(lhfn(ic)%j2+lhfn(idd)%j2)/2))
            dnorm=1.d0
            if(ia.eq.ib) dnorm=dnorm/dsqrt(2.d0)
            if(ic.eq.idd) dnorm=dnorm/dsqrt(2.d0)
!            if(i.eq.j) TDA_matrix(i,j)=lhfn(ia)%qei+lhfn(ib)%qei
            if(ib.eq.idd) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                              +dnorm*H11n(ia,ic)
            if(ia.eq.ic) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                              +dnorm*H11n(ib,idd)
            if(ia.eq.idd) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                         -phas*dnorm*H11n(ib,ic)
            if(ib.eq.ic) TDA_matrix(i,j)=TDA_matrix(i,j)
     &                                         -phas*dnorm*H11n(ia,idd)
          TDA_matrix(i,j)=TDA_matrix(i,j)+Vnn(ia,ib,ic,idd,Jp)*dnorm
     &        *(lhfn(ia)%ui*lhfn(ib)%ui*lhfn(ic)%ui*lhfn(idd)%ui
     &            +lhfn(ia)%vi*lhfn(ib)%vi*lhfn(ic)%vi*lhfn(idd)%vi)
     &        +Fnn(ia,ib,ic,idd)*dnorm
     &        *(lhfn(ia)%ui*lhfn(ib)%vi*lhfn(ic)%vi*lhfn(idd)%ui
     &            +lhfn(ia)%vi*lhfn(ib)%ui*lhfn(ic)%ui*lhfn(idd)%vi)
     &        -phas*Fnn(ia,ib,idd,ic)*dnorm
     &        *(lhfn(ia)%ui*lhfn(ib)%vi*lhfn(ic)%ui*lhfn(idd)%vi
     &            +lhfn(ia)%vi*lhfn(ib)%ui*lhfn(ic)%vi*lhfn(idd)%ui)
           endif
           if(tz1.eq.-1.and.tz2.eq.1) then
            phas=dble((-1)**(-Jp+(lhfp(ia)%j2+lhfp(ib)%j2)/2))
            dnorm=1.d0
            if(ia.eq.ib) dnorm=dnorm/dsqrt(2.d0)
            if(ic.eq.idd) dnorm=dnorm/dsqrt(2.d0)
            TDA_matrix(i,j)=TDA_matrix(i,j)+Fpn(ia,ib,ic,idd)*dnorm
     &        *(lhfp(ia)%vi*lhfp(ib)%ui*lhfn(ic)%ui*lhfn(idd)%vi
     &            +lhfp(ia)%ui*lhfp(ib)%vi*lhfn(ic)%vi*lhfn(idd)%ui)
     &        -phas*Fpn(ib,ia,ic,idd)*dnorm
     &        *(lhfp(ia)%ui*lhfp(ib)%vi*lhfn(ic)%ui*lhfn(idd)%vi
     &            +lhfp(ia)%vi*lhfp(ib)%ui*lhfn(ic)%vi*lhfn(idd)%ui)
           endif
           if(tz1.eq.1.and.tz2.eq.-1) then
            phas=dble((-1)**(-Jp+(lhfp(ic)%j2+lhfp(idd)%j2)/2))
            dnorm=1.d0
            if(ia.eq.ib) dnorm=dnorm/dsqrt(2.d0)
            if(ic.eq.id) dnorm=dnorm/dsqrt(2.d0)
            TDA_matrix(i,j)=TDA_matrix(i,j)+Fpn(ic,idd,ia,ib)*dnorm
     &        *(lhfp(ic)%vi*lhfp(idd)%ui*lhfn(ia)%ui*lhfn(ib)%vi
     &            +lhfp(ic)%ui*lhfp(idd)%vi*lhfn(ia)%vi*lhfn(ib)%ui)
     &        -phas*Fpn(idd,ic,ia,ib)*dnorm
     &        *(lhfp(ic)%ui*lhfp(idd)%vi*lhfn(ia)%ui*lhfn(ib)%vi
     &            +lhfp(ic)%vi*lhfp(idd)%ui*lhfn(ia)%vi*lhfn(ib)%ui)
           endif
          enddo
         enddo
        endif

        if(.not.(if_ort.eq.1.and.((ipar.eq.-1.and.Jp.eq.1).or.
     &        ((ipar.eq.1.and.Jp.eq.0).and.if_QTDA.eq.1)))) then


       write(91,*)'------------------------------------------------'
       do ii=1,i1
         write(91,'(1000f10.5)')(TDA_matrix(ii,jj),jj=1,i1)
       enddo


        call diagonalization(TDA_matrix,w,i1)

         write(1,*) ' pi=',ipar,'J=',Jp
         write(1,*) w(1:i1) 

        call transit_calc(ipar,Jp,ph,TDA_matrix,w,i1)
        call phonon_density_calc(ipar,Jp,ph,TDA_matrix,w,i1)

        write(3,*) ' pi=',ipar,'J=',Jp
        write(3,*)
        do i=1,i1
         write(3,*) 'phonon  E_i=',w(i)
         write(3,*) 
         do j=1,i1
          m1=ph(j)%q1
          m2=ph(j)%q2
          m3=ph(j)%tz
          if(TDA_matrix(j,i)**2.d0.gt.1.d-2) then
           if(m3.eq.-1) then
            write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))') 
     &             'p',lhfp(m1)%l,lhfp(m1)%j2,lhfp(m1)%ei,lhfp(m2)%l,
     &                 lhfp(m2)%j2,lhfp(m2)%ei,TDA_matrix(j,i)**2.d0
           endif
           if(m3.eq.+1) then
            write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))') 
     &             'n',lhfn(m1)%l,lhfn(m1)%j2,lhfn(m1)%ei,lhfn(m2)%l,
     &                 lhfn(m2)%j2,lhfn(m2)%ei,TDA_matrix(j,i)**2.d0
           endif
          endif
         enddo
         write(3,*) 
        enddo

        elseif(ipar.eq.-1.and.Jp.eq.1) then
         call diag_subspace(TDA_matrix,w,spur_CM,i1)
         write(1,*) ' pi=',ipar,'J=',Jp
         write(1,*) w(2:i1)
         call transit_calc_spur(ipar,Jp,ph,TDA_matrix,w,i1)
         call phonon_density_calc(ipar,Jp,ph,TDA_matrix,w,i1)
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
     &             'p',lhfp(m1)%l,lhfp(m1)%j2,lhfp(m1)%ei,lhfp(m2)%l,
     &                 lhfp(m2)%j2,lhfp(m2)%ei,TDA_matrix(j,i)**2.d0
            endif
            if(m3.eq.+1) then
             write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
     &             'n',lhfn(m1)%l,lhfn(m1)%j2,lhfn(m1)%ei,lhfn(m2)%l,
     &                 lhfn(m2)%j2,lhfn(m2)%ei,TDA_matrix(j,i)**2.d0
            endif
           endif
          enddo
          write(3,*)
         enddo
        else
         call diag_subspace(TDA_matrix,w,spur_N,i1)
         write(1,*) ' pi=',ipar,'J=',Jp
         write(1,*) w(2:i1)
         call transit_calc_spur(ipar,Jp,ph,TDA_matrix,w,i1)
         call phonon_density_calc(ipar,Jp,ph,TDA_matrix,w,i1)
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
     &             'p',lhfp(m1)%l,lhfp(m1)%j2,lhfp(m1)%ei,lhfp(m2)%l,
     &                 lhfp(m2)%j2,lhfp(m2)%ei,TDA_matrix(j,i)**2.d0
            endif
            if(m3.eq.+1) then
             write(3,'((a,1x),2(i2,1x),1(f9.3,1x),2(i2,1x),2(f9.3))')
     &             'n',lhfn(m1)%l,lhfn(m1)%j2,lhfn(m1)%ei,lhfn(m2)%l,
     &                 lhfn(m2)%j2,lhfn(m2)%ei,TDA_matrix(j,i)**2.d0
            endif
           endif
          enddo
          write(3,*)
         enddo
        endif

       write(4) ipar,Jp,i1
       write(4) (ph(i),i=1,i1)
       write(4) (w(i),i=1,i1)

       do i=1,i1
!        list_phon(i)%noscmax=nosc_TDA
!        list_phon(i)%ipar=ipar
!        list_phon(i)%jproj=Jp
!        list_phon(i)%dimp=id
!        list_phon(i)%dimn=id
!        list_phon(i)%ener=w(i)
!        list_phon(i)%cph_p=0.d0
!        list_phon(i)%cph_n=0.d0
!        do k1=1,i1
!         m1=ph(k1)%q1
!         m2=ph(k1)%q2
!         m3=ph(k1)%tz
!         if(m3.eq.-1) list_phon(i)%cph_p(m1,m2)=TDA_matrix(k1,i)
!         if(m3.eq.+1) list_phon(i)%cph_n(m1,m2)=TDA_matrix(k1,i)
!        enddo
        write(4) (TDA_matrix(j,i),j=1,i1)  !list_phon(i)
       enddo

       endif

        if(i1.gt.0) deallocate(ph,TDA_matrix,w)  !,list_phon)

        if(if_ort.eq.1) then
        if(ipar.eq.-1.and.Jp.eq.1) deallocate(spur_CM)
        if((ipar.eq.1.and.Jp.eq.0).and.if_QTDA.eq.1)deallocate(spur_N)
        endif

        enddo
       enddo

!       deallocate(sixj1)
       deallocate(Fpp,Fnn,Fpn)


       close(1)
       close(2)
       close(3)
       close(4)
       write(*,*) 'TDA calculation has finished.'
        

       return
      end
