      Subroutine NatOrb

       USE technical
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: hpp(id,id)
       double precision :: hnn(id,id)
       double precision, dimension(:,:),allocatable :: h_temp
       double precision :: hll(id,id)

       double precision :: Op(id,id),On(id,id),OY(id,id)
       double precision :: wnop(id),wnon(id)
       double precision, dimension(:),allocatable :: w_temp

       integer :: jlp(id),jln(id)

       hpp=0.d0
       hnn=0.d0
       hll=0.d0

       Op=0.d0
       On=0.d0
       OY=0.d0

       wnop=0.d0
       wnon=0.d0

!       valp=0.d0

!       do i=1,id
!        if(dabs(lhfp(i)%vi).gt.0.5d0) then
!        do j=1,id
!         if(dabs(lhfp(j)%vi).gt.0.5d0) then
!         do k=1,id
!          if(dabs(lhfp(k)%vi).lt.0.5d0) then
!          do l=1,id
!           if(dabs(lhfp(l)%vi).lt.0.5d0) then
!           do Jp=0,jmax
!           valp=valp+0.25d0*Vpp(k,l,i,j,Jp)**2.d0/(lhfp(i)%ei+lhfp(j)%ei
!     &                                   -lhfp(k)%ei-lhfp(l)%ei)**2.d0
!           enddo
!           endif
!          enddo
!          endif
!         enddo
!         endif
!        enddo
!        endif
!       enddo
!       do i=1,id
!        if(dabs(lhfn(i)%vi).gt.0.5d0) then
!        do j=1,id
!         if(dabs(lhfn(j)%vi).gt.0.5d0) then
!         do k=1,id
!          if(dabs(lhfn(k)%vi).lt.0.5d0) then
!          do l=1,id
!           if(dabs(lhfn(l)%vi).lt.0.5d0) then
!           do Jp=0,jmax
!           valp=valp+0.25d0*Vnn(k,l,i,j,Jp)**2.d0/(lhfn(i)%ei+lhfn(j)%ei
!     &                                   -lhfn(k)%ei-lhfn(l)%ei)**2.d0
!           enddo
!           endif
!          enddo
!          endif
!         enddo
!         endif
!        enddo
!        endif
!       enddo
!       do i=1,id
!        if(dabs(lhfp(i)%vi).gt.0.5d0) then
!        do j=1,id
!         if(dabs(lhfn(j)%vi).gt.0.5d0) then
!         do k=1,id
!          if(dabs(lhfp(k)%vi).lt.0.5d0) then
!          do l=1,id
!           if(dabs(lhfn(l)%vi).lt.0.5d0) then
!           do Jp=0,jmax
!           valp=valp+Vpn(k,l,i,j,Jp)**2.d0/(lhfp(i)%ei+lhfn(j)%ei
!     &                                   -lhfp(k)%ei-lhfn(l)%ei)**2.d0
!           enddo
!           endif
!          enddo
!          endif
!         enddo
!         endif
!        enddo
!        endif
!       enddo

       do j=1,id
        do i=j,id
         if((lhfp(i)%j2.eq.lhfp(j)%j2).and.(lhfp(i)%l.eq.lhfp(j)%l)) 
     &                                                          then
          if(dabs(lhfp(i)%vi).gt.0.5d0.and.dabs(lhfp(j)%vi).gt.0.5d0)
     &                                                           then

           if(i.eq.j) then
            hpp(j,i)=hpp(j,i)+1.d0!+valp
           endif

           do ip1=1,id
            if(dabs(lhfp(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)-0.5d0*dble(2*Jp+1)
     &  *Vpp(ip1,ip2,i,ih2,Jp)*Vpp(ip1,ip2,j,ih2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(i)%ei+lhfp(ih2)%ei-lhfp(ip1)%ei
     &   -lhfp(ip2)%ei)*(lhfp(j)%ei+lhfp(ih2)%ei-lhfp(ip1)%ei
     &                                           -lhfp(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ip1=1,id
            if(dabs(lhfp(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)-dble(2*Jp+1)
     &  *Vpn(ip1,ip2,i,ih2,Jp)*Vpn(ip1,ip2,j,ih2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(i)%ei+lhfn(ih2)%ei-lhfp(ip1)%ei
     &   -lhfn(ip2)%ei)*(lhfp(j)%ei+lhfn(ih2)%ei-lhfp(ip1)%ei
     &                                           -lhfn(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo
          endif

          if(dabs(lhfp(i)%vi).lt.0.5d0.and.dabs(lhfp(j)%vi).lt.0.5d0)
     &                                                           then

           do ih1=1,id
            if(dabs(lhfp(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)+0.5d0*dble(2*Jp+1)
     &  *Vpp(j,ip2,ih1,ih2,Jp)*Vpp(i,ip2,ih1,ih2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(ih1)%ei+lhfp(ih2)%ei-lhfp(j)%ei
     &   -lhfp(ip2)%ei)*(lhfp(ih1)%ei+lhfp(ih2)%ei-lhfp(i)%ei
     &                                           -lhfp(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ih1=1,id
            if(dabs(lhfp(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)+dble(2*Jp+1)
     &  *Vpn(j,ip2,ih1,ih2,Jp)*Vpn(i,ip2,ih1,ih2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(ih1)%ei+lhfn(ih2)%ei-lhfp(j)%ei
     &   -lhfn(ip2)%ei)*(lhfp(ih1)%ei+lhfn(ih2)%ei-lhfp(i)%ei
     &                                           -lhfn(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

          endif

          if(dabs(lhfp(i)%vi).lt.0.5d0.and.dabs(lhfp(j)%vi).gt.0.5d0)
     &                                                           then

           do ip1=1,id
            if(dabs(lhfp(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)+0.5d0*dble(2*Jp+1)
     &  *Vpp(i,ih2,ip1,ip2,Jp)*Vpp(ip1,ip2,j,ih2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(j)%ei+lhfp(ih2)%ei-lhfp(ip1)%ei
     &   -lhfp(ip2)%ei)*(lhfp(j)%ei-lhfp(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ih1=1,id
            if(dabs(lhfp(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)-0.5d0*dble(2*Jp+1)
     &  *Vpp(i,ip2,ih1,ih2,Jp)*Vpp(ih1,ih2,j,ip2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(ih1)%ei+lhfp(ih2)%ei-lhfp(i)%ei
     &   -lhfp(ip2)%ei)*(lhfp(j)%ei-lhfp(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ip1=1,id
            if(dabs(lhfp(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)+dble(2*Jp+1)
     &  *Vpn(i,ih2,ip1,ip2,Jp)*Vpn(ip1,ip2,j,ih2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(j)%ei+lhfn(ih2)%ei-lhfp(ip1)%ei
     &   -lhfn(ip2)%ei)*(lhfp(j)%ei-lhfp(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ih1=1,id
            if(dabs(lhfp(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hpp(j,i)=hpp(j,i)-dble(2*Jp+1)
     &  *Vpn(i,ip2,ih1,ih2,Jp)*Vpn(ih1,ih2,j,ip2,Jp)
     &  /(dble(lhfp(i)%j2+1)*(lhfp(ih1)%ei+lhfn(ih2)%ei-lhfp(i)%ei
     &   -lhfn(ip2)%ei)*(lhfp(j)%ei-lhfp(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

          endif
         endif
        enddo
       enddo

       do j=1,id
        do i=j+1,id
         hpp(i,j)=hpp(j,i)
        enddo
       enddo

       do j=1,id
        do i=j,id
         if((lhfn(i)%j2.eq.lhfn(j)%j2).and.(lhfn(i)%l.eq.lhfn(j)%l)) 
     &                                                          then
          if(dabs(lhfn(i)%vi).gt.0.5d0.and.dabs(lhfn(j)%vi).gt.0.5d0)
     &                                                           then

           if(i.eq.j) then
            hnn(j,i)=hnn(j,i)+1.d0!+valp
           endif

           do ip1=1,id
            if(dabs(lhfn(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)-0.5d0*dble(2*Jp+1)
     &  *Vnn(ip1,ip2,i,ih2,Jp)*Vnn(ip1,ip2,j,ih2,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(i)%ei+lhfn(ih2)%ei-lhfn(ip1)%ei
     &   -lhfn(ip2)%ei)*(lhfn(j)%ei+lhfn(ih2)%ei-lhfn(ip1)%ei
     &                                           -lhfn(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ip1=1,id
            if(dabs(lhfn(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)-dble(2*Jp+1)
     &  *Vpn(ip2,ip1,ih2,i,Jp)*Vpn(ip2,ip1,ih2,j,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(i)%ei+lhfp(ih2)%ei-lhfn(ip1)%ei
     &   -lhfp(ip2)%ei)*(lhfn(j)%ei+lhfp(ih2)%ei-lhfn(ip1)%ei
     &                                           -lhfp(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo
          endif

          if(dabs(lhfn(i)%vi).lt.0.5d0.and.dabs(lhfn(j)%vi).lt.0.5d0)
     &                                                           then

           do ih1=1,id
            if(dabs(lhfn(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)+0.5d0*dble(2*Jp+1)
     &  *Vnn(j,ip2,ih1,ih2,Jp)*Vnn(i,ip2,ih1,ih2,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(ih1)%ei+lhfn(ih2)%ei-lhfn(j)%ei
     &   -lhfn(ip2)%ei)*(lhfn(ih1)%ei+lhfn(ih2)%ei-lhfn(i)%ei
     &                                           -lhfn(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ih1=1,id
            if(dabs(lhfn(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)+dble(2*Jp+1)
     &  *Vpn(ip2,j,ih2,ih1,Jp)*Vpn(ip2,i,ih2,ih1,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(ih1)%ei+lhfp(ih2)%ei-lhfn(j)%ei
     &   -lhfp(ip2)%ei)*(lhfn(ih1)%ei+lhfp(ih2)%ei-lhfn(i)%ei
     &                                           -lhfp(ip2)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

          endif

          if(dabs(lhfn(i)%vi).lt.0.5d0.and.dabs(lhfn(j)%vi).gt.0.5d0)
     &                                                           then

           do ip1=1,id
            if(dabs(lhfn(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)+0.5d0*dble(2*Jp+1)
     &  *Vnn(i,ih2,ip1,ip2,Jp)*Vnn(ip1,ip2,j,ih2,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(j)%ei+lhfn(ih2)%ei-lhfn(ip1)%ei
     &   -lhfn(ip2)%ei)*(lhfn(j)%ei-lhfn(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ih1=1,id
            if(dabs(lhfn(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfn(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfn(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)-0.5d0*dble(2*Jp+1)
     &  *Vnn(i,ip2,ih1,ih2,Jp)*Vnn(ih1,ih2,j,ip2,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(ih1)%ei+lhfn(ih2)%ei-lhfn(i)%ei
     &   -lhfn(ip2)%ei)*(lhfn(j)%ei-lhfn(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ip1=1,id
            if(dabs(lhfn(ip1)%vi).lt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)+dble(2*Jp+1)
     &  *Vpn(ih2,i,ip2,ip1,Jp)*Vpn(ip2,ip1,ih2,j,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(j)%ei+lhfp(ih2)%ei-lhfn(ip1)%ei
     &   -lhfp(ip2)%ei)*(lhfn(j)%ei-lhfn(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

           do ih1=1,id
            if(dabs(lhfn(ih1)%vi).gt.0.5d0) then
            do ip2=1,id
             if(dabs(lhfp(ip2)%vi).lt.0.5d0) then
             do ih2=1,id
              if(dabs(lhfp(ih2)%vi).gt.0.5d0) then
               do Jp=0,jmax
                hnn(j,i)=hnn(j,i)-dble(2*Jp+1)
     &  *Vpn(ip2,i,ih2,ih1,Jp)*Vpn(ih2,ih1,ip2,j,Jp)
     &  /(dble(lhfn(i)%j2+1)*(lhfn(ih1)%ei+lhfp(ih2)%ei-lhfn(i)%ei
     &   -lhfp(ip2)%ei)*(lhfn(j)%ei-lhfn(i)%ei))
               enddo
              endif
             enddo
             endif
            enddo
            endif
           enddo

          endif
         endif
        enddo
       enddo

       do j=1,id
        do i=j+1,id
         hnn(i,j)=hnn(j,i)
        enddo
       enddo

!       do i=1,id
!        do j=1,id
!         if(dabs(hpp(i,j)).gt.1.d-7) write(401,*) i,j,hpp(i,j)
!         if(dabs(hnn(i,j)).gt.1.d-7) write(402,*) i,j,hnn(i,j)
!        enddo
!       enddo
!       sump=0.d0
!       sumn=0.d0
!       do i=1,id
!        sump=sump+hpp(i,i)*dble(lhfp(i)%j2+1)
!        sumn=sumn+hnn(i,i)*dble(lhfn(i)%j2+1)
!       enddo
!       write(403,*) sump,sumn

       call diagonalization(hpp,wnop,id)
       call diagonalization(hnn,wnon,id)

!   reordering of eigenvalues in descending order 
       allocate(w_temp(id))
       allocate(h_temp(id,id))

       do i=1,id
          w_temp(i)=wnop(id-i+1)
          do j=1,id
            h_temp(j,i)=hpp(j,id-i+1)
          enddo
       enddo

       wnop=w_temp
       hpp=h_temp

       do i=1,id
          w_temp(i)=wnon(id-i+1)
          do j=1,id
           h_temp(j,i)=hnn(j,id-i+1)
          enddo
       enddo

       wnon=w_temp
       hnn=h_temp

       
       
       tran_p=hpp
       tran_n=hnn

!       write(501,*) wnop
!       write(502,*) wnon

!       do i=1,id
!        do j=1,id
!         pij=0.d0
!         do k=1,id
!          pij=pij+hpp(k,i)*hpp(k,j)
!         enddo
!         write(503,*) i,j,pij
!        enddo
!       enddo

!       do i=1,id
!        do j=1,id
!         pij=0.d0
!         do k=1,id
!          pij=pij+hnn(k,i)*hnn(k,j)
!         enddo
!         write(504,*) i,j,pij
!        enddo
!       enddo

       if (.not.allocated(lnop)) allocate(lnop(id),lnon(id))


       do i=1,id
        d_max=0.d0
        do j=1,id
         if(dabs(hpp(j,i)).gt.d_max) then
          d_max=dabs(hpp(j,i))
          j_pointer=j
         endif
        enddo
        jlp(i)=1000*lhfp(j_pointer)%l+lhfp(j_pointer)%j2
       enddo
       do i=1,id
        llll=jlp(i)/1000
        jjjj=mod(jlp(i),1000)
        lnop(i)%ipar=(-1)**llll
        lnop(i)%l=llll
        lnop(i)%j2=jjjj
        lnop(i)%ei=wnop(i)
       enddo

       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=0
         do i=1,id
          if(lhfp(i)%j2.eq.jj.and.lhfp(i)%l.eq.ll) then
           nn=nn+1
           mm=0
           do j=1,id!id,1,-1
            if(lnop(j)%j2.eq.jj.and.lnop(j)%l.eq.ll) then
             mm=mm+1
             if(mm.eq.nn) lnop(j)%ui=lhfp(i)%ui
             if(mm.eq.nn) lnop(j)%vi=lhfp(i)%vi
            endif
           enddo
          endif
         enddo
        enddo
       enddo

       do i=1,id
        d_max=0.d0
        do j=1,id
         if(dabs(hnn(j,i)).gt.d_max) then
          d_max=dabs(hnn(j,i))
          j_pointer=j
         endif
        enddo
        jlp(i)=1000*lhfn(j_pointer)%l+lhfn(j_pointer)%j2
       enddo
       do i=1,id
        llll=jlp(i)/1000
        jjjj=mod(jlp(i),1000)
        lnon(i)%ipar=(-1)**llll
        lnon(i)%l=llll
        lnon(i)%j2=jjjj
        lnon(i)%ei=wnon(i)
       enddo

       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=0
         do i=1,id
          if(lhfn(i)%j2.eq.jj.and.lhfn(i)%l.eq.ll) then
           nn=nn+1
           mm=0
           do j=1,id!,1,-1
            if(lnon(j)%j2.eq.jj.and.lnon(j)%l.eq.ll) then
             mm=mm+1
             if(mm.eq.nn) lnon(j)%ui=lhfn(i)%ui
             if(mm.eq.nn) lnon(j)%vi=lhfn(i)%vi
            endif
           enddo
          endif
         enddo
        enddo
       enddo

       open(1,file='NAT_p.out',status='unknown',form='formatted')
        write(1,*) 'i,     l,    2*j,      ei,     vi'
        do ii = 1,id!,1,-1
         write(1,*) ii,lnop(ii)%l,lnop(ii)%j2,
     &    lnop(ii)%ei,lnop(ii)%vi
        end do
       close(1)

       open(1,file='NAT_n.out',status='unknown',form='formatted')
        write(1,*) 'i,     l,    2*j,      ei,     vi'
        do ii = 1,id!,1,-1
         write(1,*) ii,lnon(ii)%l,lnon(ii)%j2,
     &    lnon(ii)%ei,lnon(ii)%vi
        end do
       close(1)

   

       open(4,file='unitar_NATp.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (hpp(i1,j1),j1=1,id) 
        enddo
       close(4)

       open(4,file='unitar_NATn.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (hnn(i1,j1),j1=1,id) 
        enddo
       close(4)

       return
      end
