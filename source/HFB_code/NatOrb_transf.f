      Subroutine NatOrb_transf

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpn_HFB(:,:,:,:,:)

       double precision :: Tp(id,id),Tn(id,id)
       double precision :: bp(id,id),bn(id,id)

       Tp=tran_p
       Tn=tran_n

       allocate(Vpp_HFB(id,id,id,id,0:jmax))
       allocate(Vnn_HFB(id,id,id,id,0:jmax))
       allocate(Vpn_HFB(id,id,id,id,0:jmax))
       Vpp_HFB=0.d0
       Vnn_HFB=0.d0
       Vpn_HFB=0.d0

       write(*,*) 'Transf. interaction in 1st index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            valpY=0.d0
            valnY=0.d0
            do i2=1,id
             valp=valp+Vpp(i2,j1,k1,l1,Jp)*Tp(i2,i1)
             valn=valn+Vnn(i2,j1,k1,l1,Jp)*Tn(i2,i1)
             val0=val0+Vpn(i2,j1,k1,l1,Jp)*Tp(i2,i1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB

       write(*,*) 'Transf. interaction in 2nd index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            do j2=1,id
             valp=valp+Vpp(i1,j2,k1,l1,Jp)*Tp(j2,j1)
             valn=valn+Vnn(i1,j2,k1,l1,Jp)*Tn(j2,j1)
             val0=val0+Vpn(i1,j2,k1,l1,Jp)*Tn(j2,j1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB

       write(*,*) 'Transf. interaction in 3rd index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            do k2=1,id
             valp=valp+Vpp(i1,j1,k2,l1,Jp)*Tp(k2,k1)
             valn=valn+Vnn(i1,j1,k2,l1,Jp)*Tn(k2,k1)
             val0=val0+Vpn(i1,j1,k2,l1,Jp)*Tp(k2,k1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB

       write(*,*) 'Transf. interaction in 4th index'
       do Jp=0,jmax
        do i1=1,id
         do j1=1,id
          do k1=1,id
           do l1=1,id
            valp=0.d0
            valn=0.d0
            val0=0.d0
            valpY=0.d0
            valnY=0.d0
            do l2=1,id
             valp=valp+Vpp(i1,j1,k1,l2,Jp)*Tp(l2,l1)
             valn=valn+Vnn(i1,j1,k1,l2,Jp)*Tn(l2,l1)
             val0=val0+Vpn(i1,j1,k1,l2,Jp)*Tn(l2,l1)
            enddo
            Vpp_HFB(i1,j1,k1,l1,Jp)=valp
            Vnn_HFB(i1,j1,k1,l1,Jp)=valn
            Vpn_HFB(i1,j1,k1,l1,Jp)=val0
           enddo
          enddo
         enddo
        enddo
       enddo
       Vpp=Vpp_HFB
       Vnn=Vnn_HFB
       Vpn=Vpn_HFB

       deallocate(Vpp_HFB,Vnn_HFB,Vpn_HFB)

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+kin_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+kin_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       kin_p=bp
       kin_n=bn

 !  unformatted output of kinetic energy operators in natural orbital basis      
       open(165,file='kin_nat_orb.dat',status='unknown'
     &,form='unformatted')
       write(165)id
       write(165)((kin_p(i,j),i=1,id),j=1,id)
       write(165)((kin_n(i,j),i=1,id),j=1,id)
       close(165)

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE0_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE0_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE0_p=bp
       trE0_n=bn

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE1_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE1_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE1_p=bp
       trE1_n=bn

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE2_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE2_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE2_p=bp
       trE2_n=bn

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trE3_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trE3_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trE3_p=bp
       trE3_n=bn

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trEN_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trEN_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trEN_p=bp
       trEN_n=bn

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trS1_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trS1_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trS1_p=bp
       trS1_n=bn

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trM1s_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trM1s_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trM1s_p=bp
       trM1s_n=bn

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+trM1l_p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+trM1l_n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       trM1l_p=bp
       trM1l_n=bn

       open(165,file='vlk_nat_orb.dat',status='unknown'
     &,form='unformatted')
       open(1,file='vlk_NO.dat',status='unknown',form='formatted')
       write(1,*) '  Tz Par  2J   a   b   c   d          <ab|V|cd> '
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1-1
            j=2*i2-1
            k=2*i3-1
            l=2*i4-1
            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              factab=1.d0
              factcd=1.d0
              if(i.eq.j) factab=dsqrt(2.d0)
              if(k.eq.l) factcd=dsqrt(2.d0)
              xnorm=factab*factcd
              if(dabs(Vpp(i1,i2,i3,i4,Jp)).gt.precis) then 
              write(1,*)
     &        -1,0,2*Jp,i,j,k,l,Vpp(i1,i2,i3,i4,Jp)/xnorm
              write(165)int(-1,1),int(2*Jp,1),
     & int(i,2),int(j,2),int(k,2),int(l,2),(Vpp(i1,i2,i3,i4,Jp)/xnorm)
              endif 
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1-1
            j=2*i2
            k=2*i3-1
            l=2*i4
!            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              if(dabs(Vpn(i1,i2,i3,i4,Jp)).gt.precis) then 
              write(1,*) 0,0,2*Jp,i,j,k,l,Vpn(i1,i2,i3,i4,Jp)
              write(165)int(0,1),int(2*Jp,1),
     &  int(i,2),int(j,2),int(k,2),int(l,2),Vpn(i1,i2,i3,i4,Jp)
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       do Jp=0,jmax
        do i1=1,id
         do i2=1,id
          do i3=1,id
           do i4=1,id
            i=2*i1
            j=2*i2
            k=2*i3
            l=2*i4
            if(i.le.j.and.k.le.l) then
             if(1000*i+j.le.1000*k+l) then
              factab=1.d0
              factcd=1.d0
              if(i.eq.j) factab=dsqrt(2.d0)
              if(k.eq.l) factcd=dsqrt(2.d0)
              xnorm=factab*factcd
              if(dabs(Vnn(i1,i2,i3,i4,Jp)).gt.precis) then 
              write(1,*)
     &        1,0,2*Jp,i,j,k,l,Vnn(i1,i2,i3,i4,Jp)/xnorm
              write(165)int(1,1),int(2*Jp,1),
     & int(i,2),int(j,2),int(k,2),int(l,2),(Vnn(i1,i2,i3,i4,Jp)/xnorm)
             endif
            endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       close(1)
       close(165)

       return
      end
