       subroutine transf_interaction
       
       USE technical
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn_HFB(:,:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpn_HFB(:,:,:,:,:)

       double precision :: Tp(id,id),Tn(id,id)
       double precision :: bp(id,id),bn(id,id)
       double precision, allocatable :: vv(:,:),uu(:,:),vt(:,:) 

       integer :: point_p(id),point_n(id)
       integer :: pinv_p(id),pinv_n(id)

!       Vpp=Vpp+3.d0*Vpp_DD
!       Vnn=Vnn+3.d0*Vnn_DD
!       Vpn=Vpn+3.d0*Vpn_DD

       Tp=tran_p
       Tn=tran_n


       allocate(tran_HFp(id,id),tran_HFn(id,id))

       tran_HFp = tran_p
       tran_HFn = tran_n

       allocate(tran_NOp(id,id),tran_NOn(id,id))

       tran_NOp = 0.d0
       tran_NOn = 0.d0
       
       
!  toto tu doplnujem ja
 !      do i=1,id
       
 !       do j=1,id
 !        if (dabs(Tp(j,i)).gt.0.01) jlarge=j
 !       enddo
        
 !       if (Tp(jlarge,i).lt.0.0d0) then
 !        do j=1,id
 !         Tp(j,i)=-1.0d0*Tp(j,i)
 !        enddo
 !       endif 
       
 !      enddo

 !      do i=1,id
       
 !       do j=1,id
 !        if (dabs(Tn(j,i)).gt.0.01) jlarge=j
 !       enddo
       
 !       if (Tn(jlarge,i).lt.0.0d0) then
 !        do j=1,id
 !         Tn(j,i)=-1.0d0*Tn(j,i)
 !        enddo
 !       endif 
       
 !      enddo

!     transformation of kinetic energy ops


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

       open(61,file='kin_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='kin_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


       

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

       open(61,file='r2Y0_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r2Y0_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)





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

       open(61,file='r1Y1_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r1Y1_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


!************************
!       open(61,file='E1p.out',status='unknown',form='formatted')
!        do i=1,id
!         write(61,*) (trE1_p(i,j),j=1,id)
!        enddo
!       close(61)
!       open(62,file='E1n.out',status='unknown',form='formatted')
!        do i=1,id
!         write(62,*) (trE1_n(i,j),j=1,id)
!        enddo
!       close(62)
!************************

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

       open(61,file='r2Y2_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r2Y2_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)



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

       open(61,file='r3Y3_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r3Y3_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


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

       open(61,file='EN_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='EN_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)



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

       open(61,file='r3Y1_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r3Y1_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)



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

       open(61,file='M1s_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='M1s_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)



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

       open(61,file='M1l_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='M1l_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)


       do m=0,igrid2
        bp=0.d0
        bn=0.d0
        do i=1,id
         do j=1,id
          val1=0.d0
          val2=0.d0
          do k=1,id
           do l=1,id
            val1=val1+trE1_p_dens(k,l,m)*Tp(k,i)*Tp(l,j)
            val2=val2+trE1_n_dens(k,l,m)*Tn(k,i)*Tn(l,j)
           enddo
          enddo
          bp(i,j)=val1
          bn(i,j)=val2
         enddo
        enddo
        do i=1,id
         do j=1,id
          trE1_p_dens(i,j,m)=bp(i,j)
          trE1_n_dens(i,j,m)=bn(i,j)
         enddo
        enddo
       enddo

       bp=0.d0
       bn=0.d0
       do i=1,id
        do j=1,id
         val1=0.d0
         val2=0.d0
         do k=1,id
          do l=1,id
           val1=val1+H11p(k,l)*Tp(k,i)*Tp(l,j)
           val2=val2+H11n(k,l)*Tn(k,i)*Tn(l,j)
          enddo
         enddo
         bp(i,j)=val1
         bn(i,j)=val2
        enddo
       enddo
       H11p=bp
       H11n=bn

       do i=1,id
        do j=1,id
         H11p(i,j)=-(lhfp(i)%ui*lhfp(j)%vi+lhfp(i)%vi*lhfp(j)%ui)*H11p(i,j)
         H11n(i,j)=-(lhfn(i)%ui*lhfn(j)%vi+lhfn(i)%vi*lhfn(j)%ui)*H11n(i,j)
        enddo
       enddo
       do i=1,id
        H11p(i,i)=H11p(i,i)+(lhfp(i)%ui**2.d0-lhfp(i)%vi**2.d0)*(lhfp(i)%ei-ferp)
        H11n(i,i)=H11n(i,i)+(lhfn(i)%ui**2.d0-lhfn(i)%vi**2.d0)*(lhfn(i)%ei-fern)
       enddo

       open(4,file='H11p.dat',status='unknown',form='formatted')
       write(4,*) ferp
        do i1=1,id
         write(4,*) (H11p(i1,j1),j1=1,id) 
        enddo
       close(4)
       open(4,file='H11n.dat',status='unknown',form='formatted')
       write(4,*) fern
        do i1=1,id
         write(4,*) (H11n(i1,j1),j1=1,id) 
        enddo
       close(4)
       
       
       open(1,file='vlk_hfb.dat',status='unknown',form='unformatted')
       
       allocate(vv(id*id,id*id),uu(id*id,id*id),vt(id*id,id*id))
       
!   Proton-proton part       
       vv=0.0d0
       uu=0.0d0
       vt=0.0d0
       
       do i=1,id
         do j=1,id 
           ii=id*(i-1)+j    
          do k=1,id
            do l=1,id
             jj=id*(k-1)+l
              uu(jj,ii)=Tp(k,i)*Tp(l,j)
            enddo
          enddo
         enddo
       enddo
       
       do Jp=0,jmax
       write(*,*)'Transformation of PP interaction for J =',Jp
        vv=0.d0
        do i=1,id 
         do j=1,id
          ii=id*(i-1)+j
          do k=1,id
           do l=1,id
             jj=id*(k-1)+l
              if (if_dd.eq.1) then 
              vv(ii,jj)=Vpp(i,j,k,l,Jp)+3.0d0*Vpp_dd(i,j,k,l,Jp)
              else
              vv(ii,jj)=Vpp(i,j,k,l,Jp)
              endif 
           enddo
          enddo
         enddo
        enddo
        
       call dgemm('N','N',id*id,id*id,id*id,1.d0,vv,id*id,uu,id*id,0.d0,vt,id*id)
       vv=0.d0
       call dgemm('T','N',id*id,id*id,id*id,1.d0,uu,id*id,vt,id*id,0.d0,vv,id*id)
       
      
       do ii=1,id*id
          i=int(ii/id)+1          
          j=mod(ii,id) 
              
          if (mod(ii,id).eq.0) then 
           i=ii/id
           j=id
          endif
        
 
        do jj=1,id*id
           k=int(jj/id)+1
           l=mod(jj,id) 
           
          if (mod(jj,id).eq.0) then 
           k=jj/id
           l=id
          endif
 
         if((i.le.j.and.k.le.l).and.(100000*i+j.le.100000*k+l)) then
             
         if(i.eq.j) vv(ii,jj)=vv(ii,jj)/dsqrt(2.d0)
         if(k.eq.l) vv(ii,jj)=vv(ii,jj)/dsqrt(2.d0)
 !        if (dabs(vv(ii,jj)).gt.1.d-10) write(1000,*)-1,2*Jp,2*i-1,2*j-1,2*k-1,2*l-1,vv(ii,jj)
         if (dabs(vv(ii,jj)).gt.1.d-10) write(1)int(-1,1),int(2*Jp,1),int(2*i-1,2),int(2*j-1,2),int(2*k-1,2),int(2*l-1,2),vv(ii,jj)
         
         endif
        enddo
       enddo         
       enddo
       
!   Proton-neutron part       
       vv=0.0d0
       uu=0.0d0
       vt=0.0d0
       
       do i=1,id
         do j=1,id 
           ii=id*(i-1)+j    
          do k=1,id
            do l=1,id
             jj=id*(k-1)+l
              uu(jj,ii)=Tp(k,i)*Tn(l,j)
            enddo
          enddo
         enddo
       enddo
       
       do Jp=0,jmax
       write(*,*)'Transformation of PN interaction for J =',Jp
        vv=0.d0
        do i=1,id 
         do j=1,id
          ii=id*(i-1)+j
          do k=1,id
           do l=1,id
             jj=id*(k-1)+l
              
              if (if_dd.eq.1) then 
              vv(ii,jj)=Vpn(i,j,k,l,Jp)+3.0d0*Vpn_dd(i,j,k,l,Jp)
              else 
              vv(ii,jj)=Vpn(i,j,k,l,Jp)
              endif
           enddo
          enddo
         enddo
        enddo
        
       call dgemm('N','N',id*id,id*id,id*id,1.d0,vv,id*id,uu,id*id,0.d0,vt,id*id)
       vv=0.d0
       call dgemm('T','N',id*id,id*id,id*id,1.d0,uu,id*id,vt,id*id,0.d0,vv,id*id)
       
      
       do ii=1,id*id
          i=int(ii/id)+1          
          j=mod(ii,id) 
              
          if (mod(ii,id).eq.0) then 
           i=ii/id
           j=id
          endif
        
 
        do jj=1,id*id
           k=int(jj/id)+1
           l=mod(jj,id) 
           
          if (mod(jj,id).eq.0) then 
           k=jj/id
           l=id
          endif
 
         if(100000*i+j.le.100000*k+l) then
             
!         if(i.eq.j) vv(ii,jj)=vv(ii,jj)/dsqrt(2.d0)
!         if(k.eq.l) vv(ii,jj)=vv(ii,jj)/dsqrt(2.d0)
!         if (dabs(vv(ii,jj)).gt.1.d-10) write(1000,*)0,0,2*Jp,2*i-1,2*j,2*k-1,2*l,vv(ii,jj)
!          if (dabs(vv(ii,jj)).gt.1.d-10) write(1)0,0,2*Jp,2*i-1,2*j,2*k-1,2*l,vv(ii,jj)
         if (dabs(vv(ii,jj)).gt.1.d-10) write(1)int(0,1),int(2*Jp,1),int(2*i-1,2),int(2*j,2),int(2*k-1,2),int(2*l,2),vv(ii,jj)          
         
         endif
        enddo
       enddo         
       enddo       
       
       
!   Neutron-neutron part       
       vv=0.0d0
       uu=0.0d0
       vt=0.0d0
       
       do i=1,id
         do j=1,id 
           ii=id*(i-1)+j    
          do k=1,id
            do l=1,id
             jj=id*(k-1)+l
              uu(jj,ii)=Tn(k,i)*Tn(l,j)
            enddo
          enddo
         enddo
       enddo
       
       do Jp=0,jmax
       write(*,*)'Transformation of NN interaction for J =',Jp
        vv=0.d0
        do i=1,id 
         do j=1,id
          ii=id*(i-1)+j
          do k=1,id
           do l=1,id
             jj=id*(k-1)+l
              if (if_dd.eq.1) then 
              vv(ii,jj)=Vnn(i,j,k,l,Jp)+3.0d0*Vnn_dd(i,j,k,l,Jp)
              else 
              vv(ii,jj)=Vnn(i,j,k,l,Jp)
              endif
           enddo
          enddo
         enddo
        enddo
        
       call dgemm('N','N',id*id,id*id,id*id,1.d0,vv,id*id,uu,id*id,0.d0,vt,id*id)
       vv=0.d0
       call dgemm('T','N',id*id,id*id,id*id,1.d0,uu,id*id,vt,id*id,0.d0,vv,id*id)
       
      
       do ii=1,id*id
          i=int(ii/id)+1          
          j=mod(ii,id) 
              
          if (mod(ii,id).eq.0) then 
           i=ii/id
           j=id
          endif
        
 
        do jj=1,id*id
           k=int(jj/id)+1
           l=mod(jj,id) 
           
          if (mod(jj,id).eq.0) then 
           k=jj/id
           l=id
          endif
 
         if((i.le.j.and.k.le.l).and.(100000*i+j.le.100000*k+l)) then
             
         if(i.eq.j) vv(ii,jj)=vv(ii,jj)/dsqrt(2.d0)
         if(k.eq.l) vv(ii,jj)=vv(ii,jj)/dsqrt(2.d0)
!         if (dabs(vv(ii,jj)).gt.1.d-10) write(1000,*)-1,0,2*Jp,2*i,2*j,2*k,2*l,vv(ii,jj)
!         if (dabs(vv(ii,jj)).gt.1.d-10) write(1)-1,0,2*Jp,2*i,2*j,2*k,2*l,vv(ii,jj)
         if (dabs(vv(ii,jj)).gt.1.d-10) write(1)int(1,1),int(2*Jp,1),int(2*i,2),int(2*j,2),int(2*k,2),int(2*l,2),vv(ii,jj)         
         
         endif
        enddo
       enddo         
       enddo
       
       close(1)
       
       
       
       

       
       deallocate(vv,uu,vt)
       
       

!       allocate(Vpp_HFB(id,id,id,id,0:jmax))
!       allocate(Vnn_HFB(id,id,id,id,0:jmax))
!       allocate(Vpn_HFB(id,id,id,id,0:jmax))
!       Vpp_HFB=0.d0
!       Vnn_HFB=0.d0
!       Vpn_HFB=0.d0

   


       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=-1
         do i=1,id
          if(lhfp(i)%j2.eq.jj.and.lhfp(i)%l.eq.ll) then
           nn=nn+1
           lhfp(i)%nn=nn
          endif
         enddo
        enddo
       enddo

       do ll=0,(jmax+1)/2
        do jj=1,jmax,2
         nn=-1
         do i=1,id
          if(lhfn(i)%j2.eq.jj.and.lhfn(i)%l.eq.ll) then
           nn=nn+1
           lhfn(i)%nn=nn
          endif
         enddo
        enddo
       enddo

       i1p=0
       i1n=0
       do i=1,id
        if(lhfp(i)%vi**2.d0.gt.0.98) then
         i1p=i1p+1
         point_p(i)=id+1-i1p
         pinv_p(id+1-i1p)=i
        else
         point_p(i)=i-i1p
         pinv_p(i-i1p)=i
        endif
        if(lhfn(i)%vi**2.d0.gt.0.98) then
         i1n=i1n+1
         point_n(i)=id+1-i1n
         pinv_n(id+1-i1n)=i
        else
         point_n(i)=i-i1n
         pinv_n(i-i1n)=i
        endif
       enddo

       open(2,file='proton_HF.dat',status='unknown',form='formatted')
        write(2,'(1x,a61)')'n,     l,    2*j,      Tz,     qei,  ui,      vi'
        itz=-1
        do ii = 1, id
         jj=ii !pinv_p(ii)
         write(2,'(1x,4(i4,1x),3(f12.5,1x))') lhfp(jj)%nn,lhfp(jj)%l,lhfp(jj)%j2,itz,lhfp(jj)%qei,DSQRT(1.D0-lhfp(jj)%vi**2.d0),lhfp(jj)%vi
        end do
       close(2)
       open(2,file='neutron_HF.dat',status='unknown',form='formatted')
        write(2,'(1x,a61)')'n,     l,    2*j,      Tz,     qei,  ui,      vi'
        itz=1
        do ii = 1, id
         jj=ii !pinv_n(ii)
         write(2,'(1x,4(i4,1x),3(f12.5,1x))') lhfn(jj)%nn,lhfn(jj)%l,lhfn(jj)%j2,itz,lhfn(jj)%qei,DSQRT(1.D0-lhfn(jj)%vi**2.d0),lhfn(jj)%vi
        end do
       close(2)

!       open(3,file='gmat_HF_pp.dat',status='unknown',form='formatted')
!       itz=-2
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_p(i1)
!             j=i2 !pinv_p(i2)
!             k=i3 !pinv_p(i3)
!             l=i4 !pinv_p(i4)
!             if(i.le.j.and.k.le.l) then
!              if(1000*i+j.le.1000*k+l) then
!               if(dabs(Vpp_HFB(i,j,k,l,Jp)).gt.precis) then
!                factab=1.d0
!                factcd=1.d0
!                if(i.eq.j) factab=dsqrt(2.d0)
!                if(k.eq.l) factcd=dsqrt(2.d0)
!                xnorm=factab*factcd
!                write(3,'(1x,14(i4,1x),1(f12.5,1x))') 
!     &            lhfp(i)%nn,lhfp(i)%l,lhfp(i)%j2,
!     &            lhfp(j)%nn,lhfp(j)%l,lhfp(j)%j2,
!     &            lhfp(k)%nn,lhfp(k)%l,lhfp(k)%j2,
!     &            lhfp(l)%nn,lhfp(l)%l,lhfp(l)%j2,
!     &            itz,2*Jp,Vpp_HFB(i,j,k,l,Jp)/xnorm
!               endif
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

!       open(3,file='gmat_HF_nn.dat',status='unknown',form='formatted')
!       itz=2
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_n(i1)
!             j=i2 !pinv_n(i2)
!             k=i3 !pinv_n(i3)
!             l=i4 !pinv_n(i4)
!             if(i.le.j.and.k.le.l) then
!              if(1000*i+j.le.1000*k+l) then
!               if(dabs(Vnn_HFB(i,j,k,l,Jp)).gt.precis) then
!                factab=1.d0
!                factcd=1.d0
!                if(i.eq.j) factab=dsqrt(2.d0)
!                if(k.eq.l) factcd=dsqrt(2.d0)
!                xnorm=factab*factcd
!                write(3,'(1x,14(i4,1x),1(f12.5,1x))')
!     &            lhfn(i)%nn,lhfn(i)%l,lhfn(i)%j2,
!     &            lhfn(j)%nn,lhfn(j)%l,lhfn(j)%j2,
!     &            lhfn(k)%nn,lhfn(k)%l,lhfn(k)%j2,
!     &            lhfn(l)%nn,lhfn(l)%l,lhfn(l)%j2,
!     &            itz,2*Jp,Vnn_HFB(i,j,k,l,Jp)/xnorm
!               endif
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

!       open(3,file='gmat_HF_pn.dat',status='unknown',form='formatted')
!       itz=0
!        do Jp=0,jmax
!         do i1=1,id
!          do i2=1,id
!           do i3=1,id
!            do i4=1,id
!             i=i1 !pinv_p(i1)
!             j=i2 !pinv_n(i2)
!             k=i3 !pinv_p(i3)
!             l=i4 !pinv_n(i4)
!             if(1000*i+j.le.1000*k+l) then
!              if(dabs(Vpn_HFB(i,j,k,l,Jp)).gt.precis) then
!               write(3,'(1x,14(i4,1x),1(f12.5,1x))')
!     &           lhfp(i)%nn,lhfp(i)%l,lhfp(i)%j2,
!     &           lhfn(j)%nn,lhfn(j)%l,lhfn(j)%j2,
!     &           lhfp(k)%nn,lhfp(k)%l,lhfp(k)%j2,
!     &           lhfn(l)%nn,lhfn(l)%l,lhfn(l)%j2,
!     &           itz,2*Jp,Vpn_HFB(i,j,k,l,Jp)
!              endif
!             endif
!            enddo
!           enddo
!          enddo
!         enddo
!        enddo
!       close(3)

       open(4,file='unitar_p.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (Tp(i1,j1),j1=1,id) !(Tp(i1,pinv_p(j1)),j1=1,id)
        enddo
       close(4)

       open(4,file='unitar_n.dat',status='unknown',form='formatted')
        do i1=1,id
         write(4,*) (Tn(i1,j1),j1=1,id) !(Tn(i1,pinv_n(j1)),j1=1,id)
        enddo
       close(4)

!       deallocate(Vpp_HFB,Vnn_HFB,Vpn_HFB)

       return
      end
