      Subroutine NatOrb_transf

       USE technical
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

!       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp_HFB(:,:,:,:,:)
!       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn_HFB(:,:,:,:,:)
!       DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpn_HFB(:,:,:,:,:)


       double precision, allocatable :: vv(:,:),uu(:,:),vt(:,:)
       double precision :: Tp(id,id),Tn(id,id)
       double precision :: bp(id,id),bn(id,id)

       Tp=tran_p
       Tn=tran_n

       open(1,file='vlk_nat_orb.dat',status='unknown',form='unformatted')
       
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
              vv(ii,jj)=Vpp(i,j,k,l,Jp)!+3.0d0*Vpp_dd(i,j,k,l,Jp)
              else
              vv(ii,jj)=Vpp(i,j,k,l,Jp)
              endif 
           enddo
          enddo
         enddo
        enddo
        
       vt=0.d0
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
         if (dabs(vv(ii,jj)).gt.1.d-7) write(1)int(-1,1),int(2*Jp,1),int(2*i-1,2),int(2*j-1,2),int(2*k-1,2),int(2*l-1,2),vv(ii,jj)
!         if (dabs(vv(ii,jj)).gt.1.d-10) write(999,'(6i5,f15.10)')-1,2*Jp,2*i-1,2*j-1,2*k-1,2*l-1,vv(ii,jj) 
!         write(999,*)int(-1,1),int(2*Jp,1),int(2*i-1,2),int(2*j-1,2),int(2*k-1,2),int(2*l-1,2),vv(ii,jj)
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
              vv(ii,jj)=Vpn(i,j,k,l,Jp)!+3.0d0*Vpn_dd(i,j,k,l,Jp)
              else 
              vv(ii,jj)=Vpn(i,j,k,l,Jp)
              endif
           enddo
          enddo
         enddo
        enddo

!        write(994,*) '  Transf.   -----------------'        
!       do i=1,id*id 
!        write(994,'(100f8.3)')(uu(i,j),j=1,id*id)
!       enddo

!       write(994,*) '  V   -----------------'

!       do i=1,id*id 
!        write(994,'(100f8.3)')(vv(i,j),j=1,id*id)
!       enddo

        
       vt=0.d0
       call dgemm('N','N',id*id,id*id,id*id,1.d0,vv,id*id,uu,id*id,0.d0,vt,id*id)
       vv=0.d0
       call dgemm('T','N',id*id,id*id,id*id,1.d0,uu,id*id,vt,id*id,0.d0,vv,id*id)


!       write(994,*) ' transformed V-----------------'

!       do i=1,id*id 
!        write(994,'(100f8.3)')(vv(i,j),j=1,id*id)
!       enddo
       
      
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
         if (dabs(vv(ii,jj)).gt.1.d-7) write(1)int(0,1),int(2*Jp,1),int(2*i-1,2),int(2*j,2),int(2*k-1,2),int(2*l,2),vv(ii,jj)
!         if (dabs(vv(ii,jj)).gt.1.d-1) write(999,'(6i5,f15.10,2i5)')0,2*Jp,2*i-1,2*j,2*k-1,2*l,vv(ii,jj),ii,jj 
!         write(999,*)int(0,1),int(2*Jp,1),int(2*i-1,2),int(2*j,2),int(2*k-1,2),int(2*l,2),vv(ii,jj)          
         
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
              vv(ii,jj)=Vnn(i,j,k,l,Jp)!+3.0d0*Vnn_dd(i,j,k,l,Jp)
              else 
              vv(ii,jj)=Vnn(i,j,k,l,Jp)
              endif
           enddo
          enddo
         enddo
        enddo
      
       vt=0.d0
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
         if (dabs(vv(ii,jj)).gt.1.d-7) write(1)int(1,1),int(2*Jp,1),int(2*i,2),int(2*j,2),int(2*k,2),int(2*l,2),vv(ii,jj)         
  !       if (dabs(vv(ii,jj)).gt.1.d-10) write(999,'(6i5,f15.10)')1,2*Jp,2*i,2*j,2*k,2*l,vv(ii,jj) 
  !       write(999,*)int(1,1),int(2*Jp,1),int(2*i,2),int(2*j,2),int(2*k,2),int(2*l,2),vv(ii,jj)
         
         endif
        enddo
       enddo         
       enddo
       
       close(1)
       
       
       
       !!!!!

       
       deallocate(vv,uu,vt)





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
       open(165,file='kin_nat_orb.dat',status='unknown',form='unformatted')
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

       open(61,file='E0_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='E0_NO_p.dat',status='unknown',form='formatted')
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

       open(61,file='r1Y1_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r1Y1_NO_p.dat',status='unknown',form='formatted')
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

       open(61,file='r2Y2_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r2Y2_NO_p.dat',status='unknown',form='formatted')
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

       open(61,file='r3Y3_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='r3Y3_NO_p.dat',status='unknown',form='formatted')
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

       open(61,file='EN_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='EN_NO_p.dat',status='unknown',form='formatted')
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

       open(61,file='S1_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='S1_NO_p.dat',status='unknown',form='formatted')
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

       open(61,file='M1s_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='M1s_NO_p.dat',status='unknown',form='formatted')
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

       open(61,file='M1l_NO_n.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
         if (dabs(bn(i,j)).gt.1.d-10) then
          write(61,'(2i5,5x,2f30.15)')i,j,bn(i,j)
         endif
        enddo
       enddo
       close(61)

       open(61,file='M1l_NO_p.dat',status='unknown',form='formatted')
       do i=1,id
        do j=1,id
        if (dabs(bp(i,j)).gt.1.d-10) then
         write(61,'(2i5,5x,2f30.15)')i,j,bp(i,j)
        endif
        enddo
       enddo
       close(61)




       return
      end
