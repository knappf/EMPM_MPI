      subroutine kinetic

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       allocate(kin_p(id,id),kin_n(id,id))
       kin_p=0.d0
       kin_n=0.d0

       aprot = 0.5d0*(mp+mn)/mp
       aneut = 0.5d0*(mp+mn)/mn

       do j=1,id
        kin_p(j,j) = aprot*0.5d0*hbarom*(dble(levp(j)%N)+1.5d0)
!     &                                        *dsqrt(dble(levp(j)%j2+1))
        kin_n(j,j) = aneut*0.5d0*hbarom*(dble(levn(j)%N)+1.5d0)
!     &                                        *dsqrt(dble(levn(j)%j2+1))
        do k=1,id
         if(levp(j)%l.eq.levp(k)%l) then
          if(levp(j)%j2.eq.levp(k)%j2) then
           if(levp(j)%nn.eq.levp(k)%nn+1) then
            kin_p(j,k) = aprot*0.5d0*hbarom*dsqrt(dble(levp(j)%nn)*     
     &     (dble(levp(j)%nn+levp(j)%l)+0.5d0))!*dsqrt(dble(levp(j)%j2+1))
           endif
           if(levp(j)%nn+1.eq.levp(k)%nn) then
            kin_p(j,k) = aprot*0.5d0*hbarom*dsqrt(dble(levp(k)%nn)*
     &     (dble(levp(k)%nn+levp(k)%l)+0.5d0))!*dsqrt(dble(levp(j)%j2+1))
           endif
          endif
         endif
         if(levn(j)%l.eq.levn(k)%l) then
          if(levn(j)%j2.eq.levn(k)%j2) then
           if(levn(j)%nn.eq.levn(k)%nn+1) then
            kin_n(j,k) = aneut*0.5d0*hbarom*dsqrt(dble(levn(j)%nn)*
     &     (dble(levn(j)%nn+levn(j)%l)+0.5d0))!*dsqrt(dble(levn(j)%j2+1))
           endif
           if(levn(j)%nn+1.eq.levn(k)%nn) then
            kin_n(j,k) = aneut*0.5d0*hbarom*dsqrt(dble(levn(k)%nn)*
     &     (dble(levn(k)%nn+levn(k)%l)+0.5d0))!*dsqrt(dble(levn(j)%j2+1))
           endif
          endif
         endif
        enddo
       enddo

       do i=1,id
        do j=1,id
         kin_p(i,j)=kin_p(i,j)*(dble(AZ+AN-1)/dble(AZ+AN))
         kin_n(i,j)=kin_n(i,j)*(dble(AZ+AN-1)/dble(AZ+AN))
        enddo
       enddo

       return
      end
