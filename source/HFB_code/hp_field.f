      Subroutine hp_field(h1)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: h1(id,id)

       h1=0.d0
       do i=1,id
        do j=1,id
         val=0.d0
         do k=1,id
          do l=1,id
           if(levp(i)%j2.eq.levp(j)%j2
     &                .and.levp(k)%j2.eq.levp(l)%j2) then
            val=val+(Fpp(j,i,k,l)+3.d0*Fpp_DD(j,i,k,l)/2.d0)
     &               *rhop_HFB(l,k)*dsqrt(dble(levp(k)%j2+1))
     &             +(Fpn(j,i,k,l)+3.d0*Fpn_DD(j,i,k,l)/2.d0)
     &               *rhon_HFB(l,k)*dsqrt(dble(levn(k)%j2+1))
           endif
          enddo
         enddo
!         h1(i,j)=kin_p(j,i)*dsqrt(dble(levp(i)%j2+1))+val
         h1(i,j)=kin_p(j,i)+val/dsqrt(dble(levp(i)%j2+1))
        enddo
       enddo

       return
      end
