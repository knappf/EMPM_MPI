      subroutine Dn_field(D1)

       USE technical
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: D1(id,id)

       D1=0.d0
       do i=1,id
        do j=1,id
         if(levn(i)%j2.eq.levn(j)%j2.and.levn(i)%l.eq.levn(j)%l) then
          val=0.d0
          do k=1,id
           do l=1,id
           
        if (if_dd.eq.1) then            
           
            val=val+(Vnn(k,l,i,j,0)+3.d0*Vnn_DD(k,l,i,j,0)/2.d0)
     &                *kapn_HFB(l,k)*dsqrt(dble(levn(k)%j2+1))
        else        
            val=val+Vnn(k,l,i,j,0)
     &                *kapn_HFB(l,k)*dsqrt(dble(levn(k)%j2+1))
        endif
        
           enddo
          enddo
          D1(i,j)=0.5d0*val/dsqrt(dble(levn(i)%j2+1))
         endif
        enddo
       enddo

       return
      end
