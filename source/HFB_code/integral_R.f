!**************************************************************************
!  The numerical integration based on GAUSS' MECHANICAL QUADRATURE FORMULA*
!  coded by Giovanni De Gregorio on 19 September 2013                     *
!**************************************************************************
      function R_integral(i1,j1,k1,l1,alpha,rad_den)

       USE technical

       include 'define.inc'
       include 'parameters.inc' 
       include 'commons.inc'

       double precision :: rad_den(igrid)
       double precision, dimension(igrid) :: argu
!       double precision, dimension(9) :: zcross, weight   !zero-crossing and weight
       
       bos=dsqrt(0.5d0*(mp+mn)*hbarom/(hbarc**2.d0))

!       allocate(argu(0:igrid))

       argu=0.d0
       xcross=0.d0                          
   
       do i=1,igrid
        argu(i)=R_field(levp(i1)%nn,levp(i1)%l,levp(i1)%j2,i)
     & *R_field(levp(j1)%nn,levp(j1)%l,levp(j1)%j2,i)
     & *R_field(levp(k1)%nn,levp(k1)%l,levp(k1)%j2,i)
     & *R_field(levp(l1)%nn,levp(l1)%l,levp(l1)%j2,i)
     &  *zcross(i)**2.d0*rad_den(i)**alpha
       enddo
     
       R_integral=0.d0
       
       do i=1,igrid
!       write(*,*)'arg',argu(i)
!       write(*,*) 'weight', weight(i)
          R_integral=R_integral+weight(i)*argu(i)
       enddo

!       deallocate(argu)

       return
      end
