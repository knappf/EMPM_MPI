!  last update 28.7. 2011

module cmconst
!
 contains 
!
!
 subroutine cmcreate(idphp,idphn,iphp,iphn,tbase)

  implicit double precision (a-h,o-z)

  include 'types_tda_cp.inc'
 
  double precision, dimension (:,:,:), allocatable :: muel
  double precision, dimension (:,:), allocatable :: tbase
  type(ph_typ),dimension(:), allocatable :: iphn,iphp
 
  call read_mul(muel)

  ndim=idphp+idphn
  allocate(tbase(ndim,ndim))
  tbase=0.0d0

  do i=1,idphp
      ip=iphp(i)%par
      ih=iphp(i)%hol
      tbase(i,1)=muel(ip,ih,1)
  enddo

  do i=1,idphn
      ip=iphn(i)%par
      ih=iphn(i)%hol
      tbase(i+idphp,1)=muel(ip,ih,2)
  enddo

  do j=2,ndim
   tbase(j,j)=1.0d0
  enddo

 end subroutine cmcreate
!
 subroutine read_mul(muel)
 
  implicit double precision (a-h,o-z)
   
  include 'formats_tda_cp.inc'


  double precision, dimension (:,:,:), allocatable :: muel
  character*30 :: fname 

  imax=200
  allocate(muel(imax,imax,2))
  muel=0.d0

  fname='r1Y1_p.dat'
  open(21,file=fname,status='unknown',form='formatted')

    do while (.not.eof(21))
      read(21,1002)i,j,rs
      muel(i,j,1)=rs
    enddo

    close(21)

  fname='r1Y1_n.dat'
  open(21,file=fname,status='unknown',form='formatted')

    do while (.not.eof(21))
      read(21,1002)i,j,rs
      muel(i,j,2)=rs
    enddo
     
    close(21)
 

 end subroutine read_mul
!
 subroutine ortog(ndim,tbase)
 
  implicit double precision (a-h,o-z)

  double precision, dimension (:,:), allocatable :: tbase 
  double precision, dimension (:), allocatable :: a
 
 
  allocate(a(ndim))
  a=0.0d0
! normalization of first vector      

      xnorm=0.d0

      do i=1,ndim
        xnorm=xnorm+tbase(i,1)**2.0d0
      enddo
      
      do i=1,ndim
        tbase(i,1)=tbase(i,1)/dsqrt(xnorm)
      end do

      do i=2,ndim
      
      do j=1,i-1
          xovrl=0.d0
        do k=1,ndim
          xovrl=xovrl+tbase(k,j)*tbase(k,i)
        end do
          a(j)=-1.d0*xovrl      
      end do
      
      do jj=1,i-1
         do ii=1,ndim
           tbase(ii,i)=tbase(ii,i)+a(jj)*tbase(ii,jj)
         end do
      end do
       xnorm=0.d0
      do ii=1,ndim
           xnorm=xnorm+tbase(ii,i)*tbase(ii,i)
      end do
      do ii=1,ndim
           tbase(ii,i)=tbase(ii,i)/dsqrt(xnorm)
      end do
      
      end do    



!      if (nf.ne.1) then
!      xnorm=0.d0
!      do i=1,no
!      do j=1,no
!      xnorm=xnorm+dd(i,1)*d1(i,j)*dd(j,1)
!      enddo
!      enddo
      
!      do i=1,no
!      dd(i,1)=dd(i,1)/dsqrt(xnorm)
!      end do

!      do i=2,no
      
!      do j=1,i-1
!      xovrl=0.d0
!      do k=1,no
!      do l=1,no
!      xovrl=xovrl+dd(k,j)*d1(k,l)*dd(l,i)
!      end do
!      end do
!      a(j)=-1.d0*xovrl      
!      end do
      
!      do jj=1,i-1
!      do ii=1,no
!      dd(ii,i)=dd(ii,i)+a(jj)*dd(ii,jj)
!      end do
!      end do
!      xnorm=0.d0
!      do ii=1,no
!      do jj=1,no
!      xnorm=xnorm+dd(ii,i)*d1(ii,jj)*dd(jj,i)
!      end do
!      end do
!      do ii=1,no
!      dd(ii,i)=dd(ii,i)/dsqrt(xnorm)
!      end do
      
!      end do    
!      endif 
           
      return 
      end subroutine ortog
!
!
   subroutine ch_fmat(fp,fn,fpn,idphp,idphn,iphp,iphn)

      include 'types_tda_cp.inc'
      include 'input_tda_cp.inc'

      double precision, dimension(:,:,:,:), allocatable :: fp,fn,fpn
      double precision, dimension (:,:), allocatable :: tbase
      type(ph_typ),dimension(:), allocatable :: iphn,iphp

      call cmcreate(idphp,idphn,iphp,iphn,tbase)
  
       do i=1,idphp
            ip=iphp(i)%par
            ih=iphp(i)%hol
            do j=1,idphp
              ipp=iphp(j)%par
              ihp=iphp(j)%hol

      
!      amat=0.d0

      
        
      
!      ifz=(-1)**((levp(ipp)%j+levp(ihp)%j)/2+ijj)
!      amat=amat-fp(ip,ih,ihp,ipp)*dfloat(ifz)
      

!              amtr(i,j)=amat
            enddo
          enddo





   end subroutine ch_fmat
      
end module cmconst







