!     Cholesky    procedure      
!     LAPACK implementation 
!     last modification 25.3.2020

module choleski

 use read_admat      

 contains
       
 subroutine cholesk(i_spur,ndim,no,noo,dd,nx,mxt)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
     
      
      double precision, dimension(:,:), allocatable :: dd
    
      double precision, dimension(:), allocatable :: r,work 

      integer, dimension(:), allocatable :: nx,mxt
 
      if (i_spur.eq.0) then
       allocate(dd(ndim,ndim))
       dd=0.d0
       call read_dmat(dd)
      endif  
 
      allocate(work(10*ndim))
      allocate(mxt(ndim),nx(ndim))

      mxt=0
      nx=0 

      
      tol=0.0001d0

      write(*,*)' dimension = ',ndim
            
      call dpstrf('U', ndim, dd, ndim, mxt, no, tol, work, info )

      write(*,*)'Cholesky info =', info


      do i=1,no
        nx(mxt(i))=1
      enddo

      open(6,file='mxt.dat',status='unknown',form='unformatted')

      write(6)ndim,no
      do i=1,ndim
        write(6)i,mxt(i)
      enddo
      close(6)

      
      write(*,*)' Number of linearly independent states ',no


 
      return 
 end subroutine cholesk

end module choleski
     
