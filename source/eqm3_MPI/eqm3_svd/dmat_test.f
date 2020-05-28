c     last modification 28.1.2010

      subroutine dmat_test(ndim)

      use choleski

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'

      double precision, dimension(:,:), allocatable :: dmatr,amatr
      double precision, dimension(:), allocatable :: work,e

      
      allocate(dmatr(ndim,ndim))
      dmatr=0.d0

  99  format(f15.6)    

      open(999,file='dmat.log',status='unknown',form='formatted')
c      open(6,file='dmat.dat',status='old',form='unformatted')

c      do while (.not.eof(6))
c      read(6)i,j,dd
c      dmatr(i,j)=dd
c      enddo


       call read_dmat(dmatr)

      do i=1,ndim

c      write(997,102)(dmatr(i,j),j=1,ndim)

       do j=i,ndim

       if (dabs(dmatr(i,j)-dmatr(j,i)).gt.1.d-6) 
     *write(*,*)'D nonsymmetric ',i,j,dmatr(i,j),dmatr(j,i)
       
       enddo
      enddo 

      lwork=26*ndim         
      allocate(work(26*ndim))
      allocate(e(ndim))
          
          
      call DSYEV('V','L',ndim,dmatr,ndim,e,WORK,LWORK,INFO )

      write(999,99)(e(i),i=1,ndim)

      close(999)
      
      deallocate(dmatr,work,e)

      end subroutine dmat_test
