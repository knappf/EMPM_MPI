
!     Cholesky    procedure      
!     LAPACK implementation 
!     last modification 25.3.2020

   module choleski_decomp

      contains
       
      subroutine cholesk_dec(ndim,dd,nx,mxt)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
     
      
      double precision, dimension(:,:), allocatable :: dd,cq,d1,cdu
    
      double precision, dimension(:), allocatable :: r,rcp 

      integer, dimension(:), allocatable :: nx,mxt


      allocate(d1(ndim,ndim))
      d1=0.d0
      allocate(cdu(ndim,ndim))
      cdu=0.d0
      allocate(rcp(ndim+1),r(ndim),nx(ndim),mxt(ndim))
      rcp=0.d0
      r=0.d0
      nx=0
      mxt=0
      
      norder=ndim

      

      open(6,file='mxt.dat',status='unknown',form='unformatted')

      write(6)ndim,no
      do i=1,ndim
        write(6)i,mxt(i)
      enddo
      close(6)

      write(*,*)' Number of linearly independent states ',no


      deallocate(cdu)
 
      return 
      end subroutine cholesk_dec

      subroutine read_sub_dmat(dmatr,ndim)

            implicit double precision (a-h,o-z)
      
      
            double precision, dimension(:,:), allocatable :: dmatr
            double precision, dimension(:), allocatable :: dmm
      

            if (.not.allocated(dmatr))  allocate(dmatr(ndim,ndim))
            dmatr=0.d0
      
            open(66,file='d_mat.dat',status='old',form='unformatted')
      
            read(66)ndimrt,ndimt
            allocate(dmm(ndimt))
            dmm=0.0d0
      
      !      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
      !        write(*,*)' Dimensions does not match in d_m file'
      !        stop
      !      endif
            do iii=1,ndim
             read(66)(dmm(jjj),jjj=1,ndim)
               do jjj=1,ndim
                 dmatr(iii,jjj)=dmm(jjj)
               enddo 
            enddo
      
            close(66)
            deallocate(dmm)
      
            return
      
            end subroutine read_sub_dmat

 

      end module choleski_decomp
     
