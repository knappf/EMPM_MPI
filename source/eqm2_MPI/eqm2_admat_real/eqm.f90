!     last modification 11.5.2018      
      
      program eqm 

      use types_eqm  
      use phonon_base
!      use metricmat
      use admatr

      implicit double precision (a-h,o-z)

      include 'mpif.h'
!      include 'types_eqm.inc'

      type(phonbase_typ), dimension (:), allocatable :: phonbs
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      integer, dimension (:), allocatable :: phonus, phonmus
      type (amp2_typ), dimension(:,:), allocatable :: xcc

!     choleski arrays
      double precision, dimension(:,:), allocatable :: dd,cq,d1,cdu,xr,vr,h_corr
      double precision, dimension(:), allocatable :: wr
      
      integer, dimension(:), allocatable :: nx,mxt,irow,mxtr
      integer, dimension (:,:), allocatable ::ndbx
      integer, dimension (:,:,:), allocatable ::  ipozx
      integer, dimension (:), allocatable :: ndcc,i_resh
      
      character*30 namex,names,namec,namex2
 
      integer :: myid,ierr,numprocs


!      CALL OMP_SET_NUM_THREADS(8)

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      
!      if (myid.eq.0) open(992,file='phonon_base.dat',status='unknown',form='formatted')
      nf=2  
       
!      namex2='2phonon/2f_x.dat'
!      call readx(namex2,xcc,ndcc,nd1f,nd2f)
!      call redrsumx(nd2f,nd1f,xcc,ndcc,ipozx,ndbx)

      idim1=500000
      idim2=500000
      idimbs=50000000
      
      open(23,file='AD_J_Pi.dat',status='old',form='formatted')     
      read(23,*)ipar,jcal    
      close(23)

      if (myid.eq.0) then 
        write(*,*)
        write(*,*)'----------------------------------------------'
        write(*,*)' Parity = ',ipar,'   J = ',jcal
      endif

!      if (myid.eq.0) then 
!        write(992,*)
!        write(992,*)'----------------------------------------------'
!        write(992,*)' Parity = ',ipar,'   J = ',jcal
!      endif 

      call phonbase(nf,ipar,jcal,phonus,phonmus,idim1,idim2,idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr,n_spur,myid)

      if (myid.eq.0) then 
        write(*,*)' Dimension = ',idphon
        write(*,*)' Truncated dimension = ',idphontr
      endif

      idphontot=idphon

!      if (idphon.gt.0) then 

!        call MPI_BARRIER(  MPI_COMM_WORLD, ierr)
! calculation of matrices A and D
!       allocate(i_resh(idphontr))
!       do i=1,idphontr
!         i_resh(i)=i
!       enddo

!      call rperm(idphontr, i_resh)

        call admat(nf,ipar,jcal,phonbs,idphon,no,phon1,phon2,mxt,myid,numprocs,i_resh)

!       CALL execute_command_line('./run_admat.sh' )
!       stop


      close(99)
      close(12)
      close(22)
!      close(992)
    
      call MPI_FINALIZE(irc)
      end
