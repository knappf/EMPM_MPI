!     
!     Program phon_dens1 computes 1phonon densities in proton-neutron 
!     J-coupled formalism.

!     last update 22.2.2018

program phon_dens1
 
   use input_sp
   use eofmod
   use dens1 

   implicit double precision (a-h,o-z)


   include 'mpif.h'
   include 'input_phon_dens1.inc'
   include 'formats_phon_dens1.inc'
   include 'types_phon_dens1.inc'

   type(level_typ),dimension(:), allocatable :: levn,levp
   integer, dimension(:), allocatable :: jphon
   integer io
   double precision, dimension(:,:,:), allocatable :: cph_p,cph_n

   integer :: myid,ierr,numprocs
      
!     loading of input data 
call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      xrotrunc=1.d-8
      
      if (myid.eq.0) write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      close(1)

      allocate(levn(2*ipnmx),levp(2*ippmx))

      call inp_sp(levn,levp,myid)

      jmaxn=0
      do i=1,ippmx
       if (levp(i)%j.gt.jmaxn) jmaxn=levp(i)%j
      enddo

      do i=1,ipnmx
       if (levn(i)%j.gt.jmaxn) jmaxn=levn(i)%j
      enddo

      if (myid.eq.0) write(*,*)'jmax =',jmaxn

      allocate(jphon(20000))

      open (4,file='1phonon/1f_states.dat',status='old',form='unformatted')

!      do while (.not.eof(4)) 
       do while (io.ne.-1) 
       read(4,iostat=io)i,ipar,ijj,en       
       if (i.gt.20000) then 
          write(*,*)' Increase dimension of jphon '
          stop
       endif
       jphon(i)=ijj
      enddo

      close(4)



      ifmx=i ! treba upravit aby sa to nacitalo 
      jmax=jmaxn
      if (myid.eq.0) then 
          write(*,*)' # of 1-phononon states ',ifmx
          write(*,*)' Jmax= ',jmax
      endif 

!      write(*,*)' Energy threshodl for 1 phonon states?'
!      read(*,*)xthres
!      write(*,*)xthres
      xthres=100000.0

      call read_cph(ifmx,cph_p,cph_n,ippmn,ippmx,ihpmn, &
&ihpmx,ipnmn,ipnmx,ihnmn,ihnmx)

!    call MPI_BARRIER(  MPI_COMM_WORLD, ierr)

    call rop(-1,ifmx,jmax,levp,xthres,jphon,cph_p,cph_n,myid,numprocs)

    call MPI_BARRIER(  MPI_COMM_WORLD, ierr)


    call rop(1,ifmx,jmax,levn,xthres,jphon,cph_p,cph_n,myid,numprocs)

    call MPI_BARRIER(  MPI_COMM_WORLD, ierr)


    call roh(-1,ifmx,jmax,levp,xthres,jphon,cph_p,cph_n,myid,numprocs)


    call MPI_BARRIER(  MPI_COMM_WORLD, ierr)

    call roh(1,ifmx,jmax,levn,xthres,jphon,cph_p,cph_n,myid,numprocs)



call MPI_FINALIZE(irc)      

end 
!     
!     END of the main program 
! 

      
