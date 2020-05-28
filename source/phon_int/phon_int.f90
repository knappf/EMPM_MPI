!*     
!*     Program phon_int computes phonon interaction in proton-neutron 
!*     J-coupled formalism.

!*     last update 14.3.2018

program phon_int 

use input_sp
use phoninterac

  implicit double precision (a-h,o-z)

  include 'mpif.h'


  include 'input_phon_int.inc'
  include 'formats_phon_int.inc'
  include 'types_phon_int.inc'

  type(level_typ),dimension(:), allocatable :: levn,levp
  character*2 :: type_phon


   integer :: myid,ierr,numprocs

!     loading of input data
call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      
!     loading of input data 

      xrotrunc=1.d-10
      
      if (myid.eq.0) write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
!      read(1,26)alfa,beta
!      read(1,*)
!      read(1,15)iparmn,iparmx
!      read(1,15)jminn,jmaxn
      close(1)

      allocate(levn(ipnmx+1000),levp(ippmx+1000))

      call inp_sp(levn,levp,myid)

      jmaxn=0
      do i=1,ippmx
       if (levp(i)%j.gt.jmaxn) jmaxn=levp(i)%j
      enddo

      do i=1,ipnmx
       if (levn(i)%j.gt.jmaxn) jmaxn=levn(i)%j
      enddo

      if (myid.eq.0) write(*,*)'jmax =',jmaxn


      jmax=jmaxn
      
!c      call checkf
!c      call checkfmat
!c      call checkro
      
!c      stop
      
      
      ifmx=5000

      if (myid.eq.0) write(*,*)'energy threshold for 1 phonon states'
!     read(*,*)xthrun_min,xthrun_max
      xthrun_min=0.0
      xthrun_max=1000000000.0
      if (myid.eq.0) write(*,*)xthrun_min,xthrun_max
 
!      write(*,*)'Type of phonon-phonon intetaction: pp,np,ph,nh'
!      read(*,*)type_phon
!      write(*,*)type_phon


call vintp(1,ifmx,jmax,levn,xthrun_min,xthrun_max,myid,numprocs)  ! neutron particle

call MPI_BARRIER(  MPI_COMM_WORLD, ierr)

call vintp(-1,ifmx,jmax,levp,xthrun_min,xthrun_max,myid,numprocs)

call MPI_BARRIER(  MPI_COMM_WORLD, ierr)

call vinth(1,ifmx,jmax,levn,xthrun_min,xthrun_max,myid,numprocs)

call MPI_BARRIER(  MPI_COMM_WORLD, ierr)

call vinth(-1,ifmx,jmax,levp,xthrun_min,xthrun_max,myid,numprocs)

!      if (type_phon.eq.'np') call vintp(1,ifmx,jmax,levn,xthrun_min,xthrun_max)  ! neutron particle 
!      if (type_phon.eq.'pp') call vintp(-1,ifmx,jmax,levp,xthrun_min,xthrun_max) ! proton particle
!      if (type_phon.eq.'ph') call vinth(-1,ifmx,jmax,levp,xthrun_min,xthrun_max) ! proton hole
!      if (type_phon.eq.'nh') call vinth(1,ifmx,jmax,levn,xthrun_min,xthrun_max)  ! neutron hole


call MPI_FINALIZE(irc)      
end 
!*     
!*     END of the main program 

!* 

