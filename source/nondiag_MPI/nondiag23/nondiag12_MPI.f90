!     last update 14.6.2013

      program nondiag12
      
      use phoninteracn
      use input_sp

      implicit double precision (a-h,o-z)
      
      include 'types_ndgi_int.inc'
      include 'input_ndgi_int.inc'
      include 'formats_ndgi_int.inc'

      include 'mpif.h'
      
      type(level_typ),dimension(:), allocatable :: levn,levp 
      
      integer :: myid,ierr,numprocs

      !     loading of input data
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      
      
      write(*,*)'Parity J ?'
      open(myid,file='input_nondg',status='old',form='formatted')
      read(myid,*)ipar,jcal
      if (myid ==0) write(*,*)ipar,jcal
      close(myid)
            
      call loadsp(levn,levp,myid)
      call inp_sp(levn,levp,myid)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call vintn23(3,ipar,jcal,levn,levp,myid,numprocs)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (myid ==0 ) then
            CALL execute_command_line('cat scratch/Vint23_* > Vint_phon23.dat' )
            CALL execute_command_line('rm scratch/Vint23_*' )      
      endif

      call MPI_FINALIZE(ierr)



      end
