!     last modification 11.5.2018      
      
      program eqm 

      use types_eqm  
      use phonon_base
!      use metricmat
      use admatr
      use input_sp

      implicit double precision (a-h,o-z)

      include 'mpif.h'
!      include 'types_eqm.inc'
      include 'formats_eqm.inc'

      type(phonbase_typ), dimension (:), allocatable :: phonbs
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      integer, dimension (:), allocatable :: phonus, phonmus
      type (amp2_typ), dimension(:,:), allocatable :: xcc
      type(rho2_typ), dimension(:,:), allocatable :: c_np,c_pp,c_nh,c_ph
      integer,dimension(:),allocatable :: idimnp,idimnh,idimph,idimpp
!     choleski arrays
      double precision, dimension(:,:), allocatable :: dd,cq,d1,cdu,xr,vr,h_corr
      double precision, dimension(:), allocatable :: wr
      type(level_typ),dimension(:), allocatable :: levn,levp
      
      integer, dimension(:), allocatable :: nx,mxt,irow,mxtr
      integer, dimension (:,:), allocatable ::ndbx
      integer, dimension (:,:,:), allocatable ::  ipozx
      integer, dimension (:), allocatable :: ndcc,i_resh
      
      character*30 namex,names,namec,namex2
      character*2:: ch2
 
      integer :: myid,ierr,numprocs


!      CALL OMP_SET_NUM_THREADS(8)

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            read(1,15)ia,iz
            read(1,15)ihnmn,ihnmx
            read(1,15)ihpmn,ihpmx
            read(1,15)ipnmn,ipnmx
            read(1,15)ippmn,ippmx
      close(1)


!      if (myid.eq.0) open(992,file='phonon_base.dat',status='unknown',form='formatted')
      nf=3  
       
      namex2='2phonon/2f_x.dat'
!      call readx(namex2,xcc,ndcc,nd1f,nd2f)
      call read_amp_all(namex2,xcc,ndcc,nd1f,nd2f)
      call redrsumx(nd2f,nd1f,xcc,ndcc,ipozx,ndbx)

      idim1=nd2f
      idim2=nd2f
      idimbs=5000000
      
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

      jmax2=0
      do ia=1,idim2
      if(phon2(ia)%j.gt.jmax2) jmax2=phon2(ia)%j
      end do
      write(*,*) 'Jmax 2ph=',jmax2  

      allocate(idimnp(0:jmax2),idimnh(0:jmax2),idimpp(0:jmax2),idimph(0:jmax2))
      idimnp=0
      idimnh=0
      idimph=0
      idimpp=0
      ndro=15000000
      allocate(c_pp(0:jmax2,ndro),c_nh(0:jmax2,ndro),c_np(0:jmax2,ndro),c_ph(0:jmax2,ndro))
      allocate(levn(200),levp(200))
      call inp_sp(levn,levp)

      jmaxn=0
      do i=1,ippmx
       if (levp(i)%j.gt.jmaxn) jmaxn=levp(i)%j
      enddo

      do i=1,ipnmx
       if (levn(i)%j.gt.jmaxn)jmaxn=levn(i)%j
      enddo

       write(*,*)'jmax=',jmaxn

       jmax=jmaxn
       isimax=jmax

      do jj=0,jmax2

      do ib=1,idim2
      
      jb=phon2(ib)%j
       do isi=0,isimax
        if(abs(jb-jj).le.isi.and.isi.le.(jb+jj))then
         !write(19,*) !iaaa,ib,isi,ja,jb!ib,iaaa,isi,ja,jb
         do i_2=ipnmn,ipnmx
          do i_1=ipnmn,ipnmx
           if(abs(levn(i_1)%j-levn(i_2)%j).le.2*isi.and.2*isi.le.(levn(i_1)%j+levn(i_2)%j)) then
           idimnp(jj)=idimnp(jj)+1
           if(idimnp(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
           if(idimnp(jj).gt.ndro) stop
           c_np(jj,idimnp(jj))%ilap=ib
           c_np(jj,idimnp(jj))%j=isi
           c_np(jj,idimnp(jj))%i1=i_1
           c_np(jj,idimnp(jj))%i2=i_2
           !write(58)ib,isi,i_2,i_1
           !write(588,*) idimnp(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ihnmn,ihnmx
          do i_1=ihnmn,ihnmx
           if(abs(levn(i_1)%j-levn(i_2)%j).le.2*isi.and.2*isi.le.(levn(i_1)%j+levn(i_2)%j)) then         
           idimnh(jj)=idimnh(jj)+1
           if(idimnh(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
           if(idimnh(jj).gt.ndro) stop
           c_nh(jj,idimnh(jj))%ilap=ib
           c_nh(jj,idimnh(jj))%j=isi
           c_nh(jj,idimnh(jj))%i1=i_1
           c_nh(jj,idimnh(jj))%i2=i_2
!           write(59)ib,isi,i_2,i_1
           !write(590,*) idimnh(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ippmn,ippmx
          do i_1=ippmn,ippmx
           if(abs(levp(i_1)%j-levp(i_2)%j).le.2*isi.and.2*isi.le.(levp(i_1)%j+levp(i_2)%j)) then
           idimpp(jj)=idimpp(jj)+1
          if(idimpp(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
          if(idimpp(jj).gt.ndro) stop
          c_pp(jj,idimpp(jj))%ilap=ib
          c_pp(jj,idimpp(jj))%j=isi
          c_pp(jj,idimpp(jj))%i1=i_1
          c_pp(jj,idimpp(jj))%i2=i_2
!           write(60)ib,isi,i_2,i_1
           !write(600,*) idimpp(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ihpmn,ihpmx
          do i_1=ihpmn,ihpmx
           if(abs(levp(i_1)%j-levp(i_2)%j).le.2*isi.and.2*isi.le.(levp(i_1)%j+levp(i_2)%j)) then         
           idimph(jj)=idimph(jj)+1
           if(idimph(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
           if(idimph(jj).gt.ndro) stop
           c_ph(jj,idimph(jj))%ilap=ib
           c_ph(jj,idimph(jj))%j=isi
           c_ph(jj,idimph(jj))%i1=i_1
           c_ph(jj,idimph(jj))%i2=i_2
           !write(61)ib,isi,i_2,i_1
           !write(610,*) idimph(jj),ib,isi,i_2,i_1
           end if
          end do
         end do
         
        end if!jj
       end do!isi
      end do !ib
     end do
        deallocate(levn,levp)
        call admat(nf,ipar,jcal,phonbs,idphon,no,phon1,phon2,mxt,xcc,ndcc,ipozx,ndbx,myid,numprocs,i_resh,c_pp,c_ph,c_np,c_nh,idimpp,idimph,idimnp,idimnh)

!       CALL execute_command_line('./run_admat.sh' )
!       stop


      close(99)
      close(12)
      close(22)
!      close(992)
    
      call MPI_FINALIZE(irc)
      end
