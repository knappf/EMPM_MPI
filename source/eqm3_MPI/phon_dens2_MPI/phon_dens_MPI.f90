!*     Program phon_dens computes 2-phonon densities in proton-neutron 
!*     J-coupled formalism.
!     
!*     last update 7.7.2020

      program phon_dens

      use rdens
      use input_sp
!c      use rdens_al

      implicit double precision (a-h,o-z)

!      include 'input_phon_dens.inc'

      include 'mpif.h'
      include 'formats_phon_dens.inc'
!      include 'types_phon_dens.inc'

      type(level_typ),dimension(:), allocatable :: levn,levp

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,iphous2

      character(len=16) fname
      character*1 :: ch1
      character*2 :: ch2

      integer :: myid,ierr,numprocs
      integer,dimension(:),allocatable :: idimnp,idimnh,idimph,idimpp
      
!*     loading of input data 

call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )



      xrotrunc=1.d-10
      
!      write(*,*)'Loading of input '

!      open(63,file='2_phon_dens_calc.dat',status='old',form='formatted',position='append')

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
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


      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
!c       jphon(i)=ijj
      enddo

      close(3)


      ifmx=i 
      jmax=jmaxn

      jamax=2*jmax
      jbmax=jamax

      isimax=jmax


      call loadphon(phonbs,nphon)

      if (myid.eq.0) then 
       write(*,*)' Number of 1phonon states ',nphon(1)
       write(*,*)' Number of 2phonon states ',nphon(2)
      endif 

      nff=2

      call selphon(nff,phonbs,nphon,iphous,ns1,myid)
      call selphon2(nff,phonbs,nphon,iphous2,ns2,myid)

!'*******************Modifications-1/05/2020*******************************************************

      jmax2=0
      do ia=1,nphon(2)
        if(phonbs(2,ia)%jj.gt.jmax2) jmax2=phonbs(2,ia)%jj
      end do
      if (myid.eq.0) write(*,*) 'Jmax 2ph=',jmax2
!****************************************************************************************************
!**************idim:arrays containing the number of combinations of ib isi i_2 i_1 for each Ja*******
!**************Nedeed for the check at the end of the 2ph densities calculation**********************
!****************************************************************************************************
      allocate(idimnp(0:jmax2),idimnh(0:jmax2),idimpp(0:jmax2),idimph(0:jmax2))
      idimnp=0
      idimnh=0
      idimph=0
      idimpp=0
516   format(i2.2)
      do jj=0,jmax2
      write(ch2,516)jj

!open(58,file='scratch/2f_rnp_comb.dat_'//ch2,status='unknown',form='unformatted')
!open(59,file='scratch/2f_rnh_comb.dat_'//ch2,status='unknown',form='unformatted')
!open(60,file='scratch/2f_rpp_comb.dat_'//ch2,status='unknown',form='unformatted')
!open(61,file='scratch/2f_rph_comb.dat_'//ch2,status='unknown',form='unformatted')
      write(588,*) jj
      write(599,*) jj
      write(600,*) jj
      write(610,*) jj
      do ib=1,nphon(2)
      
      jb=phonbs(2,ib)%jj
       do isi=0,isimax
        if(abs(jb-jj).le.isi.and.isi.le.(jb+jj))then
         !write(19,*) !iaaa,ib,isi,ja,jb!ib,iaaa,isi,ja,jb
         do i_2=ipnmn,ipnmx
          do i_1=ipnmn,ipnmx
           if(iabs(levn(i_1)%j-levn(i_2)%j).le.2*isi.and.2*isi.le.(levn(i_1)%j+levn(i_2)%j)) then
           idimnp(jj)=idimnp(jj)+1
!           write(58)ib,isi,i_2,i_1
!           write(588,*) idimnp(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ihnmn,ihnmx
          do i_1=ihnmn,ihnmx
           if(iabs(levn(i_1)%j-levn(i_2)%j).le.2*isi.and.2*isi.le.(levn(i_1)%j+levn(i_2)%j)) then         
           idimnh(jj)=idimnh(jj)+1
!           write(59)ib,isi,i_2,i_1
!           write(590,*) idimnh(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ippmn,ippmx
          do i_1=ippmn,ippmx
           if(iabs(levp(i_1)%j-levp(i_2)%j).le.2*isi.and.2*isi.le.(levp(i_1)%j+levp(i_2)%j)) then
           idimpp(jj)=idimpp(jj)+1
!           write(60)ib,isi,i_2,i_1
!           write(600,*) idimpp(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ihpmn,ihpmx
          do i_1=ihpmn,ihpmx
           if(iabs(levp(i_1)%j-levp(i_2)%j).le.2*isi.and.2*isi.le.(levp(i_1)%j+levp(i_2)%j)) then         
           idimph(jj)=idimph(jj)+1
!           write(61)ib,isi,i_2,i_1
!           write(610,*) idimph(jj),ib,isi,i_2,i_1
           end if
          end do
         end do
         
        end if!jj
       end do!isi
      end do !ib
!      close(58)
!      close(59)
!      close(60)
!      close(61)

     end do


call roo(jamax,jbmax,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax,iphous,iphous2,phonbs,nphon,ns1,ns2,myid,numprocs,levn,levp,idimnp,idimnh,idimpp,idimph)

call MPI_FINALIZE(irc)


      end 
!*     
!*     END of the main program 
!* 

      




      
      
