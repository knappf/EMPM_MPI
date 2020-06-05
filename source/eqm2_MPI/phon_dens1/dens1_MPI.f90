module dens1

use eofmod

contains 

subroutine read_cph(ifmx,cph_p,cph_n,ippmin,ippmax,ihpmin,ihpmax,ipnmin,ipnmax,ihnmin,ihnmax)

implicit none

include 'types_phon_dens1.inc'

type(amp_typ), dimension(:), allocatable :: camp

double precision, dimension(:,:,:), allocatable :: cph_p,cph_n

!character*9

integer :: ifmx,ippmin,ippmax,ihpmin,ihpmax,ipnmin,ipnmax,ihnmin,ihnmax
integer :: i,ijj,ilam,ndc,ndamp,ipar

ndamp=20000
allocate (camp(ndamp))

open(11,file='1phonon/1f_cp.dat',status='old',form='unformatted')

allocate (cph_p(ifmx,ippmin:ippmax,ihpmin:ihpmax))
cph_p=0.d0

ilam=0

!do while (.not.eof(11))
do ilam=1,ifmx 
!   ilam=ilam+1
     read(11)ipar,ijj,ndc
!       jphon(ilam)=ijj
        if (ndc.gt.ndamp) then
                           write(*,*)'small dimension of array camp'
                           stop
                   endif

read(11)(camp(i)%par,camp(i)%hol,camp(i)%am,i=1,ndc)

      do i=1,ndc
         cph_p(ilam,camp(i)%par,camp(i)%hol)=camp(i)%am
      enddo

enddo


close(11)

open(11,file='1phonon/1f_cn.dat',status='old',form='unformatted')

allocate (cph_n(ifmx,ipnmin:ipnmax,ihnmin:ihnmax))
cph_n=0.d0

ilam=0

!do while (.not.eof(11))
!   ilam=ilam+1
do ilam=1,ifmx
     read(11)ipar,ijj,ndc
!       jphon(ilam)=ijj
        if (ndc.gt.ndamp) then
                           write(*,*)'small dimension of array camp'
                           stop
                   endif
read(11)(camp(i)%par,camp(i)%hol,camp(i)%am,i=1,ndc)

      do i=1,ndc
         cph_n(ilam,camp(i)%par,camp(i)%hol)=camp(i)%am
      enddo

enddo


close(11)

deallocate(camp)

return
end subroutine read_cph

!************************************************************************
subroutine rop(ipart,ifmx,isimax,lev,xthres,jphon,cph_p,cph_n,myid,numprocs)

  use anglib   ! angular momentum staff
  use input_sp
  use eofmod

      implicit double precision(a-h,o-z)

!      include '/usr/lib/openmpi/include/mpif.h'
      include 'types_phon_dens1.inc'
      include 'input_phon_dens1.inc'

      double precision, dimension(:,:,:), allocatable, target :: cph_p,cph_n
      double precision, pointer,dimension(:,:,:) :: came

      integer, dimension(:), allocatable :: jphon,ius,ig_resh

      type(amp_typ), dimension(:), allocatable :: camp
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh

      integer ierr,myid,numprocs,irc

      integer :: n_seg

      character*10 namer
      character*9  namec
      character*4 nlam


!call MPI_INIT( ierr )
!--Who am I? --- get my rank=myid
!call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
!--How many processes in the global group?
!call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

if (myid.eq.0) then
  if (ipart.eq.1) write(*,*)'Calculation of 1-phonon particle neutron densities'
  if (ipart.eq.-1) write(*,*)'Calculation of 1-phonon particle proton densities'
  write(*,*)' # of MPI processes = ',numprocs
endif

      ndamp=10000
      ndrho=10000000     

      allocate (camp(ndamp))
!      allocate (jphon(ifmx))
!      jphon=0

      allocate(rh(ndrho))

      if (ipart.eq.1) then 
              namer='1f_rnp.dat'
              namec='1f_cn.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
      endif

      if (ipart.eq.-1) then 
              namer='1f_rpp.dat'
              namec='1f_cp.dat'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
      endif

!      allocate (came(ifmx,ipmin:ipmax,ihmin:ihmax))
!      came=0.d0

      allocate(ius(ifmx))

      ius=1  ! all by default 

      allocate(ig_resh(ifmx))
      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) then 
                 icount=icount+1
                 ig_resh(icount)=i
       endif
      enddo


      if (myid.eq.0) write(*,*)' Number of 0hom+1hom+2hom+3hom phonons ',icount
      
       call rperm(icount, ig_resh)


!      icount=0
!      do i=1,ifmx
!       if (ius(i).eq.1) then
!                 icount=icount+1
!                 ig_resh(icount)=i
!       endif
!      enddo

 
      if (ipart.eq.1) then 
!        if (myid.eq.0) write(*,*)'Calculation of 1-phonon particle neutron densities'
        nullify(came)
        came => cph_n
      endif 

      if (ipart.eq.-1) then
         nullify(came)
!        if (myid.eq.0) write(*,*)'Calculation of 1-phonon particle proton densities'
        came => cph_p
      endif 


!      open(33,file='1phonon/'//namer,status='unknown'
!     *,form='unformatted')


!      do ig=1,ifmx

!      n_seg=icount/numprocs+1   
      if (mod(icount,numprocs).eq.0) then 
              n_seg=icount/numprocs
      else
              n_seg=icount/numprocs+1
      endif 
!      n_rem=mod(icont,numprocs)

if (myid.eq.0) write(*,*) ' size of segment = ',n_seg

!      write(*,*)'Proces ',myid,'range:', myid*n_seg+1,(myid+1)*n_seg
      do igg=myid*n_seg+1,min((myid+1)*n_seg,icount) 

!      do igg=1,icount


       ig=ig_resh(igg)
  
!      write(*,*)' ig =', ig

      write(*,*) ' Process #',myid,'  calculating ig = ',ig

      write(nlam,'(i4.4)')ig

      open(33,file='scratch/'//namer//'_'//nlam,status='unknown',form='unformatted')


       iirg=0

!c       write(33)ig
!c       write(997,*)ig
       jig=jphon(ig)


!      if (ius(ig).ne.0) then 

      do ib=1,ifmx
       jib=jphon(ib)

       isi_min=iabs(jib-jig)
       isi_max=jib+jig

!      do isi=0,isimax
 
      do ip1=ipmin,ipmax
       jp1=lev(ip1)%j 
      do ip2=ipmin,ipmax
       jp2=lev(ip2)%j

!         if (((jp1+jp2)/2).lt.isi_max) isi_max=(jp1+jp2)/2
!         if ((iabs(jp1-jp2)/2).gt.isi_min) isi_min=iabs(jp1-jp2)/2
   
         do isi=isi_min,isi_max

         ronp=0.d0
      
          do ih=ihmin,ihmax
            jh=lev(ih)%j
            ifaz=(-1)**(jig+isi+(jh+jp1)/2)
            xsixj=sixj(2*jib,2*isi,2*jig,jp2,jh,jp1) 
            fact=(dfloat((2*jib+1)*(2*jig+1)*(2*isi+1)))**0.5d0
            ronp=ronp+dfloat(ifaz)*fact*came(ig,ip2,ih)*came(ib,ip1,ih)*xsixj
          enddo
             
       if (dabs(ronp).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in rop!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ip1
               rh(iirg)%i2=ip2
               rh(iirg)%rho=real(ronp)
       endif

        enddo ! loop ip2
      enddo ! loop ip1
      enddo ! loop isi
      enddo ! loop ib

!      endif

      write(33)ig,iirg
      if (iirg.gt.0) then
       write(33)(rh(iii)%ib,iii=1,iirg)
       write(33)(rh(iii)%isi,iii=1,iirg)
       write(33)(rh(iii)%i1,iii=1,iirg)
       write(33)(rh(iii)%i2,iii=1,iirg)
       write(33)(rh(iii)%rho,iii=1,iirg)
      endif

      if (iirg.eq.0) then       
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)float(iirg)
      endif

!      write(997,*)ig,iirg

!c      write(977,*)0,0,0,0,0.d0
!      endif

      close(33)
      enddo ! loop ig      
      

      deallocate(camp,rh)
!c      write(33)10000000
!c      write(33)0,0,0,0,0.d0

close(33)

!call MPI_FINALIZE(irc)

!nullify(came)
return
end subroutine  
!************************************************************************
subroutine roh(ipart,ifmx,isimax,lev,xthres,jphon,cph_p,cph_n,myid,numprocs)

   use anglib   ! angular momentum staff
   use input_sp    
   use eofmod

      implicit double precision(a-h,o-z)

      include 'types_phon_dens1.inc'
      include 'input_phon_dens1.inc'

      double precision, dimension(:,:,:), allocatable, target :: cph_p,cph_n
      double precision, pointer,dimension(:,:,:) :: came

      integer, dimension(:), allocatable :: jphon,ius,ig_resh

      type(amp_typ), dimension(:), allocatable :: camp
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh

      integer ierr,myid,numprocs,irc

      integer :: n_seg


      character*10 namer
      character*9  namec
      character*4 nlam

if (myid.eq.0) then
  if (ipart.eq.1) write(*,*)'Calculation of 1-phonon hole neutron densities'
  if (ipart.eq.-1) write(*,*)'Calculation of 1-phonon hole proton densities'
  write(*,*)' # of MPI processes = ',numprocs
endif



      ndamp=10000    
      ndrho=10000000  
      allocate (camp(ndamp))
!      allocate (jphon(ifmx))
!      jphon=0

      allocate(rh(ndrho))


      if (ipart.eq.1) then 
              namer='1f_rnh.dat'
!c              namerf='1f_rnhf.dat'
              namec='1f_cn.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
      endif

      if (ipart.eq.-1) then 
              namer='1f_rph.dat'
!c              namerf='1f_rphf.dat'
              namec='1f_cp.dat'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
      endif

!      allocate (came(ifmx,ipmin:ipmax,ihmin:ihmax))
!      came=0.d0

      allocate(ius(ifmx))

      ius=1  ! all by default

      allocate(ig_resh(ifmx))
      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) then
                 icount=icount+1
                 ig_resh(icount)=i
       endif
      enddo

      if (myid.eq.0) write(*,*)' Number of 0hom+1hom+2hom+3hom phonons ',icount

      call rperm(icount, ig_resh)


      if (ipart.eq.1) then
!        write(*,*)'Calculation of 1-phonon particle neutron densities'
        nullify(came)
        came => cph_n
      endif

      if (ipart.eq.-1) then
         nullify(came)
!        write(*,*)'Calculation of 1-phonon particle proton densities'
        came => cph_p
      endif


      if (mod(icount,numprocs).eq.0) then
              n_seg=icount/numprocs
      else
              n_seg=icount/numprocs+1
      endif
!      n_rem=mod(icont,numprocs)

     if (myid.eq.0) write(*,*) ' size of segment = ',n_seg


     do igg=myid*n_seg+1,min((myid+1)*n_seg,icount)
      
!      do igg=1,ifmx

      ig=ig_resh(igg)
!      ig=igg
      write(*,*) ' Process #',myid,'  calculating ig = ',ig

!      write(*,*)' ig =', ig


      write(nlam,'(i4.4)')ig

      open(33,file='scratch/'//namer//'_'//nlam,status='unknown',form='unformatted')


       iirg=0

!c       write(33)ig
!c       write(433)ig
!c       write(998,*)ig
       jig=jphon(ig)

!       if (ius(ig).ne.0) then


      do ib=1,ifmx
       jib=jphon(ib)
       isi_min=iabs(jib-jig)
       isi_max=jib+jig


!      do isi=0,isimax
      do ih1=ihmin,ihmax
       jh1=lev(ih1)%j 
      do ih2=ihmin,ihmax
       jh2=lev(ih2)%j
 
       do isi=isi_min,isi_max

          ronh=0.d0
      
          do ip=ipmin,ipmax
            jp=lev(ip)%j
            ifaz=(-1)**(jib+(jp+jh2)/2)
            xsixj=sixj(2*jib,2*isi,2*jig,jh1,jp,jh2) 
            fact=(dfloat((2*jib+1)*(2*jig+1)*(2*isi+1)))**0.5d0
            ronh=ronh-dfloat(ifaz)*fact*came(ig,ip,ih1)*came(ib,ip,ih2)*xsixj
          enddo
          
!c       if (dabs(ronh).gt.xrotrunc) write(33)ib,isi,ih1,ih2,ronh

       if (dabs(ronh).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in roh!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ih1
               rh(iirg)%i2=ih2
               rh(iirg)%rho=real(ronh)
!c        write(33)ib,isi,ip1,ip2,ronp
       endif


!c       if (ib.eq.ig.and.ih1.eq.ih2.and.isi.eq.0) 
!c     *ronh=ronh+dfloat((jh1+1)*(2*jig+1))**0.5d0
!c       if (dabs(ronh).gt.xrotrunc) write(433)ib,isi,ih1,ih2,ronh
!c       if (dabs(ronh).gt.xrotrunc) 
!c     *write(998,*)ib,isi,ih1,ih2,ronh

!  55    continue      
        enddo ! loop ih2
      enddo ! loop ih1
      enddo ! loop isi
      enddo ! loop ib

!      endif

      write(33)ig,iirg
      if (iirg.gt.0) then
       write(33)(rh(iii)%ib,iii=1,iirg)
       write(33)(rh(iii)%isi,iii=1,iirg)
       write(33)(rh(iii)%i1,iii=1,iirg)
       write(33)(rh(iii)%i2,iii=1,iirg)
       write(33)(rh(iii)%rho,iii=1,iirg)
      endif

      if (iirg.eq.0) then
       write(33)iirg
       write(33)iirg
       write(33)iirg
       write(33)iirg
       write(33)float(iirg)
      endif

!      endif
!      write(998,*)ig,iirg

!c      write(33)0,0,0,0,0.d0
!c      write(433)0,0,0,0,0.d0
!c      write(998,*)0,0,0,0,0.d0

      close(33)
      enddo ! loop ig      

      
      deallocate(camp,rh)
!c      write(33)10000000
!c      write(33)0,0,0,0,0.d0
!      close(33)
!c      write(433)10000000
!c      write(433)0,0,0,0,0.d0
!c      close(433)

!nullify(came)
return      
end subroutine  
!
subroutine rperm(N, p)

 integer(kind=4), intent(in) :: N
 integer(kind=4), dimension(:), intent(out) :: p

 integer(kind=4) :: i
 integer(kind=4) :: k, j, ipj, itemp, m
 real(kind=4), dimension(100) :: u

p = (/ (i, i=1,N) /)

! Generate up to 100 U(0,1) numbers at a time.
do i=1,N,100
m = min(N-i+1, 100)
call random_number(u)
do j=1,m
ipj = i+j-1
k = int(u(j)*(N-ipj+1)) + ipj
itemp = p(ipj)
p(ipj) = p(k)
p(k) = itemp
end do
end do
return

end subroutine rperm


!************************************************************************
end module dens1      
      
