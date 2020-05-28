!   lat update 11.6.2013
!***********************************************  

      module fullmat

      contains

      subroutine bstot(nfon,ipar,jcal,nbs,ndmpho,phonbs)

      implicit double precision(a-h,o-z)

      include 'types_fullham.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension(:), allocatable :: nphon,ndmpho
      integer, dimension(:,:), allocatable :: nbs

      nffmx=500000

      if (.not.allocated(phonbs)) then
            allocate(phonbs(0:4,nffmx))
      endif 
      if (.not.allocated(nphon)) then
         allocate(nphon(0:4))
         nphon=0
      endif

      if (.not.allocated(ndmpho)) then      
         allocate(ndmpho(0:4))
         ndmpho=0
      endif

      if (.not.allocated(nbs)) then       
         allocate(nbs(0:4,nffmx))
         nbs=0
      endif

      if (nfon.eq.0) then
              phonbs(0,1)%par=1
              phonbs(0,1)%jj=0
              phonbs(0,1)%enf=0.d0
              nphon(0)=1
      endif


      if (nfon.ne.0) then 

        if (nfon.eq.1) namef='1phonon/1f_states.dat'
        if (nfon.eq.2) namef='2phonon/2f_states.dat'
        if (nfon.eq.3) namef='3phonon/3f_states.dat'

      
      open(1,file=namef,status='old',form='unformatted')

      do while (.not.eof(1))
       read(1)i,ipartt,ijj,en
         phonbs(nfon,i)%par=ipartt
         phonbs(nfon,i)%jj=ijj
         phonbs(nfon,i)%enf=en
      enddo
      nphon(nfon)=i

      close(1)
      
      endif

      ii=0

      do i=1,nphon(nfon)
      ipart=phonbs(nfon,i)%par
      jjt=phonbs(nfon,i)%jj
      ene=phonbs(nfon,i)%enf

      if (ene.lt.0.d0) goto 661

      if (ipart.eq.ipar.and.jjt.eq.jcal) then 
              ii=ii+1
              nbs(nfon,ii)=i
      endif

 661  continue

      enddo

      ndmpho(nfon)=ii
      
      return
      end subroutine bstot

!************************************************************************
      subroutine nondiag1r(nff,nbs,ndmpho,phonbs,ipcal,jcal)

      implicit double precision (a-h,o-z)

      include 'types_fullham.inc'

      type (phon_typ), dimension (:,:), allocatable :: phonbs
 
      double precision, dimension(:,:,:), allocatable :: 
     *vred

      double precision, dimension (:), allocatable :: vint,vintm

      integer, dimension(:,:), allocatable :: nbs
      integer, dimension(:), allocatable :: ipoza,ipozinv
      type(ampr_typ),dimension(:,:),allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc,ndmpho

      character*30 fname,fnamev,fnamevm,fnameo



      if (nff.eq.2) then 
              fnamev='Vint_phon12.dat'
              fnamevm='Vint_phon01.dat'
              fnameo='nondiag12.dat'
      endif

      if (nff.eq.3) then 
              fnamev='Vint_phon23.dat'
              fnamevm='Vint_phon01.dat'
              fnameo='nondiag23.dat'
      endif


      nfona=5000  ! nff-1 phonon dimension   ! redukovana dimenzia nemusi obsahovat vsetky 1 fononove stavy
      nfong=50000  ! nff-1 phonon dimension
      nfons=200   ! 1 phonon dimension sigma index
      ndxcc=10000

      open(12,file=fnameo,status='unknown',form='unformatted')

      open(3,file=fnamev,status='old',form='unformatted')
!      open(4,file=fnamevm,status='old',form='unformatted')

!      allocate(vintm(nfons))
!      vintm=0.d0
!      do while (.not.eof(4))
!        read(4)isi,vintm(isi)
!      enddo

  
      allocate (ipoza(nfong),ipozinv(nfong))
      ipoza=0
      ipozinv=0

      write(*,*)'aaaa'

      allocate(vred(nfong,nfons,nfona))

      write(*,*)'bbb'

      vred=0.d0

 
      iii=0

!    POZOR na nezrovnalosti s idexami ig ia

      do while (.not.eof(3))
        read(3)ig,ia,is,vr


!c        write(*,*)ia,ig,is,vr

        if (ipoza(ig).eq.0) then
                iii=iii+1
!c                write(*,*)ipozic,ia
                ipoza(ig)=iii
                ipozinv(iii)=ig
!c                write(*,*)ipozz,ia
        endif

!c        write(9999,*)ig,il,ib,vr
!c        ili=phonuspov(il)
!c        igi=phonuspov(ig)
!c        ibi=phonuspov(ib)
!c        if (igi.eq.0.or.ili.eq.0.or.ibi.eq.0) goto 11

      vred(ia,is,ipoza(ig))=vr
!      vred(ia,ig,is)=vr

 11   continue
      enddo

      close(3)

      write(*,*)' Ipozzmax =', iii
c      stop


      if (nff.eq.2) then 
              fname='2phonon/2f_x.dat'
              ndimj=1610000
                    endif

      if (nff.eq.3) then
              fname='3phonon/3f_x.dat'
              ndimj=1800000
                    endif

      ndimi=2000
      call readx(fname,xcc,ndcc,ipcal,jcal,ndimj)

      nf2=ndmpho(nff)
      nf1=ndmpho(nff-1)

c      do i=1,nf2 ! cycle over alpha
c      write(*,*)i,nf2
c        ii=nbs(nff,i)     

        allocate(vint(nf2))
        vint=0.d0

        do j=1,nf1

         write(*,*)j,nf1

          vint=0.d0
          jj=nbs(nff-1,j)
c         jjj=ipoza(jj)
c          ib=ipbaset(1,j)  
c         ndmx=ndcc(ii)
         ndmx=ndcc(nbs(nff,1))

          do kk=1,ndmx
           do i=1,nf2
            
           ii=nbs(nff,i)
           isi=xcc(ii,kk)%is
           igi=xcc(ii,kk)%ig

            vint(i)=vint(i)+xcc(ii,kk)%am*
     *vred(igi,isi,ipoza(jj))
!          if (igi.eq.jj) then 
!            vint(i)=vint(i)+xcc(ii,kk)%am*
!     *vintm(isi)       
!          endif


            enddo
          enddo

          do i=1,nf2     
            if (dabs(vint(i)).gt.1.d-10) write(12)i,j,vint(i)
          enddo
c          if (dabs(vint).gt.1.d-10) write(998,*)i,j,vint
        enddo

c      enddo

      deallocate(vred,xcc)

      close(2)
      close(3)

      return
      end subroutine nondiag1r
!************************************************************************

      subroutine nondiag2(nff,nbs,ndmpho,phonbs,ipcal,jcal)

      implicit double precision (a-h,o-z)

      include 'types_fullham.inc'

      type (phon_typ), dimension (:,:), allocatable :: phonbs
 
      double precision, dimension(:,:), allocatable :: 
     *vred
      double precision, dimension(:), allocatable ::
     *vinti

      integer, dimension(:,:), allocatable :: nbs
      type(ampr_typ),dimension(:,:),allocatable :: xcc,xaa
      integer, dimension (:), allocatable :: ndcc,ndaa,ndmpho

      character*30 fname,fname2,fname3


      if (nff.eq.2) then 

      open(12,file='nondiag02.dat',status='unknown',form='unformatted')

      open(3,file='Vint_phon2.dat',status='old',form='unformatted')

      nfon1=1600
      nfon1t=1600
      ndxcc=100000

      allocate(vred(nfon1,nfon1t))
      vred=0.d0

      do while (.not.eof(3))
        read(3)ig,is,vr
c        write(9999,*)ig,il,ib,vr
c        ili=phonuspov(il)
c        igi=phonuspov(ig)
c        ibi=phonuspov(ib)
c        if (igi.eq.0.or.ili.eq.0.or.ibi.eq.0) goto 11
        vred(ig,is)=vr
 11   continue
      enddo

      close(3)


      fname='2phonon/2f_x.dat'

      ndimi=2000
      ndimj=16000
      call readx(fname,xcc,ndcc,ipcal,jcal,ndimj)


      nf2=ndmpho(nff)
      nf0=ndmpho(nff-2)

      do i=1,nf2 ! cycle over alpha

        ii=nbs(nff,i)         
        ndmx=ndcc(ii)
        jibi=phonbs(nff,ii)%jj
        xfact=dfloat(2*jibi+1)**(-1.d0)


        do j=1,nf0
          vint=0.d0
          jj=nbs(nff-2,j)


c          ib=ipbaset(1,j)  

          do kk=1,ndmx
           isi=xcc(ii,kk)%is
           igi=xcc(ii,kk)%ig

           jisi=phonbs(nff-1,isi)%jj
           jigi=phonbs(nff-1,igi)%jj
           
           ifaz=(-1)**(jisi+jigi+jibi)

            vint=vint+xcc(ii,kk)%am*
     *vred(igi,isi)*dfloat(ifaz)*(dfloat(2*jisi+1)**0.5d0)*xfact 
          enddo
      
          if (dabs(vint).gt.1.d-10) write(12)i,j,vint
c          if (dabs(vint).gt.1.d-10) write(998,*)i,j,vint
        enddo

      enddo

      deallocate(vred,xcc)

      close(2)
      close(3)

      endif

      if (nff.eq.3) then 

      open(12,file='nondiag13.dat',status='unknown',form='unformatted')

      open(3,file='Vint_phon2.dat',status='old',form='unformatted')

      nfon1=200
      nfon1t=200
      ndxcc=100000

      allocate(vred(nfon1,nfon1t))
      vred=0.d0

      do while (.not.eof(3))
        read(3)ig,is,vr
c        write(9999,*)ig,il,ib,vr
c        ili=phonuspov(il)
c        igi=phonuspov(ig)
c        ibi=phonuspov(ib)
c        if (igi.eq.0.or.ili.eq.0.or.ibi.eq.0) goto 11
        vred(ig,is)=vr
c 11   continue
      enddo

      close(3)


      fname2='2phonon/2f_x.dat'
      fname3='3phonon/3f_x.dat'


      ndimi=2000
      ndimj=2000

c      write(*,*) 'bbbb '

      ndimj=500000
      call readx(fname3,xcc,ndcc,ipcal,jcal,ndimj)



      nf2=ndmpho(nff)
      nf0=ndmpho(nff-2)

      allocate(vinti(nf2))
      vinti=0.d0

      do i=1,1 !nf2 ! cycle over alpha

      write(*,*) 'i= ',i,nf2

        ii=nbs(nff,i)         
        ndmx=ndcc(ii)
        jibi=phonbs(nff,ii)%jj
        xfact=dfloat(2*jibi+1)**(-1.d0)


        do j=1,nf0
          
          write(998,*) 'j= ',j,nf0
          vint=0.d0
          vinti=0.d0
          jj=nbs(nff-2,j)

          jold=-1000
          ipold=-10

          do kk=1,ndmx
c           write(*,*)kk,ndmx
           isi=xcc(ii,kk)%is
           igi=xcc(ii,kk)%ig

           jisi=phonbs(1,isi)%jj
           jigi=phonbs(nff-1,igi)%jj
           ipg=phonbs(nff-1,igi)%par

           if ((jigi.ne.jold).or.(ipg.ne.ipold)) then
            ndimj=15000
            call readx(fname2,xaa,ndaa,ipg,jigi,ndimj)
             jold=jigi
             ipold=ipg
           endif

           
           ifaz=(-1)**(jisi+jigi+jibi)
         
           ndma=ndaa(igi)

           do ll=1,ndma

           iaaa=xaa(igi,ll)%ig
           isii=xaa(igi,ll)%is

           if (iaaa.eq.jj) then

            do iii=1,nf2       
            iiii=nbs(nff,iii)

c           vint=vint+xcc(ii,kk)%am*xaa(igi,ll)%am*
c    *vred(isi,isii)*dfloat(ifaz)*xfact 
            vinti(iii)=vinti(iii)+xcc(iiii,kk)%am*xaa(igi,ll)%am*
     *vred(isi,isii)*dfloat(ifaz)*xfact

            enddo
c            if (dabs(vint).gt.1.d-10) write(917,*)i,j,vint

            endif

            enddo

     
          enddo
      
          do iii=1,nf2
          if (dabs(vinti(iii)).gt.1.d-10) then 
c                  write(*,*)vint
                  write(12)iii,j,vinti(iii)
          endif
          enddo
c          if (dabs(vint).gt.1.d-10) write(998,*)i,j,vint
        enddo

      enddo

      deallocate(vred,xcc)

      close(2)
      close(3)

      endif





      return
      end subroutine nondiag2
!************************************************************************

      subroutine nondiag0(nff,nbs,ndmpho,phonbs,ipcal,jcal)

        implicit double precision (a-h,o-z)
  
        include 'types_fullham.inc'
  
        type (phon_typ), dimension (:,:), allocatable :: phonbs
   
        double precision, dimension(:), allocatable :: 
     *vred
        double precision, dimension(:), allocatable ::
     *vinti
  
        integer, dimension(:,:), allocatable :: nbs
        type(ampr_typ),dimension(:,:),allocatable :: xcc,xaa
        integer, dimension (:), allocatable :: ndcc,ndaa,ndmpho
  
        character*30 fname,fname2,fname3
  
  
        if (nff.eq.1) then 
  
        open(12,file='nondiag01.dat',status='unknown'
     *,form='unformatted')
  
        open(3,file='Vint_phon01.dat',status='old',form='unformatted')
  
        nfon1=1600
      
        allocate(vred(nfon1))
        vred=0.d0
  
        do while (.not.eof(3))
          read(3)is,vr
          vred(is)=vr
         continue
        enddo
  
        close(3)
  
        nf1=ndmpho(nff)
  
        do i=1,nf1 ! cycle over alpha     
          ii=nbs(nff,i)                   
          vint=vred(ii)
            if (dabs(vint).gt.1.d-10) write(12)i,1,vint
        enddo
  
      
  
        deallocate(vred)
  
        close(2)
        close(3)
  
        endif

  
        return

        end subroutine nondiag0

!c************************************************************************

      subroutine readx(fname,xccr,ndcc,ipcal,jcal,ndimj)

      implicit double precision (a-h,o-z)

c      include 'formats_phon_int.inc'
      include 'types_fullham.inc'

      type (amp_typ), dimension(:,:), allocatable :: xcc
      type (ampr_typ), dimension(:,:), allocatable :: xccr

      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

c      allocate(xcc(ndimi,ndimj))
c      allocate(ndcc(ndimi))

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

c      ndimj=533000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      ndimj=idphon

c      write(*,*)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc)
      if (allocated(xccr).eq..TRUE.) deallocate(xccr)
      if (allocated(ndcc).eq..TRUE.) deallocate(ndcc)

c      write(*,*)'ccc'
      allocate(xcc(1,idphon))
      allocate(xccr(ilamp:ilamp+no,idphon))
      allocate(ndcc(ilamp:ilamp+no))
c      write(*,*)'ddd'

c      stop


      do ilam=1,no

c      write(*,*)ilam,ilamp

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  

       read(2)(xcc(1,i)%ig,xcc(1,i)%is
     *,xcc(1,i)%am,i=1,idphon)
       xccr(ilam+ilamp,:)%ig=xcc(1,:)%ig
       xccr(ilam+ilamp,:)%is=xcc(1,:)%is
       xccr(ilam+ilamp,:)%am=real(xcc(1,:)%am)

       ndcc(ilam+ilamp)=idphon
c       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is
c     *,xcc(ilam+ilamp,i)%am,i=1,idphon)
c       ndcc(ilam+ilamp)=idphon
      enddo

      ilamp=ilamp+no

      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11

      deallocate(xcc,ndcc)
      deallocate(xccr)

      enddo

  11  continue     

      close(2)

      return
      end subroutine readx
****************************************
!c************************************************************************

  subroutine read_ampl(fname,xccr,ndcc,ipcal,jcal,ndimj)

  implicit double precision (a-h,o-z)

c      include 'formats_phon_int.inc'
  include 'types_fullham.inc'

  type (amp_typ), dimension(:,:), allocatable :: xcc
  type (ampr_typ), dimension(:,:), allocatable :: xccr

  integer, dimension (:), allocatable :: ndcc

  character(len=30)fname

c      allocate(xcc(ndimi,ndimj))
c      allocate(ndcc(ndimi))

  open(2,file=fname,status='old',form='unformatted')

  ilamp=0

c      ndimj=533000

 

  do while (.not.eof(2))

  read(2)ipar,ijj,no,idphon

  ndimj=idphon

c      write(*,*)ipar,ijj,no,idphon

  if (allocated(xcc).eq..TRUE.) deallocate(xcc)
  if (allocated(xccr).eq..TRUE.) deallocate(xccr)
  if (allocated(ndcc).eq..TRUE.) deallocate(ndcc)

c      write(*,*)'ccc'
  allocate(xcc(1,idphon))
  allocate(xccr(ilamp:ilamp+no,idphon))
  allocate(ndcc(ilamp:ilamp+no))
c      write(*,*)'ddd'

c      stop


  do ilam=1,no

c      write(*,*)ilam,ilamp

   
     if (idphon.gt.ndimj) then 
  write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
           stop
       endif


   read(2)(xcc(1,i)%ig,xcc(1,i)%is
 *,xcc(1,i)%am,i=1,idphon)
   xccr(ilam+ilamp,:)%ig=xcc(1,:)%ig
   xccr(ilam+ilamp,:)%is=xcc(1,:)%is
   xccr(ilam+ilamp,:)%am=real(xcc(1,:)%am)

   ndcc(ilam+ilamp)=idphon
c       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is
c     *,xcc(ilam+ilamp,i)%am,i=1,idphon)
c       ndcc(ilam+ilamp)=idphon
  enddo

  ilamp=ilamp+no

  if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11

  deallocate(xcc,ndcc)
  deallocate(xccr)

  enddo

11  continue     

  close(2)

  return
  end subroutine read_ampl
****************************************

      subroutine fulhamt(ndmpho,nbs)

      implicit double precision (a-h,o-z)

      double precision, dimension(:,:), allocatable :: amtr
      double precision, dimension(:), allocatable :: wr,work

      integer, dimension(:), allocatable :: ndmpho
      integer,  dimension (:,:), allocatable :: nbs
      character*30 :: nameen,namewf,namecon


c     construction of total H matrix


      ndimtot=ndmpho(0)+ndmpho(1)+ndmpho(2)+ndmpho(3)
     
      allocate(amtr(ndimtot,ndimtot))
      amtr=0.d0


      ndwr=300000
      allocate(wr(ndwr))

      if (ndmpho(0).gt.0) amtr(1,1)=0.d0

      if (ndmpho(1).gt.0) then 

      open(2,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')
     
      do while (.not.eof(2))
        read(2)nttt,iptt,jttt,enttt 
          if (nttt.gt.ndwr) then
                  write(*,*)'1f_states.dat: Small dimension of array wr'
                 stop
          endif
        wr(nttt)=enttt 
      enddo

      close(2)

    
      do i=1,ndmpho(1)
        iii=nbs(1,i)
        amtr(i+ndmpho(0),i+ndmpho(0))=wr(iii)
      enddo
      endif


c      close(2)

       if (ndmpho(2).gt.0) then

      open(2,file='2phonon/2f_states.dat',status='old',
     *form='unformatted')
     
      do while (.not.eof(2))
        read(2)nttt,iptt,jttt,enttt 
          if (nttt.gt.ndwr) then
                  write(*,*)'2f_states.dat: Small dimension of array wr'
                 stop
          endif
        wr(nttt)=enttt 
c        write(9995,*)wr(nttt)
      enddo

      close(2)

      write(*,*)'OK1'

      do i=1,ndmpho(2)
        iii=nbs(2,i)
        amtr(i+ndmpho(1)+ndmpho(0),i+ndmpho(1)+ndmpho(0))=wr(iii)
      enddo

      write(*,*)'OK2'

      endif

      if (ndmpho(3).gt.0) then

      open(2,file='3phonon/3f_states.dat',status='old',
     *form='unformatted')
     
      do while (.not.eof(2))
        read(2)nttt,iptt,jttt,enttt 
          if (nttt.gt.ndwr) then
                  write(*,*)'3f_states.dat: Small dimension of array wr'
                 stop
          endif
        wr(nttt)=enttt 
c        write(9995,*)wr(nttt)
      enddo

      close(2)

      write(*,*)'OK3'

      do i=1,ndmpho(3)
        iii=nbs(3,i)
        amtr(i+ndmpho(2)+ndmpho(1)+ndmpho(0),
     *ndmpho(2)+i+ndmpho(1)+ndmpho(0))=wr(iii)
      enddo

      write(*,*)'OK4'

      endif


 
      if (ndmpho(2).gt.0.and.ndmpho(1).gt.0) then         
      open(3,file='nondiag12.dat',status='old',form='unformatted')


      do while (.not.eof(3))
        read(3)ii,jj,vv
        amtr(ii+ndmpho(1)+ndmpho(0),jj+ndmpho(0))=vv
      enddo 

      close(3)
      endif

      if (ndmpho(3).gt.0.and.ndmpho(3).gt.0) then
         open(3,file='nondiag23.dat',status='old',form='unformatted')


       do while (.not.eof(3))
        read(3)ii,jj,vv
        amtr(ii+ndmpho(2)+ndmpho(1)+ndmpho(0),jj+ndmpho(1)+ndmpho(0))=vv
       enddo

       close(3)
       endif


!       if (ndmpho(0).gt.0.and.ndmpho(1).gt.0) then 
!        open(3,file='nondiag01.dat',status='old',form='unformatted')
  
  
!        do while (.not.eof(3))
!          read(3)ii,jj,vv
!          amtr(ii+ndmpho(0),jj)=vv
!        enddo 
  
!        close(3)
!        endif



      if (ndmpho(0).gt.0.and.ndmpho(2).gt.0) then 
      open(3,file='nondiag02.dat',status='old',form='unformatted')


      do while (.not.eof(3))
        read(3)ii,jj,vv
        amtr(ii+ndmpho(1)+ndmpho(0),jj)=vv
      enddo 

      close(3)
      endif

      if (ndmpho(1).gt.0.and.ndmpho(3).gt.0) then 
      open(3,file='nondiag13.dat',status='old',form='unformatted')


      do while (.not.eof(3))
        read(3)ii,jj,vv
        amtr(ii+ndmpho(2)+ndmpho(1)+ndmpho(0),jj+ndmpho(0))=vv
      enddo 

      close(3)
      endif



      deallocate(wr)

      iout=0

      if (iout.eq.1) then 

  302 format(10000f13.4)
      do i=1,ndimtot
        write(203,302)(amtr(i,j),j=1,ndimtot)
      enddo

      endif


      allocate(wr(ndimtot))

      wr=0.d0
      lwork=26*ndimtot
      allocate(work(lwork))
      work=0.d0

      call DSYEV('V','L',ndimtot,amtr,ndimtot,wr,WORK,LWORK,INFO) 

      write(*,*)' info =', info

      open(23,file='eigenenergies.dat',status='unknown',
     *form='formatted')
      do i=1,ndimtot
        write(23,*)wr(i)
      enddo

      close(23)


      iphomax=3
      if (ndmpho(3).eq.0) iphomax=2
      if (ndmpho(2).eq.0.and.ndmpho(3).eq.0) iphomax=1
      if (ndmpho(1).eq.0.and.ndmpho(2).eq.0.and.ndmpho(3).eq.0) 
     *iphomax=0

      write(*,*)' Phonmax ',iphomax

      if (iphomax.eq.2) then
              nameen='ef_2f.dat'
              namewf='wf_2f.dat'
              namecon='eigst_cont_2f.dat'
      endif

      if (iphomax.eq.1) then
             nameen='ef_1f.dat'
             namewf='wf_1f.dat'
             namecon='eigst_cont_1f.dat'
      endif

      if (iphomax.eq.3) then
        nameen='ef_3f.dat'
        namewf='wf_3f.dat'
        namecon='eigst_cont_3f.dat'
      endif



      open(23,file=namewf,status='unknown',
     *form='unformatted')
      write(23)iphomax
      write(23)(ndmpho(ipho),ipho=0,iphomax)
      write(23)ndimtot
      write(23)((nbs(ipho,i),i=1,ndmpho(ipho)),ipho=0,iphomax)
      write(23)((amtr(i,j),i=1,ndimtot),j=1,ndimtot)


      write(201,302)(wr(i),i=1,ndimtot)
      write(201,*)
      do i=1,ndimtot
         write(201,302)(amtr(i,j),j=1,10)
      enddo


      close(23)


      
      open(23,file=nameen,status='unknown',
     *form='unformatted')
      write(23)ndimtot
      write(23)(wr(i),i=1,ndimtot)
      close(23)


 274  format(i5,6f12.5) 

      open(2,file=namecon,status='unknown',form='formatted')
      do i=1,ndimtot
        con0=0.d0
        con1=0.d0
        con2=0.d0
        con3=0.d0

        do k=0,iphomax

        aaa=0.d0
         
         if (k.eq.0) then
              jj=0
         else
                 
           jj=0  
           do kkk=0,k-1       
            jj=jj+ndmpho(kkk)
           enddo
         endif

         do j=1+jj,ndmpho(k)+jj
           eee=amtr(j,i)
           aaa=aaa+amtr(j,i)**2.0d0
         enddo
         
         if (k.eq.0) con0=aaa
         if (k.eq.1) con1=aaa
         if (k.eq.2) con2=aaa
         if (k.eq.3) con3=aaa

         enddo 

        write(2,274)i,wr(i),con0,con1,con2,con3 
      enddo

      close(2)

      deallocate(amtr,wr,work)


c      do infon=0,2
c      write(2005,*)'phonon subspace ',infon
c      do i=1,ndmpho(infon)
c      write(2005,*)i,ipbaset(infon,i)
c      enddo
c      enddo

c      do i=1,ndimtot
c        write(2006,*)amtr(i,1)
c      enddo


      end subroutine fulhamt

******      

      end module fullmat

