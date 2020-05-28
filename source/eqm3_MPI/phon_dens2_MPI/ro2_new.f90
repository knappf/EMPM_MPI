!     last update 17.1.2011
      module rdens

      contains 
!************************************************************
! calculation of R 

      subroutine calcr(ia,roc,iror,csixj,phonbs,nphon,iphous,xx,isimax,cc,nb,nbb)

!      subroutine roo(jamax,jbmax,isimax)


      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'
!      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (ro_typ), dimension(:), allocatable :: roc
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      integer, dimension (:), allocatable :: ndx,ndc

      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:), allocatable :: rog

      character(len=30) fnamex,fnamec



      ndimroc=1000000000
      allocate(roc(ndimroc))

      roc%rho=0.d0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0


      nff=2

!      call selphon(nff,phonbs,nphon,iphous)

!c      stop

      ndamp=10000
      ndrho=1000000     

      allocate(rog(nphon(1)))
      rog=0.d0

 
!      call readc2(fnamex,xx,nphon(nff),nphon(1),nphon(1))

           roc%rho=0.d0
           roc%ib=0
           roc%is=0
           roc%i1=0
           roc%i2=0


         do ibbb=1,nphon(nff)

         write(*,*)'ibbb ',ibbb

          jb=phonbs(nff,ibbb)%jj

          do isi=0,isimax

           do il=1,nphon(1) ! lambda'

             rog=0.d0
             jilap=phonbs(1,il)%jj

             do i=1,nbb


            il1=cc(iaaa,i)%is 
            ila=cc(iaaa,i)%ig
            jila=phonbs(1,ila)%jj
            jil1=phonbs(1,il1)%jj 

            ifaz=(-1)**(isi+jilap+ja+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0
             
   
       rog(ila)=rog(ila)+fact*cc(iaaa,i)%am*xx(ibbb,il1,il)*csixj(isi,jilap,jila,jil1,ja,jb)

            ifaz=(-1)**(isi+jila+jb+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0
             
       rog(il1)=rog(il1)+fact*cc(iaaa,i)%am*xx(ibbb,il,ila)*csixj(isi,jilap,jil1,jila,ja,jb)
      
            enddo ! over i

            do ila=1,nphon(1)

            if (dabs(rog(ila)).gt.1.d-8) then 
             roc(iroc)%rho=rog(ila)       
             roc(iroc)%ib=ibbb
             roc(iroc)%is=isi
             roc(iroc)%i1=ila
             roc(iroc)%i2=il

!c             write(998,17)iaaa,ibbb,isi,ila,il,roc(iroc)%rho

            iroc=iroc+1

            if (iroc.gt.ndimroc) then 
            write(*,*)' Increase dimension of array roc '
            stop
            endif
            endif

            enddo

           enddo  !over il

           enddo   ! over isi

         enddo   !over iaaa

          
!         write(*,*)'iroc', iroc-1
!c         write(998,*)(roc(iii)%rho,iii=1,iroc)
!         write(22)iaaa,iroc-1
!         write(22)(roc(iii)%ib,roc(iii)%is,roc(iii)%i1,roc(iii)%i2,roc(iii)%rho,iii=1,iroc-1)

!          if ((iroc-1).gt.0) then 
!           roc%rho=0.d0
!           roc%ib=0
!           roc%is=0
!           roc%i1=0
!           roc%i2=0
!          endif

         continue


      return
      end subroutine calcr

!************************************************************************

      subroutine ropar(imin,imax,ityp,iph,jamax,jbmax,isimax)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (ro_typ), dimension(:), allocatable :: roc,rom,rh
      type(rho_typ), dimension(:), allocatable :: rop

      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:,:,:,:), allocatable :: rho,rhon
      double precision, dimension (:,:), allocatable :: rr



      character(len=30) fnamer1,fnamer2,fnamer,fnamern,fnameo,fnamex,fnamec

!      open(22,file='ro1.dat',status='unknown',form='unformatted')

      jjmx=14
      isimx=isimax


      call loadphon(phonbs,nphon)

      ndamp=10000
      ndrho=10000000     
      nff=2

      write(*,*)' Number of 1phonon states ',nphon(1)
      write(*,*)' Number of 2phonon states ',nphon(2)
      write(*,*)' Imin, Imax ',imin,imax

      if (nff.eq.2) then 
!c       fnamer1='2phonon/ro1.dat'
!       fnamer2='2phonon/roo.dat'

       if (ityp.eq.-1.and.iph.eq.1) then
        write(*,*)' Calculating proton particle densities'  
        fnamer='1phonon/1f_rpp.dat'
        fnamern='1phonon/1f_rpp.dat'
        fnameo='2phonon/2f_rpp.dat'
       endif

       if (ityp.eq.-1.and.iph.eq.-1) then
        write(*,*)' Calculating proton hole densities'  
        fnamer='1phonon/1f_rph.dat'
        fnamern='1phonon/1f_rph.dat'
        fnameo='2phonon/2f_rph.dat'
       endif


       if (ityp.eq.1.and.iph.eq.1) then  
        write(*,*)' Calculating neutron particle densities' 
        fnamer='1phonon/1f_rnp.dat'
        fnamern='1phonon/1f_rnp.dat'
        fnameo='2phonon/2f_rnp.dat'
       endif

       if (ityp.eq.1.and.iph.eq.-1) then
        write(*,*)' Calculating neutron hole densities'  
        fnamer='1phonon/1f_rnh.dat'
        fnamern='1phonon/1f_rnh.dat'
        fnameo='2phonon/2f_rnh.dat'
       endif



      endif

      call selphon(2,phonbs,nphon,iphous)


      allocate(rr(nphon(2),0:15))
      rr=0.d0


!      call readrho(fnamer,rho)
!      call readrho(fnamern,rhon)
      call readrho(fnamern,rhon,nphon(1),isimax,imin,imax)  

      call precal6j(jjmx,csixj)

      open(33,file=fnameo,status='unknown',form='unformatted')

      ndrh=100000000
      allocate(rh(ndrh))

       fnamex='2phonon/2f_x.dat'
       fnamec='2phonon/2f_c.dat'

       call readc2(fnamex,xx,nphon(2),nphon(1),nphon(1))



      do ipa=-1,1,2
       do ja=0,jamax

        write(*,*)'Ipar, J ',ipa,ja

        call readx(fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)

        do iaaa=ilax+1,ilax+nb

        write(*,*)'iaaa',iaaa


       enerf=phonbs(nff,iaaa)%enf

       write(*,*)'Phonon energy ',phonbs(nff,iaaa)%enf

       iroc=1

        if (iphous(iaaa).eq.1) then

        ibb=ilac





 
!      do iaaa=1,nphon(nff)
        
        iirg=0

        call calcr(iaaa,rom,irom,csixj,phonbs,nphon,iphous,xx,isimax,cc,nb,nbb)  ! calculates all R's for a give iaaaa
        write(*,*)'***',iaaa,iror,irom

        do i1=imin,imax
         do i2=imin,imax

        rr=0.d0


         do ii=1,irom

         isir=rom(ii)%is
         iba=rom(ii)%ib
         ilam=rom(ii)%i1
         ilamp=rom(ii)%i2
 
       rr(iba,isir)=rr(iba,isir)+rhon(ilam,ilamp,isir,i1,i2)*rom(ii)%rho

            
         enddo ! over ii


         do iba=1,nphon(nff)
          do isir=0,isimx
          iirg=iirg+1
          rh(iirg)%ib=iba
          rh(iirg)%is=isir
          rh(iirg)%i1=i1
          rh(iirg)%i2=i2
          rh(iirg)%rho=rr(iba,isir)
          enddo
         enddo

          enddo    
         enddo

      write(*,*)'iirg = ',iirg   
 666  continue
      write(33)iaaa,iirg
      if (iirg.gt.0) then
      write(33)(rh(iii)%ib,iii=1,iirg)
      write(33)(rh(iii)%is,iii=1,iirg)
      write(33)(rh(iii)%i1,iii=1,iirg)
      write(33)(rh(iii)%i2,iii=1,iirg)
      write(33)(rh(iii)%rho,iii=1,iirg)
      endif

      if (iirg.eq.0) then
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)iirg
      write(33)dfloat(iirg)
      endif


       endif

      enddo           !    over iaaa


      enddo
      enddo  

      close(33)
      close(88)
      close(77)

      deallocate(rh,rr,rom)

      return
      end subroutine ropar
!************************************************************************

      subroutine loadphon(phonbs,nphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon


      allocate (nphon(0:4))
      nphon=0
      allocate (phonbs(0:4,1:1000000))

      do nfon=1,2

      if (nfon.eq.1) namef='1phonon/1f_states.dat'
      if (nfon.eq.2) namef='2phonon/2f_states.dat'
      
      open(1,file=namef,status='old',form='unformatted')

      do while (.not.eof(1))
       read(1)i,ipartt,ijj,en
         phonbs(nfon,i)%par=ipartt
         phonbs(nfon,i)%jj=ijj
         phonbs(nfon,i)%enf=en
      enddo
      nphon(nfon)=i

      close(1)

      enddo

      end subroutine loadphon

!*********************************************************************
      subroutine selphon(nf,phonbs,nphon,iphous)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous


      allocate(iphous(nphon(nf)))
      iphous=0


      write(*,*)'energy threshold for 2 phonon states'
      read(*,*)xthre
      write(*,*)xthre

!      xthre=30.0d0


      ii=0
      do i=1,nphon(nf)
       xene=phonbs(nf,i)%enf
       if (xene.le.xthre) then
        ii=ii+1
        iphous(i)=1
       endif
      enddo

      write(*,*)' Number of selected phonons ',ii

      end subroutine selphon

!********************************************************************** 
    
      subroutine readx(fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

!c      allocate(xcc(ndimi,ndimj))
!c      allocate(ndcc(ndimi))

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(ilamp+1:ilamp+no+1,ndimj))
      allocate(ndcc(ilamp+1:ilamp+1+no))

      xcc%is=0
      xcc%ig=0
      xcc%am=0.d0
      ndcc=0


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon
      enddo

      ilamps=ilamp
      ilamp=ilamp+no


      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


      deallocate(xcc,ndcc)
      no=0
      idphon=0

      enddo

  11  continue     

      close(2)

      return
      end subroutine readx
!************************************************************
      subroutine readx2(fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

!c      allocate(xcc(ndimi,ndimj))
!c      allocate(ndcc(ndimi))

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000



      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(ilamp+1:ilamp+no+1,ndimj))
      allocate(ndcc(ilamp+1:ilamp+1+no))

      xcc%is=0
      xcc%ig=0
      xcc%am=0.d0
      ndcc=0


      do ilam=1,no


         if (idphon.gt.ndimj) then
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif

       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon
      enddo

      ilamps=ilamp
      ilamp=ilamp+no


      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


      deallocate(xcc,ndcc)
      no=0
      idphon=0

      enddo

  11  continue

      close(2)

      return
      end subroutine readx2

      subroutine cutx(fname)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=16)fname
      character*7 nlam

!c      allocate(xcc(ndimi,ndimj))
!c      allocate(ndcc(ndimi))

      open(2,file=fname,status='old',form='unformatted')
      

      ilamp=0

      ndimj=100000


      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(ilamp+1:ilamp+no+1,ndimj))
      allocate(ndcc(ilamp+1:ilamp+1+no))

      xcc%is=0
      xcc%ig=0
      xcc%am=0.d0
      ndcc=0

      do ilam=1,no
         if (idphon.gt.ndimj) then
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif

       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon

       write(nlam,'(i7.7)')ilam


       open(3,file=fname//'_'//nlam,status='unknown',form='unformatted')

       write(3)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%am,i=1,idphon)

      enddo

    

      ilamps=ilamp
      ilamp=ilamp+no

      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11

      deallocate(xcc,ndcc)
      no=0
      idphon=0

      enddo

  11  continue

      close(2)

      return
      end subroutine cutx



      subroutine readro(fname,ifile,ia,roc,iror)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (ro_typ), dimension(:), allocatable :: roc

      character(len=30)fname


!c      open(22,file=fname,status='old',form='unformatted')

     

!c      do while (.not.eof(22))


      read(ifile)iaaa,iror

      if (iaaa.ne.ia) then 
              write(*,*)' Error in reading in readro'
              stop

      endif

      if (allocated(roc)) deallocate(roc)

      allocate(roc(iror))


      read(ifile)(roc(iii)%ib,roc(iii)%is,roc(iii)%i1,roc(iii)%i2,roc(iii)%rho,iii=1,iror)


!c      if (iaaa.eq.ia) goto 11


!c      enddo

!c   11 continue   

!c      close(22)

      return
      end subroutine readro
!****************************************

      subroutine readc2(fname,cc,n1,n2,n3)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      double precision, dimension (:,:,:), allocatable :: cc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

      ndimi=2000
      ndimj=100

      allocate (cc(n1,n2,n3))
      cc=0.d0


      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(1,ndimj))
      allocate(ndcc(1))


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(1,i)%ig,xcc(1,i)%is,xcc(1,i)%am,i=1,idphon)
       ndcc(1)=idphon

       do i=1,idphon
        ii=xcc(1,i)%is
        jj=xcc(1,i)%ig
        cc(ilam+ilamp,ii,jj)=xcc(1,i)%am
       enddo
      enddo


   

      ilamps=ilamp
      ilamp=ilamp+no



      enddo

  11  continue     



      close(2)

      deallocate(xcc,ndcc)

      return
      end subroutine readc2
!************************************************************
      subroutine readr(fname,xcc,ndcc,iaaa)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname


      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(ilamp:ilamp+no,ndimj))
      allocate(ndcc(ilamp:ilamp+no))


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon
      enddo

      ilamps=ilamp
      ilamp=ilamp+no


      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


      deallocate(xcc,ndcc)
      no=0
      idphon=0

      enddo

  11  continue     

      close(2)

      return
      end subroutine readr
!************************************************************
      subroutine readrho(fname,rho,nfon1,isimax,imin,imax)

      implicit double precision (a-h,o-z)

!c      include 'formats_eqm.inc'
      include 'types_phon_dens.inc'

      double precision, dimension (:,:,:,:,:), allocatable :: rho

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=30)fname

!      allocate (rho(104,104,0:10,10,10))
      allocate (rho(nfon1,nfon1,0:isimax,imin:imax,imin:imax))

      rho=0.d0

      ndro=50000
      allocate(ron(ndro))


      open(2,file=fname,status='old',form='unformatted')

      do while (.not.eof(2))

       read(2)igg,ndgg

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readrho'
                stop
       endif

       read(2)(ron(ii)%ilap,ii=1,ndgg)
       read(2)(ron(ii)%j,ii=1,ndgg)
       read(2)(ron(ii)%i1,ii=1,ndgg)
       read(2)(ron(ii)%i2,ii=1,ndgg)
       read(2)(ron(ii)%ro,ii=1,ndgg)

       do i=1,ndgg
        ilapp=ron(i)%ilap
        jp=ron(i)%j
        i1p=ron(i)%i1
        i2p=ron(i)%i2
        rho(igg,ilapp,jp,i1p,i2p)=ron(i)%ro
       enddo

       enddo


      close(2)

      deallocate(ron)

      return
      end subroutine readrho

!**********************************************************************

      subroutine precal6j(jjmx,csixj)
      
      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj


      allocate (csixj(0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx))

      csixj=0.d0

      write(*,*)'calculation of 6-j up to j=',jjmx

      do i1=0,jjmx
!       write(*,*)'i1 =',i1
       do i2=0,jjmx
        do i3=0,jjmx
         do i4=0,jjmx
          do i5=0,jjmx
           do i6=0,jjmx
             csixj(i1,i2,i3,i4,i5,i6)=sixj(2*i1,2*i2,2*i3,2*i4,2*i5,2*i6)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo


      return 
      end subroutine precal6j 

      end module rdens
      
