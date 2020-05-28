!     last update 17.4 .2015
      module rdens

      contains 

!************************************************************************
      subroutine roo(jamax,jbmax,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax,iphous,iphous2,phonbs,nphon,ns1,ns2,myid,numprocs)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (ro_typ), dimension(:), allocatable :: roc
      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,iphous2,ipozl
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:,:), allocatable :: rog

      double precision, dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph

      character(len=30) fnamex,fnamec
      character(len=30) fnamer1,fnamer2,fnamer,fnameo
      character(len=10) fnamern,name_new
      integer*8 ndimroc,iroc

      integer myid,numprocs



      open(23,file='2phonon/2f_rnh.dat',status='unknown',form='unformatted') 
      open(32,file='2phonon/2f_rnp.dat',status='unknown',form='unformatted')
      open(34,file='2phonon/2f_rph.dat',status='unknown',form='unformatted')
      open(43,file='2phonon/2f_rpp.dat',status='unknown',form='unformatted')

      jjmx=jamax

      allocate (csixj(0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx))

      csixj=0.d0
  
      write(*,*)'CG coef. for jmx =',jjmx

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

      ndimroc=2100000000
      allocate(roc(ndimroc))

      roc%rho=0.d0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0
    
      iroc=1

!      call loadphon(phonbs,nphon)

!      write(*,*)' Number of 1phonon states ',nphon(1)
!      write(*,*)' Number of 2phonon states ',nphon(2)

      nff=2

!      call selphon(nff,phonbs,nphon,iphous,ns1)
!      call selphon2(nff,phonbs,nphon,iphous2,ns2)

      isi_mn=0
      isi_mx=isimax



      if (nff.eq.2) then 
       fnamex='2phonon/2f_x.dat'
       fnamec='2phonon/2f_c.dat'
      endif

      allocate(rog(nphon(2),nphon(1),0:isimax))
      rog=0.d0

      allocate(ipozl(nphon(2)))
      ipozl=0

      call readc2_part(fnamex,xx,ns2,nphon(1),nphon(1),iphous2,ipozl)

      write(*,*)'Loading 1-phonon densities'

!      call readro1_all(ro1_all,nphon(1),ihpmn,ipnmx,ihpmn,ipnmx,isi_mx)
      call  readro1_all(rpp,rnp,rph,rnh,nphon(1),ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax)

      do ipa=-1,1,2
       do ja=0,jamax

        write(*,*)'Ipar, J ',ipa,ja

        call readx(fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)

        do iaaa=ilax+1,ilax+nb

       enerf=phonbs(nff,iaaa)%enf

       if (iroc.ne.1) write(*,*)'iroc =', iroc

       iroc=1

        if (iphous(iaaa).eq.1) then 

        write(*,*)'ia',iaaa
        write(*,*)'Phonon energy ',phonbs(nff,iaaa)%enf

        ibb=ilac

           do il=1,nphon(1) ! lambda'

             rog=0.d0
             jilap=phonbs(1,il)%jj

             do i=1,nbb

            il1=cc(iaaa,i)%is 
            ila=cc(iaaa,i)%ig
            jila=phonbs(1,ila)%jj
            jil1=phonbs(1,il1)%jj

!      do isi=0,isimax 

        do ibbb=1,nphon(nff)
          jb=phonbs(nff,ibbb)%jj

          isi_mn=iabs(ja-jb)
          isi_mx=ja+jb

          do isi=isi_mn,isi_mx

            if (dabs(csixj(isi,jilap,jila,jil1,ja,jb)).gt.1.0d-9) then

            ifaz=(-1)**(isi+jilap+ja+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0
             
   
       rog(ibbb,ila,isi)=rog(ibbb,ila,isi)+fact*cc(iaaa,i)%am*xx(ibbb,il1,il)*csixj(isi,jilap,jila,jil1,ja,jb)

            endif

            if (dabs(csixj(isi,jilap,jil1,jila,ja,jb)).gt.1.0d-9) then

            ifaz=(-1)**(isi+jila+jb+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0
             
       rog(ibbb,il1,isi)=rog(ibbb,il1,isi)+fact*cc(iaaa,i)%am*xx(ibbb,il,ila)*csixj(isi,jilap,jil1,jila,ja,jb)
              
            endif

            enddo   ! over ibbb

            enddo  ! isi
      
            enddo ! over i

!            endif


          do ibbb=1,nphon(nff)

!          isi_mn=iabs(ja-jb)
!          isi_mx=ja+jb

            do ila=1,nphon(1)
!            do isi=isi_mn,isi_mx

              do isi=0,isimax

            if (dabs(rog(ibbb,ila,isi)).gt.1.d-8) then 
             roc(iroc)%rho=rog(ibbb,ila,isi)       
             roc(iroc)%ib=ibbb
             roc(iroc)%is=isi
             roc(iroc)%i1=ila
             roc(iroc)%i2=il

            iroc=iroc+1

            if (iroc.gt.ndimroc) then 
            write(*,*)' Increase dimension of array roc '
            stop
            endif
            endif

             enddo
            enddo ! ila
          
            enddo

           enddo  !over il

!           enddo   ! over isi


         endif
          

        do ityp=-1,1,2 
          do iph=-1,1,2

       if (ityp.eq.-1.and.iph.eq.1) then
!        fnamer='1phonon/1f_rpp.dat'
!        fnamern='1phonon/1f_rpp.dat'
        fnamern='1f_rpp.dat'
        name_new='2f_rpp.dat'

        ifile=43
        imin=ippmn
        imax=ippmx
       endif

       if (ityp.eq.-1.and.iph.eq.-1) then
!        fnamer='1phonon/1f_rph.dat'
!        fnamern='1phonon/1f_rph.dat'
        fnamern='1f_rph.dat'
        name_new='2f_rph.dat'

        ifile=34
        imin=ihpmn
        imax=ihpmx
       endif


       if (ityp.eq.1.and.iph.eq.1) then
!        fnamer='1phonon/1f_rnp.dat'
!        fnamern='1phonon/1f_rnp.dat'
        fnamern='1f_rnp.dat'
        name_new='2f_rnp.dat'

        imin=ipnmn
        imax=ipnmx

        ifile=32
       endif

       if (ityp.eq.1.and.iph.eq.-1) then
!        fnamer='1phonon/1f_rnh.dat'
!        fnamern='1phonon/1f_rnh.dat'
        fnamern='1f_rnh.dat'
        name_new='2f_rnh.dat'


        ifile=23
        imin=ihnmn
        imax=ihnmx

       endif
        
!        call ropar2(imin,imax,ityp,iph,jamax,jbmax,isimax,roc,iroc-1,ifile,phonbs,nphon,fnamern,iaaa,name_new)
        call ropar3(imin,imax,ityp,iph,jamax,jbmax,isimax,roc,iroc-1,ifile,phonbs,nphon,fnamern,iaaa,name_new,rnp,rpp,rnh,rph)  

        enddo
        enddo

         continue

         enddo

       enddo

      enddo

      return
      end subroutine roo

!************************************************************************
      subroutine ropar3(imin,imax,ityp,iph,jamax,jbmax,isimax,rom,irom,ifile,phonbs,nphon,fnamern,iaaa,name_new,rnp,rpp,rnh,rph)

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
      integer, dimension (:), allocatable :: nphon
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:,:,:,:), allocatable :: rho,rhon
      double precision, dimension (:,:), allocatable :: rr

      double precision, dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph

      character(len=30) fnamer1,fnamer2,fnamer,fnameo

      character(len=10) fnamern,name_new
      integer*8 irom

      character*6 nlam

      iirg=0

      if (irom.ne.0) then

      isimx=isimax

      nff=2
      ndla=nphon(nff-1)

      if (.not.allocated(rr)) allocate(rr(nphon(2),0:isimx))

      rr=0.d0

!      call readrho(fnamern,rhon,nphon,isimx,imin,imax)

       ndrh=100000000
       if (.not.allocated(rh)) allocate(rh(ndrh))

        iirg=0

        do i1=imin,imax
         do i2=imin,imax

        rr=0.d0

!$omp parallel default(shared) private(ii, &
!$omp isir,iba,ilam,ilamp)
!$omp do schedule (dynamic)

         do ii=1,irom

         isir=rom(ii)%is
         iba=rom(ii)%ib
         ilam=rom(ii)%i1
         ilamp=rom(ii)%i2  
         if (ityp.eq.-1.and.iph.eq.1) rrr=rpp(ilam,ilamp,isir,i1,i2)
         if (ityp.eq.-1.and.iph.eq.-1) rrr=rph(ilam,ilamp,isir,i1,i2)
         if (ityp.eq.1.and.iph.eq.1) rrr=rnp(ilam,ilamp,isir,i1,i2)
         if (ityp.eq.1.and.iph.eq.-1) rrr=rnh(ilam,ilamp,isir,i1,i2)       

         rr(iba,isir)=rr(iba,isir)+rrr*rom(ii)%rho

         enddo ! over ii

!$omp end do
!$omp end parallel


         do iba=1,nphon(nff)
          do isir=0,isimx

          if (dabs(rr(iba,isir)).gt.1.d-10) then

          iirg=iirg+1
          if (iirg.gt.ndrh) then
           write(*,*)'Increase ndrh in ropar2'
           stop
          endif
          rh(iirg)%ib=iba
          rh(iirg)%is=isir
          rh(iirg)%i1=i1
          rh(iirg)%i2=i2
          rh(iirg)%rho=rr(iba,isir)
          endif
          enddo
         enddo

          enddo
         enddo

 666  continue

      endif ! irom

!      write(ifile)iaaa,iirg
      if (iirg.gt.0) then
!      write(ifile)(rh(iii)%ib,iii=1,iirg)
!      write(ifile)(rh(iii)%is,iii=1,iirg)
!      write(ifile)(rh(iii)%i1,iii=1,iirg)
!      write(ifile)(rh(iii)%i2,iii=1,iirg)
!      write(ifile)(rh(iii)%rho,iii=1,iirg)
     write(nlam,'(i6.6)')iaaa

      open(73,file='scratch/'//name_new//'_'//nlam,status='unknown',form='unformatted')
      write(73)iaaa,iirg
      write(73)(rh(iii)%ib,iii=1,iirg)
      write(73)(rh(iii)%is,iii=1,iirg)
      write(73)(rh(iii)%i1,iii=1,iirg)
      write(73)(rh(iii)%i2,iii=1,iirg)
      write(73)(rh(iii)%rho,iii=1,iirg)
      close(73)

      endif


!      if (iirg.eq.0) then
!      write(ifile)iirg
!      write(ifile)iirg
!      write(ifile)iirg
!      write(ifile)iirg
!      write(ifile)dfloat(iirg)
!      endif


      return
      end subroutine ropar3

      subroutine ropar2(imin,imax,ityp,iph,jamax,jbmax,isimax,rom,irom,ifile,phonbs,nphon,fnamern,iaaa,name_new)

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
      integer, dimension (:), allocatable :: nphon
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:,:,:,:), allocatable :: rho,rhon
      double precision, dimension (:,:), allocatable :: rr



      character(len=30) fnamer1,fnamer2,fnamer,fnameo

      character(len=10) fnamern,name_new
      integer*8 irom

      character*6 nlam

 

      iirg=0 



      if (irom.ne.0) then  

      isimx=isimax

      nff=2
      ndla=nphon(nff-1)

      if (.not.allocated(rr)) allocate(rr(nphon(2),0:15))

      rr=0.d0

  
!      call readrho(fnamern,rhon,nphon,isimx,imin,imax)

      call readro11(fnamern,rhon,ndla,imin,imax,imin,imax,isimx)


       ndrh=100000000
       if (.not.allocated(rh)) allocate(rh(ndrh))

        
        iirg=0

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

          if (dabs(rr(iba,isir)).gt.1.d-10) then

          iirg=iirg+1
          if (iirg.gt.ndrh) then
           write(*,*)'Increase ndrh in ropar2'
           stop
          endif
          rh(iirg)%ib=iba
          rh(iirg)%is=isir
          rh(iirg)%i1=i1
          rh(iirg)%i2=i2
          rh(iirg)%rho=rr(iba,isir)
          endif
          enddo
         enddo

          enddo    
         enddo

 666  continue

      endif ! irom 

!      write(ifile)iaaa,iirg
      if (iirg.gt.0) then
!      write(ifile)(rh(iii)%ib,iii=1,iirg)
!      write(ifile)(rh(iii)%is,iii=1,iirg)
!      write(ifile)(rh(iii)%i1,iii=1,iirg)
!      write(ifile)(rh(iii)%i2,iii=1,iirg)
!      write(ifile)(rh(iii)%rho,iii=1,iirg)

      write(nlam,'(i6.6)')iaaa

      open(73,file='scratch/'//name_new//'_'//nlam,status='unknown',form='unformatted')      
      write(73)iaaa,iirg
      write(73)(rh(iii)%ib,iii=1,iirg)
      write(73)(rh(iii)%is,iii=1,iirg)
      write(73)(rh(iii)%i1,iii=1,iirg)
      write(73)(rh(iii)%i2,iii=1,iirg)
      write(73)(rh(iii)%rho,iii=1,iirg)
      close(73)

      endif


!      if (iirg.eq.0) then
!      write(ifile)iirg
!      write(ifile)iirg
!      write(ifile)iirg
!      write(ifile)iirg
!      write(ifile)dfloat(iirg)
!      endif


      return
      end subroutine ropar2
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
      subroutine selphon(nf,phonbs,nphon,iphous,ns1)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous


      allocate(iphous(nphon(nf)))
      iphous=0


      xthr=100000000.0d0
      ii=0
      do i=1,nphon(nf)
       xene=phonbs(nf,i)%enf
       jjf=phonbs(nf,i)%jj

       if (xene.le.xthr) then
        ii=ii+1
        iphous(i)=1
       endif
      enddo

      write(*,*)'Energy threshold for 2 phon. dens.',xthr
      write(*,*)' Number of selected phonons a)',ii

      ns1=ii

      end subroutine selphon

!**********************************************************************
      subroutine selphon2(nf,phonbs,nphon,iphous,ns2)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous


      allocate(iphous(nphon(nf)))
      iphous=0



      ii=0
      do i=1,nphon(nf)
       xene=phonbs(nf,i)%enf
       jjf=phonbs(nf,i)%jj

       if (xene.gt.-10.0d0) then
        ii=ii+1
        iphous(i)=1
       endif
      enddo

      write(*,*)' Number of selected phonons b) ',ii
      ns2=ii

      end subroutine selphon2

!**********************************************************************
    

      subroutine readx(fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

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

      subroutine readro(fname,ifile,ia,roc,iror)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (ro_typ), dimension(:), allocatable :: roc

      character(len=30)fname



      read(ifile)iaaa,iror

      if (iaaa.ne.ia) then 
              write(*,*)' Error in reading in readro'
              stop

      endif

      if (allocated(roc)) deallocate(roc)

      allocate(roc(iror))


      read(ifile)(roc(iii)%ib,roc(iii)%is,roc(iii)%i1,roc(iii)%i2,roc(iii)%rho,iii=1,iror)

      return
      end subroutine readro
!************************************************************
      subroutine readrho(fname,rho,nphon,isimx,imin,imax)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      double precision, dimension (:,:,:,:,:), allocatable :: rho

      type(rho_typ), dimension(:), allocatable :: ron
      integer, dimension (:), allocatable :: nphon


      character(len=30)fname

      ndimx=nphon(1)

      allocate (rho(ndimx,ndimx,0:isimx,imin:imax,imin:imax))

      rho=0.d0

      ndro=100000
      allocate(ron(ndro))

      ron%ilap=0
      ron%j=0
      ron%i1=0
      ron%i2=0
      ron%ro=0.0d0


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

!***********************************************************************
      subroutine readro11(fname,ron,ndla,n1mn,n1mx,n2mn,n2mx,nsi)

      implicit double precision (a-h,o-z)

!      include 'types_eqm.inc'
      include 'types_phon_dens.inc'

      double precision, dimension(:,:,:,:,:), allocatable :: ron
      type(rho_typ), dimension(:), allocatable :: ronn
      logical je_tam_subor


      character(len=10)fname
      character(len=4)nlam

      if (.not.allocated(ron)) allocate(ron(ndla,ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))
      ron=0.d0

      ifile=33

      do ig=1,ndla

      write(nlam,'(i4.4)')ig

      inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam_subor)

      if (je_tam_subor.eq..FALSE.) then
        ndgg=0
        write(*,*)'File scratch/'//fname//'_'//nlam,'not found! '
        return
      endif


      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')


      ndro=5000000
      ndgg=0

      if (.not.allocated(ronn)) allocate (ronn(ndro))


       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readro',ndgg,ndro
                stop
       endif

       read(ifile)(ronn(ii)%ilap,ii=1,ndgg)
       read(ifile)(ronn(ii)%j,ii=1,ndgg)
       read(ifile)(ronn(ii)%i1,ii=1,ndgg)
       read(ifile)(ronn(ii)%i2,ii=1,ndgg)
       read(ifile)(ronn(ii)%ro,ii=1,ndgg)

      close(ifile)

!      if (.not.allocated(ron)) allocate(ron(ndla,ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))

!      ron=0.d0

       do ii=1,ndgg

        ibt=ronn(ii)%ilap
        isit=ronn(ii)%j
        i1t=ronn(ii)%i1
        i2t=ronn(ii)%i2
        rot=ronn(ii)%ro

        if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
         write(*,*)' Small dimensions of ron in read11',ibt,isit,i1t,i2t
         stop
       endif

        ron(ig,ibt,isit,i1t,i2t)=rot

       enddo

       enddo

      return
      end subroutine readro11


!******************************************************************************
!******************************************************************************
      subroutine readro1_all(rpp,rnp,rph,rnh,ndla,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax)

      implicit double precision (a-h,o-z)

!      include 'types_eqm.inc'
      include 'types_phon_dens.inc'

      double precision, dimension(:,:,:,:,:), allocatable :: rpp,rnp,rph,rnh
      type(rho_typ), dimension(:), allocatable :: ronn
      logical je_tam_subor


      character(len=10)fname
      character(len=4)nlam

      if (.not.allocated(rpp)) allocate(rpp(ndla,ndla,0:isimax,ippmn:ippmx,ippmn:ippmx))
      rpp=0.d0
      if (.not.allocated(rnp)) allocate(rnp(ndla,ndla,0:isimax,ipnmn:ipnmx,ipnmn:ipnmx))
      rnp=0.d0
      if (.not.allocated(rph)) allocate(rph(ndla,ndla,0:isimax,ihpmn:ihpmx,ihpmn:ihpmx))
      rph=0.d0
      if (.not.allocated(rnh)) allocate(rnh(ndla,ndla,0:isimax,ihnmn:ihnmx,ihnmn:ihnmx))
      rnh=0.d0
      
      ndro=5000000
      if (.not.allocated(ronn)) allocate (ronn(ndro))
 
      ifile=33


      do ityp=1,4
       if (ityp.eq.1) fname='1f_rnp.dat'
       if (ityp.eq.2) fname='1f_rpp.dat'
       if (ityp.eq.3) fname='1f_rnh.dat'
       if (ityp.eq.4) fname='1f_rph.dat'

      ifile=33

      do ig=1,ndla

      write(nlam,'(i4.4)')ig

      inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam_subor)

      if (je_tam_subor.eq..FALSE.) then
        ndgg=0
        write(*,*)'File scratch/'//fname//'_'//nlam,'not found! '
        return
      endif


      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')

      ndgg=0

       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readro',ndgg,ndro
                stop
       endif

       read(ifile)(ronn(ii)%ilap,ii=1,ndgg)
       read(ifile)(ronn(ii)%j,ii=1,ndgg)
       read(ifile)(ronn(ii)%i1,ii=1,ndgg)
       read(ifile)(ronn(ii)%i2,ii=1,ndgg)
       read(ifile)(ronn(ii)%ro,ii=1,ndgg)

       close(ifile)
!      if (.not.allocated(ron)) allocate(ron(ndla,ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))

!      ron=0.d0

       do ii=1,ndgg

        ibt=ronn(ii)%ilap
        isit=ronn(ii)%j
        i1t=ronn(ii)%i1
        i2t=ronn(ii)%i2
        rot=ronn(ii)%ro

!        if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
!         write(*,*)' Small dimensions of ron in read11',ibt,isit,i1t,i2t
!         stop
!       endif
       if (isit.le.isimax) then
        if (ityp.eq.1) rnp(ig,ibt,isit,i1t,i2t)=rot
        if (ityp.eq.2) rpp(ig,ibt,isit,i1t,i2t)=rot
        if (ityp.eq.3) rnh(ig,ibt,isit,i1t,i2t)=rot
        if (ityp.eq.4) rph(ig,ibt,isit,i1t,i2t)=rot

       endif
       enddo
       enddo
       enddo

      return
      end subroutine readro1_all


!******************************************************************************

      subroutine readc2_part(fname,cc,n1,n2,n3,iphous2,ipozl)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      double precision, dimension (:,:,:), allocatable :: cc
      integer, dimension (:), allocatable :: ndcc,iphous2,ipozl
      real :: xxr
      character(len=30)fname

      ndimi=2000
      ndimj=600000

      allocate (cc(n1,n2,n3))
      cc=0.0
      illl=1


      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      do while (.not.eof(2))

       read(2)ipar,ijj,no,idphon

!      if (ipar.eq.-1.and.ijj.eq.1) then
!     read(2)ipar,ijj,no,idphon
      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(1,idphon))
      allocate(ndcc(1))

      do ilam=1,no

         if (idphon.gt.ndimj) then
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif


       if (iphous2(ilam+ilamp).eq.1) then

       read(2)(xcc(1,i)%ig,xcc(1,i)%is,xcc(1,i)%am,i=1,idphon)
       ndcc(1)=idphon


       do i=1,idphon
        ii=xcc(1,i)%is
        jj=xcc(1,i)%ig
        xxr=real(xcc(1,i)%am)
        cc(illl,ii,jj)=xxr
        ipozl(ilam+ilamp)=illl
       enddo

       illl=illl+1
       else

       read(2)

       endif
      enddo

 !     endif

      ilamps=ilamp
      ilamp=ilamp+no

      enddo
  11  continue
      close(2)

      deallocate(xcc,ndcc)
      return
      end subroutine readc2_part

      end module rdens
      
