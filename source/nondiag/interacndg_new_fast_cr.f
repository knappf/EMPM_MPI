!*     Phoninterac contains routines for calculation of redefined phonon
!*     interaction in p-h J-scheme

!*     last update 11.6.2013

!      modification used in restricted Pb208 calculation up to 3 hbar omega
!      only some elements are calculated 

      module phoninteracn

      contains
      
      subroutine vintn1(nf,ipcal,jcal,levn,levp)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_ndgi_int.inc'
      include 'input_ndgi_int.inc'
      include 'formats_ndgi_int.inc'


      double precision, dimension(:,:,:,:,:),allocatable ::
     *fp,fn,fpn

      integer, dimension(:), allocatable :: jphon,jphonm,ironp,iropp,
     *ironh,iroph,ndcamn,ndcamp,iphon,iphonm,phonus

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      type(level_typ),dimension(:), allocatable :: levn,levp

      character*10 namer
      character*30 namefp,namefn,namefpn,namecp,namecn,namerpp,
     *namernp,namerph,namernh,namev,nameff
     

            
      ndamp=3300
      ifmx=3300
      ifmmx=24000
      allocate (jphon(ifmx))  ! 1phonon 
      jphon=0
      allocate (iphon(ifmx))  ! 1phonon
      iphon=0
     

      allocate (jphonm(ifmmx)) ! n-1 phonon
      jphonm=0
      allocate (iphonm(ifmx)) 
      iphonm=0



      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
       iphon(i)=ipar
      enddo

      close(3)

      ifmx=i
      ndlam=ifmx

      if (nf.eq.2) nameff='1phonon/1f_states.dat'
      if (nf.eq.3) nameff='2phonon/2f_states.dat'


      open (3,file=nameff,
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphonm(i)=ijj
       iphonm(i)=ipar
      enddo

      close(3)

      ifmmx=i
   
      allocate(phonus(ifmx))
      phonus=0

      open(2,file='0hom_phon.dat',status='old',form='formatted')

!      phonus=0
!      phonmus=0

      do i=1,ifmx
      phonus(i)=0
      enddo


      do while (.not.eof(2))
      read(2,*)i,ee
!       write(*,*)i
       phonus(i)=1
      enddo

      close(2)

      open(2,file='1hom_phon.dat',status='old',form='formatted')

!      phonus=0
!      phonmus=0


      do while (.not.eof(2))
      read(2,*)i,ee
!       write(*,*)i
       phonus(i)=1
      enddo

      close(2)


      do i=1,ifmx
       phonus(i)=1
      enddo


              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
              namefp='f_p.dat'
              namefn='f_n.dat'
              namefpn='f_pn.dat'

              if (nf.eq.2) then 
                namerpp='1phonon/1f_rpp.dat'
                namernp='1phonon/1f_rnp.dat'
                namerph='1phonon/1f_rph.dat'
                namernh='1phonon/1f_rnh.dat'
                namev='Vint_phon12.dat'
              endif

              if (nf.eq.3) then 
                namerpp='2phonon/2f_rpp.dat'
                namernp='2phonon/2f_rnp.dat'
                namerph='2phonon/2f_rph.dat'
                namernh='2phonon/2f_rnh.dat'
                namev='Vint_phon23.dat'
              endif





!c     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihpmn,ippmx,ihpmn,ippmx,fp)
      call readfin(namefn,jmin,jmax,ihnmn,ipnmx,ihnmn,ipnmx,fn)
!c     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihpmn,ippmx,ihnmn,ipnmx,fpn)
!c     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
!c     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      
      write(*,*)'Calculation of redefined interaction'

!       open(33,file=namernp,status='old'
!     *,form='unformatted')

!       open(34,file=namerpp,status='old'
!     *,form='unformatted')

!       open(43,file=namernh,status='old'
!     *,form='unformatted')

!       open(44,file=namerph,status='old'
!     *,form='unformatted')


      open(5,file=namev,status='unknown'
     *,form='unformatted') !,access='append')



!   POZOR na vymenu indexov ia ig oproti vzoreckom v poznamkach

      do ia=1,ifmmx

      write(*,*)ia,ifmmx

       jia=jphonm(ia)
       ipa=iphonm(ia)

!      if (jia.eq.jcal.and.ipa.eq.ipcal) then 
!        if (jia.eq.jia) then 

!c       xfact=(dfloat(2*jia+1))**(-1.d0)

!c      write(5)ia

           call readro('1f_rnp.dat',ia,ronp,nronp)
           call readro('1f_rpp.dat',ia,ropp,nropp)
           call readro('1f_rnh.dat',ia,ronh,nronh)
           call readro('1f_rph.dat',ia,roph,nroph)

!      call readro(33,ia,ronp,nronp)
!      call readro(34,ia,ropp,nropp)
!      call readro(43,ia,ronh,nronh)
!      call readro(44,ia,roph,nroph)




      do ig=1,ifmmx

       jig=jphonm(ig)
       ipg=iphonm(ig)


       if (jig.eq.jcal.and.ipg.eq.ipcal) then


       xfact=(dfloat(2*jig+1))**(-1.d0)
      
       jisiprev=-100

       do isi=1,ifmx 
!
       ips=iphonm(isi)  ! parita sigma
      

        vint=0.d0

!       if (phonus(isi).ne.0) then

 
       if ((ips*ipa).eq.ipg) then   ! paritne pravidlo plynuce z X 
    

       jisi=jphon(isi)

       jmaxi=jia+jisi
       jmini=iabs(jia-jisi)

      if (jcal.ge.jmini.and.jcal.le.jmaxi) then  !  pravidlo pre uhlovy moment plynuce z X
       

       ifaz=(-1)**(jia+jig+jisi)



      if (jisi.ne.jisiprev) then 

       if (allocated(ironp).eq..TRUE.) deallocate (ironp)
       if (allocated(iropp).eq..TRUE.) deallocate (iropp)
       if (allocated(ironh).eq..TRUE.) deallocate (ironh)
       if (allocated(iroph).eq..TRUE.) deallocate (iroph)


      allocate(ironp(nronp))
      ironp=0 
      allocate(iropp(nropp))
      iropp=0 
      allocate(ironh(nronh))
      ironh=0 
      allocate(iroph(nroph))
      iroph=0 

      call rosub(jisi,ig,ropp,nropp,iropp,nropps)
      call rosub(jisi,ig,ronp,nronp,ironp,nronps)     
      call rosub(jisi,ig,roph,nroph,iroph,nrophs)
      call rosub(jisi,ig,ronh,nronh,ironh,nronhs)

      endif

      jisiprev=jisi
              
     
c      do ipp=ippmn,ippmx
c       jh1=lev(ih1)%j 
c      do ihp=ihpmn,ihpmx
c       jh2=lev(ih2)%j
             
        do ii=1,nropps
         i1=ropp(iropp(ii))%i1
         i2=ropp(iropp(ii))%i2
         ji1=levp(i1)%j
         ji2=levp(i2)%j
         iffz=(-1)**(jig-jia+(ji1-ji2)/2)
          do jj=1,ndcamp(isi)
           ip=camp(isi,jj)%par
           ih=camp(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camp(isi,jj)%am
           vint=vint+0.5d0*campn*fp(jisi,ip,ih,i2,i1)*ropp(iropp(ii))%ro
     **dfloat(iffz)

         enddo

         do jj=1,ndcamn(isi)
           ip=camn(isi,jj)%par
           ih=camn(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camn(isi,jj)%am
         vint=vint+campn*fpn(jisi,i2,i1,ip,ih)*ropp(iropp(ii))%ro
     **dfloat(iffz)

         enddo

        enddo

        do ii=1,nrophs
         i1=roph(iroph(ii))%i1
         i2=roph(iroph(ii))%i2
         ji1=levp(i1)%j
         ji2=levp(i2)%j
         iffz=(-1)**(jig-jia+(ji1-ji2)/2)
          do jj=1,ndcamp(isi)
           ip=camp(isi,jj)%par
           ih=camp(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camp(isi,jj)%am
         vint=vint+0.5d0*campn*fp(jisi,ip,ih,i2,i1)*roph(iroph(ii))%ro
     **dfloat(iffz)
         enddo

         do jj=1,ndcamn(isi)
           ip=camn(isi,jj)%par
           ih=camn(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camn(isi,jj)%am
         vint=vint+campn*fpn(jisi,i2,i1,ip,ih)*roph(iroph(ii))%ro
     **dfloat(iffz)     
         enddo

        enddo



        do ii=1,nronps
         i1=ronp(ironp(ii))%i1
         i2=ronp(ironp(ii))%i2
         ji1=levn(i1)%j
         ji2=levn(i2)%j
         iffz=(-1)**(jig-jia+(ji1-ji2)/2)
         
          do jj=1,ndcamn(isi)
           ip=camn(isi,jj)%par
           ih=camn(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camn(isi,jj)%am
         vint=vint+0.5d0*campn*fn(jisi,ip,ih,i2,i1)*ronp(ironp(ii))%ro
     **dfloat(iffz)

         enddo

         do jj=1,ndcamp(isi)
           ip=camp(isi,jj)%par
           ih=camp(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camp(isi,jj)%am
         vint=vint+campn*fpn(jisi,ip,ih,i2,i1)*ronp(ironp(ii))%ro
     **dfloat(iffz)    

         enddo

        enddo

        do ii=1,nronhs
         i1=ronh(ironh(ii))%i1
         i2=ronh(ironh(ii))%i2
         ji1=levn(i1)%j
         ji2=levn(i2)%j
         iffz=(-1)**(jig-jia+(ji1-ji2)/2)
         
          do jj=1,ndcamn(isi)
           ip=camn(isi,jj)%par
           ih=camn(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camn(isi,jj)%am
         vint=vint+0.5d0*campn*fn(jisi,ip,ih,i2,i1)*ronh(ironh(ii))%ro
     **dfloat(iffz)    
         enddo

         do jj=1,ndcamp(isi)
           ip=camp(isi,jj)%par
           ih=camp(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camp(isi,jj)%am
         vint=vint+campn*fpn(jisi,ip,ih,i2,i1)*ronh(ironh(ii))%ro
     **dfloat(iffz)
         enddo

        enddo

c        if (dabs(vint).gt.xrotrunc) 
c     *write(998,*)ia,ig,isi,xfact*dfloat(ifaz)*vint

        endif

        endif

        if (dabs(vint).gt.xrotrunc) 
     *write(5)ig,ia,isi,xfact*dfloat(ifaz)*vint

        if (dabs(vint).gt.xrotrunc)
     *write(998,*)ig,ia,isi,xfact*dfloat(ifaz)*vint



!c        enddo ! loop ih2
!c      enddo ! loop ih1


      enddo ! loop isi

      endif


      enddo ! loop ib
!      endif
!c      write(5)0,0,0,0,0.d0      
      enddo ! loop ig      

      
      deallocate(camp,camn,jphon,jphonm,fp,fpn,ronp,ronh,ropp,roph)
      deallocate (ironp,iropp,ironh,iroph)

c      write(5)10000000
c      write(5)0,0,0,0,0.d0 

!      close(33)
!      close(34)
!      close(43)
!      close(44)
      close(5)
      return
      end subroutine vintn1



***************************************************************************

      subroutine vintn2

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_ndgi_int.inc'
      include 'input_ndgi_int.inc'
      include 'formats_ndgi_int.inc'


      double precision, dimension(:,:,:,:,:),allocatable ::
     *fp,fn,fpn

      integer, dimension(:), allocatable :: jphon,jphonm,ironp,iropp,
     *ironh,iroph,ndcamn,ndcamp

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
!      type(level_typ),dimension(*) :: lev

      character*10 namer
      character*30 namefp,namefn,namefpn,namecp,namecn,namerpp,
     *namernp,namerph,namernh,namev,nameff
     

            
      ndamp=3200

      ifmx=3200
      ifmmx=24000
      allocate (jphon(ifmx))  ! 1phonon 
      jphon=0

      allocate (jphonm(ifmmx)) ! n-1 phonon
      jphonm=0



      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
      enddo

      close(3)

      ifmx=i
      ndlam=ifmx

     

              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
              namefp='f_p.dat'
              namefn='f_n.dat'
              namefpn='f_pn.dat'
              namev='Vint_phon2.dat'
c              ipmin=ipnmn
c              ipmax=ipnmx
c              ihmin=ihnmn
c              ihmax=ihnmx
c              ihminp=ihpmn
c              ipmaxp=ippmx
 
 


c     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihpmn,ippmx,ihpmn,ippmx,fp)
      call readfin(namefn,jmin,jmax,ihnmn,ipnmx,ihnmn,ipnmx,fn)
c     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihpmn,ippmx,ihnmn,ipnmx,fpn)
c     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
c     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      
      write(*,*)'Calculation of redefined interaction'



      open(5,file=namev,status='unknown'
     *,form='unformatted')

      do ia=1,ifmx

      write(*,*)ia,ifmx

       jig=jphon(ia)

c      write(5)ia


       do isi=1,ifmx 


        vint=0.d0


       jisi=jphon(isi)

       ifaz=(-1)**(jig+jisi)


       if (jig.eq.jisi) then 

              
             
        do ii=1,ndcamp(ia)
           ipa=camp(ia,ii)%par
           iha=camp(ia,ii)%hol
           campp=camp(ia,ii)%am

          do jj=1,ndcamp(isi)
           ips=camp(isi,jj)%par
           ihs=camp(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camp(isi,jj)%am

         vint=vint+0.25d0*campn*campp*fp(jisi,ipa,iha,ips,ihs)

         enddo

         enddo


          do ii=1,ndcamn(ia)
           ipa=camn(ia,ii)%par
           iha=camn(ia,ii)%hol
           campp=camn(ia,ii)%am

          do jj=1,ndcamn(isi)
           ips=camn(isi,jj)%par
           ihs=camn(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camn(isi,jj)%am

         vint=vint+0.25d0*campn*campp*fn(jisi,ipa,iha,ips,ihs)

         enddo

         enddo

         do ii=1,ndcamp(ia)
           ipa=camp(ia,ii)%par
           iha=camp(ia,ii)%hol
           campp=camp(ia,ii)%am

          do jj=1,ndcamn(isi)
           ips=camn(isi,jj)%par
           ihs=camn(isi,jj)%hol
c           jp=lev(ip)%j
           campn=camn(isi,jj)%am

         vint=vint+campn*campp*fpn(jisi,ipa,iha,ips,ihs)

         enddo

         enddo

          endif




        if (dabs(vint).gt.xrotrunc) 
     *write(997,*)ia,isi,vint

        if (dabs(vint).gt.xrotrunc) 
     *write(5)ia,isi,vint


      enddo ! loop isi
      enddo ! loop ib

      
      deallocate(camp,camn,jphon,fp,fpn)

c      write(5)10000000
c      write(5)0,0,0,0,0.d0 

      close(33)
      close(34)
      close(43)
      close(44)
      close(5)
      return
      end subroutine vintn2


************************************************************************


      subroutine readcam(fname,ndimi,ndimj,cam,ndcc)

      implicit double precision (a-h,o-z)

      include 'formats_ndgi_int.inc'
      include 'types_ndgi_int.inc'

      type (amp_typ), dimension(:,:), allocatable :: cam
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

      allocate(cam(ndimi,ndimj))
      allocate(ndcc(ndimi))
  
      open(2,file=fname,status='old',form='unformatted')

      ilam=0

      do while (.not.eof(2))
       ilam=ilam+1
       if (ilam.gt.ndimi) then 
               write(*,*)'Readcam: allocate bigger array in ndimi'
               stop
           endif
       read(2)ipar,ijj,ndc
         if (ndc.gt.ndimj) then 
            write(*,*)'Readcam: allocate bigger array in ndimj'
               stop
           endif
  
       read(2)(cam(ilam,i)%par,cam(ilam,i)%hol,cam(ilam,i)%am,i=1,ndc)
       ndcc(ilam)=ndc
      enddo

      close(2)

      return
      end subroutine readcam



!***********************************************************************
      subroutine rosub(j,ilam,rop,nrop,irop,nrops)

      implicit double precision (a-h,o-z)

      include 'types_ndgi_int.inc'

      type(rho_typ), dimension(:), allocatable :: rop
      integer, dimension(:), allocatable :: irop

      ii=0
      do i=1,nrop
      if (rop(i)%j.eq.j.and.rop(i)%ilap.eq.ilam) then
              ii=ii+1
              irop(ii)=i
             endif
 
      enddo

      nrops=ii

      end subroutine rosub

!***********************************************************************

      subroutine readfin(fname,jmin,jmax,imin,imax,kmin,kmax,fpp)

      implicit double precision (a-h,o-z)

      include 'formats_ndgi_int.inc'

      double precision, dimension(:,:,:,:,:), allocatable ::
     *fpp

      character(len=30)fname

      allocate(fpp(jmin:jmax,
     *imin:imax,imin:imax,kmin:kmax,kmin:kmax))
      fpp=0.d0

c      open(2,file=fname,status='old',form='formatted')
      open(2,file=fname,status='old',form='unformatted')

      do while (.not.eof(2))
c       read(2,10)itt,ipt,ijt,i,j,k,l,vint
       read(2)itt,ipt,ijt,i,j,k,l,vint


c        if (ipt.ne.ipar) goto 11
      if (ijt.gt.jmax.or.ijt.lt.jmin) goto 11
      
      if (i.gt.imax.or.i.lt.imin
     *.or.j.gt.imax.or.j.lt.imin
     *.or.k.gt.kmax.or.k.lt.kmin
     *.or.l.gt.kmax.or.l.lt.kmin) goto 11


      fpp(ijt,i,j,k,l)=vint
  
 11   enddo
      close(2)


      return
      end subroutine readfin

! 
!*****************************************************************************
      subroutine readro(fname,ig,ron,ndgg)

      implicit double precision (a-h,o-z)


      include 'formats_ndgi_int.inc'
      include 'types_ndgi_int.inc'

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=10)fname
      character(len=4)nlam
      logical je_tam

      ifile=33

      write(nlam,'(i4.4)')ig


      inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam)
 
      if (je_tam.eq..FALSE.) then 

       write(*,*)'WARNING: ',''//fname//'_'//nlam,' not present!'

        ndgg=0
        return        
      endif
     

      open(ifile,file='scratch/'//fname//'_'//nlam,
     *status='unknown',form='unformatted')

      ndro=5000000
      ndgg=0

      if (.not.allocated(ron)) allocate (ron(ndro))

       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readro'
                stop
       endif

       read(ifile)(ron(ii)%ilap,ii=1,ndgg)
       read(ifile)(ron(ii)%j,ii=1,ndgg)
       read(ifile)(ron(ii)%i1,ii=1,ndgg)
       read(ifile)(ron(ii)%i2,ii=1,ndgg)
       read(ifile)(ron(ii)%ro,ii=1,ndgg)


       close(ifile)
      return
      end subroutine readro

!***********************************************************************


     
      end module phoninteracn 
      
      
