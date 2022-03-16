!*     Phoninterac contains routines for calculation of redefined phonon
!*     interaction in p-h J-scheme

!*     last update 6.2.2020


module phoninteracn

contains

subroutine dim_phon(nf,dim_phon_sp)
implicit none

integer :: nf,dim_phon_sp,i
character*30 name_file 

if (nf.eq.1) name_file='1phonon/1f_states.dat'
if (nf.eq.2) name_file='2phonon/2f_states.dat'

open (3,file=name_file,status='old',form='unformatted')

do while (.not.eof(3))
 read(3)i
enddo

dim_phon_sp=i
write(*,*)' dimension of ',nf,'-phon subspace: ',dim_phon_sp

close(3)
return 
end subroutine dim_phon

  
 !     
subroutine vintn12(nf,ipcal,jcal,levn,levp,myid,numprocs)

use anglib   ! angular momentum staff

implicit double precision(a-h,o-z)

include 'types_ndgi_int.inc'
include 'input_ndgi_int.inc'
include 'formats_ndgi_int.inc'


double precision, dimension(:,:,:,:,:),allocatable ::fp,fn,fpn
integer, dimension(:), allocatable :: jphon,jphonm,ironp,iropp,ironh,iroph,ndcamn,ndcamp,iphon,iphonm,phonus
type(amp_typ), dimension(:,:), allocatable :: camp,camn
type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
type(level_typ),dimension(:), allocatable :: levn,levp
character*10 namer
character*30 namefp,namefn,namefpn,namecp,namecn,namerpp,namernp,namerph,namernh,namev,nameff
character(len=4)nlam


      ndamp=1600
      ifmx=1600
      ifmmx=24000

      call dim_phon(1,ifmx)
      allocate (jphon(ifmx))  ! 1phonon 
      jphon=0
      allocate (iphon(ifmx))  ! 1phonon
      iphon=0
     

      call dim_phon(nf-1,ifmmx)
      allocate (jphonm(ifmmx)) ! n-1 phonon
      jphonm=0
      allocate (iphonm(ifmmx)) 
      iphonm=0



      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

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


      open (3,file=nameff,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphonm(i)=ijj
       iphonm(i)=ipar
      enddo

      close(3)

      ifmmx=i
   
      allocate(phonus(ifmx))
      phonus=1

              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
              namefp='fmat_p.dat'
              namefn='fmat_n.dat'
              namefpn='fmat_pn.dat'

!              if (nf.eq.2) then 
                namerpp='1f_rpp.dat'
                namernp='1f_rnp.dat'
                namerph='1f_rph.dat'
                namernh='1f_rnh.dat'
                namev='Vint_phon12'
!              endif

       write(*,*)'jmin,jmax',jmin,jmax

!c     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihpmn,ippmx,ihpmn,ippmx,fp)
      call readfin(namefn,jmin,jmax,ihnmn,ipnmx,ihnmn,ipnmx,fn)
!c     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihpmn,ippmx,ihnmn,ipnmx,fpn)
!c     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
!c     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      
      if (myid == 0) write(*,*)'Calculation of redefined interaction'


      ndro=500000
       allocate(ronp(ndro))
       allocate(ropp(ndro))
       allocate(ronh(ndro))
       allocate(roph(ndro))


      write(nlam,'(i4.4)')myid
      ifunit=55+myid
      open(ifunit,file='scratch/Vint12_'//nlam,status='unknown',form='unformatted')      

!     boundaries for J a Pi 
      ia_max=1
      ia_min=ifmmx

      do ia=1,ifmmx 
         jia=jphonm(ia)
         ipa=iphonm(ia)
         if (jia.eq.jcal.and.ipa.eq.ipcal) then 
          if (ia > ia_max) ia_max=ia
          if (ia < ia_min) ia_min=ia
         endif 
      enddo
      
      if (mod(ia_max-ia_min+1,numprocs).eq.0) then
            n_seg=(ia_max-ia_min+1)/numprocs
      else
            n_seg=(ia_max-ia_min+1)/numprocs+1
      endif

      if (myid ==0 ) then 
            write(*,*) 'Interval ia',ia_min,ia_max
            write(*,*) ' Size of segment', n_seg
      endif
    


      do ia=myid*n_seg+ia_min,min((myid+1)*n_seg+ia_min-1,ia_max)    
      
      write(*,*)' Process  #',myid,' ia = ',ia,myid*n_seg+ia_min,min((myid+1)*n_seg+ia_min-1,ia_max)     

       jia=jphonm(ia)
       ipa=iphonm(ia)

       xfact=(dfloat(2*jia+1))**(-1.d0)



       if (jia.eq.jcal.and.ipa.eq.ipcal) then

           call readro(namernp,ia,ronp,nronp)
           call readro(namerpp,ia,ropp,nropp)
           call readro(namernh,ia,ronh,nronh)
           call readro(namerph,ia,roph,nroph)

       if (allocated(ironp).eq..TRUE.) deallocate (ironp)
       if (allocated(iropp).eq..TRUE.) deallocate (iropp)
       if (allocated(ironh).eq..TRUE.) deallocate (ironh)
       if (allocated(iroph).eq..TRUE.) deallocate (iroph)    

       allocate(ironp(nronp))
       allocate(iropp(nropp))
       allocate(ironh(nronh))
       allocate(iroph(nroph))

      do ig=1,ifmmx

       jig=jphonm(ig)
       ipg=iphonm(ig)

      
       jisiprev=-100

       do isi=1,ifmx 
!
       ips=iphon(isi)  ! parita sigma
      

        vint=0.d0

!       if (phonus(isi).ne.0) then

 
       if ((ips*ipg).eq.ipa) then   ! paritne pravidlo plynuce z X 
    

       jisi=jphon(isi)

       jmaxi=jig+jisi
       jmini=iabs(jig-jisi)

      if (jcal.ge.jmini.and.jcal.le.jmaxi) then  !  pravidlo pre uhlovy moment plynuce z X
       

       ifaz=(-1)**(jia+jig+jisi)



      ironp=0 
      iropp=0 
      ironh=0 
      iroph=0 

      call rosub(jisi,ig,ropp,nropp,iropp,nropps)
      call rosub(jisi,ig,ronp,nronp,ironp,nronps)     
      call rosub(jisi,ig,roph,nroph,iroph,nrophs)
      call rosub(jisi,ig,ronh,nronh,ironh,nronhs)

!      endif

      jisiprev=jisi
              
             
        do ii=1,nropps
         i1=ropp(iropp(ii))%i1
         i2=ropp(iropp(ii))%i2
         ji1=levp(i1)%j
         ji2=levp(i2)%j
         iffz=(-1)**(jig-jia+(ji1-ji2)/2)
          do jj=1,ndcamp(isi)
           ip=camp(isi,jj)%par
           ih=camp(isi,jj)%hol
!c           jp=lev(ip)%j
           campn=camp(isi,jj)%am
           vint=vint+0.5d0*campn*fp(jisi,ip,ih,i1,i2)*ropp(iropp(ii))%ro

         enddo

         do jj=1,ndcamn(isi)
           ip=camn(isi,jj)%par
           ih=camn(isi,jj)%hol
!c           jp=lev(ip)%j
           campn=camn(isi,jj)%am
           vint=vint+campn*fpn(jisi,i1,i2,ip,ih)*ropp(iropp(ii))%ro

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
!           jp=lev(ip)%j
           campn=camp(isi,jj)%am
           vint=vint+0.5d0*campn*fp(jisi,ip,ih,i1,i2)*roph(iroph(ii))%ro
         enddo

         do jj=1,ndcamn(isi)
           ip=camn(isi,jj)%par
           ih=camn(isi,jj)%hol
!           jp=lev(ip)%j
           campn=camn(isi,jj)%am
           vint=vint+campn*fpn(jisi,i1,i2,ip,ih)*roph(iroph(ii))%ro
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
!           jp=lev(ip)%j
           campn=camn(isi,jj)%am
           vint=vint+0.5d0*campn*fn(jisi,ip,ih,i1,i2)*ronp(ironp(ii))%ro

         enddo

         do jj=1,ndcamp(isi)
           ip=camp(isi,jj)%par
           ih=camp(isi,jj)%hol
!           jp=lev(ip)%j
           campn=camp(isi,jj)%am
           vint=vint+campn*fpn(jisi,ip,ih,i1,i2)*ronp(ironp(ii))%ro

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
!           jp=lev(ip)%j
           campn=camn(isi,jj)%am
           vint=vint+0.5d0*campn*fn(jisi,ip,ih,i1,i2)*ronh(ironh(ii))%ro
         enddo

         do jj=1,ndcamp(isi)
           ip=camp(isi,jj)%par
           ih=camp(isi,jj)%hol
!           jp=lev(ip)%j
           campn=camp(isi,jj)%am
           vint=vint+campn*fpn(jisi,ip,ih,i1,i2)*ronh(ironh(ii))%ro
         enddo

        enddo

!        if (dabs(vint).gt.xrotrunc) 
!     *write(998,*)ia,ig,isi,xfact*dfloat(ifaz)*vint

        endif

        endif

        if (dabs(vint).gt.xrotrunc) write(ifunit)ia,ig,isi,real(xfact*dfloat(ifaz)*vint)


      enddo ! loop isi


      enddo ! loop ib


      endif
      enddo ! loop ig      

      
      deallocate(camp,camn,jphon,jphonm,fp,fpn)
 !     if (allocated(ronp)) deallocate(ronp)
 !     if (allocated(ropp)) deallocate(ropp)
 !     if (allocated(roph)) deallocate(roph)
 !     if (allocated(ronh)) deallocate(ronh)
 !     if (allocated(ironp)) deallocate(ironp)


!      deallocate (ironp,iropp,ironh,iroph)


      close(ifunit)
      return
      end subroutine vintn12




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
      subroutine rosub2(j,ilam,rop,nrop,irop,nrops)

      implicit double precision (a-h,o-z)

      include 'types_ndgi_int.inc'

      type(rho2_typ), dimension(:), allocatable :: rop
      integer, dimension(:), allocatable :: irop

      ii=0
      do i=1,nrop
      if (rop(i)%j.eq.j.and.rop(i)%ilap.eq.ilam) then
              ii=ii+1
              irop(ii)=i
             endif
 
      enddo

      nrops=ii

      end subroutine rosub2

!***********************************************************************

      subroutine readfin(fname,jmin,jmax,imin,imax,kmin,kmax,fpp)

      implicit double precision (a-h,o-z)

      include 'formats_ndgi_int.inc'

      double precision, dimension(:,:,:,:,:), allocatable ::fpp

      character(len=30)fname
      integer(kind=1) :: j_f
      integer(kind=2) :: i_i,i_j,i_k,i_l


      allocate(fpp(jmin:jmax,imin:imax,imin:imax,kmin:kmax,kmin:kmax))
      fpp=0.d0

!      open(2,file=fname,status='old',form='formatted')
      open(2,file=fname,status='old',form='unformatted')

      do while (.not.eof(2))
!       read(2,10)itt,ipt,ijt,i,j,k,l,vint
!       read(2)itt,ipt,ijt,i,j,k,l,vint
      read(2)j_f,i_i,i_j,i_k,i_l,vint

      ijt=j_f
      i=i_i
      j=i_j
      k=i_k
      l=i_l

!        if (ipt.ne.ipar) goto 11
      if (ijt.gt.jmax.or.ijt.lt.jmin) goto 11
      
      if (i.gt.imax.or.i.lt.imin.or.j.gt.imax.or.j.lt.imin.or.k.gt.kmax.or.k.lt.kmin.or.l.gt.kmax.or.l.lt.kmin) goto 11

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
     

      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')

      ndro=size(ron)
      ndgg=0

!      if (.not.allocated(ron)) allocate (ron(ndro))

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

!*****************************************************************************
      subroutine readro2(fname,ig,ron,ndgg)

            implicit double precision (a-h,o-z)
      
            include 'formats_ndgi_int.inc'
            include 'types_ndgi_int.inc'
      
            Real(Kind=4), dimension(:), allocatable :: ron
      
            character(len=10)fname
            character(len=6)nlam
            logical je_tam
      
            ifile=33
      
            write(nlam,'(i6.6)')ig
      
            inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam)
      
            if (je_tam.eq..FALSE.) then
              write(*,*)'WARNING: ',''//fname//'_'//nlam,' not present!'
              ndgg=0
              return
            endif
      
      
            open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')
      
            ndro=10000000
            ndgg=0
      
            if (.not.allocated(ron)) allocate (ron(ndro))
      
             read(ifile)ndgg

      
             if (ndgg.gt.ndro) then
                      write(*,*)'WARNING: Increase dimension in readro'
                      stop
             endif
      
             read(ifile)(ron(ii),ii=1,ndgg)
      
      
             close(ifile)
            return
            end subroutine readro2
      
      !***********************************************************************
      


      subroutine loadsp(levn,levp,myid)

      use input_sp

      implicit double precision (a-h,o-z)

      include 'types_ndgi_int.inc'
      include 'formats_ndgi_int.inc'
      include 'input_ndgi_int.inc'

      type(level_typ),dimension(:), allocatable :: levn,levp


!      write(*,*)'Loading of input '


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

      jmin=0
      jmax=0
      do i=1,ippmx
       if (levp(i)%j.gt.jmax) jmax=levp(i)%j
      enddo

      do i=1,ipnmx
       if (levn(i)%j.gt.jmax) jmax=levn(i)%j
      enddo

      end subroutine loadsp

     
      end module phoninteracn 
      
      
