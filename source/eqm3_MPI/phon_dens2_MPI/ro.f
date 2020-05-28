!     last update 30.6.2010

      module rdens

      contains 

      subroutine rop(jamax,jbmax,isimax)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
      include 'input_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (ro_typ), dimension(:), allocatable :: roc
      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj


      character(len=30) fnamex,fnamec

      open(22,file='ro1.dat',status='unknown',form='unformatted')

      jjmx=12

      allocate (csixj(0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx))

      csixj=0.d0

      do i1=0,jjmx
       write(*,*)'i1 =',i1
       do i2=0,jjmx
        do i3=0,jjmx
         do i4=0,jjmx
          do i5=0,jjmx
           do i6=0,jjmx             
             csixj(i1,i2,i3,i4,i5,i6)=
     *sixj(2*i1,2*i2,2*i3,2*i4,2*i5,2*i6)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      ndimroc=10000000
      allocate(roc(ndimroc))

      roc%rho=0.d0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0

      call loadphon(phonbs,nphon)

      ndamp=10000
      ndrho=1000000     
      nff=2

      if (nff.eq.2) then 
       fnamex='2phonon/2f_x.dat'
       fnamec='2phonon/2f_c.dat'
      endif

      write(*,*)' Number of 1phonon states ',nphon(1)
      write(*,*)' Number of 1phonon states ',nphon(2)



c      call readc(fnamec,cc,ndc)
 
      call readc2(fnamex,xx,nphon(2),nphon(1),nphon(1))

      do ipa=1,-1,-2
       do ja=0,jamax

        write(*,*)ipa,ja

        call readx(fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)

        do iaaa=ilax+1,ilax+nb
        write(998,*)' ********* ',iaaa-1,iroc

        iroc=1

        write(*,*)'iaaa',iaaa,iroc

        ibb=ilac

         do ibbb=1,nphon(nff)

          jb=phonbs(nff,ibbb)%jj

          do isi=0,isimax
           do il=1,nphon(1)

             do i=1,nbb


            ila=cc(iaaa,i)%is
            iga=cc(iaaa,i)%ig

            jila=phonbs(nff-1,ila)%jj
            
            jilap=phonbs(nff-1,il)%jj
            jigi=phonbs(nff-1,iga)%jj 

            ifaz=(-1)**(isi+jila+jb+jigi)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0

             
       roc(iroc)%rho=roc(iroc)%rho+fact*
     *cc(iaaa,i)%am*xx(ibbb,il,iga)*csixj(isi,jilap,jila,jigi,ja,jb)
      
            enddo ! over i

            if (dabs(roc(iroc)%rho).gt.1.d-10) then 
             roc(iroc)%ib=ibbb
             roc(iroc)%is=isi
             roc(iroc)%i1=ila
             roc(iroc)%i2=il

c             write(998,*)iaaa,ibbb,isi,ila,il,roc(iroc)%rho

            iroc=iroc+1

            if (iroc.gt.ndimroc) then 
            write(*,*)' Increase dimension of array roc '
            stop
            endif
            endif

           enddo  !over il

           enddo   ! over isi


         enddo   !over iaaa

          
         write(*,*)'iroc', iroc
c         write(998,*)(roc(iii)%rho,iii=1,iroc)
         write(22)(roc(iii)%ib,roc(iii)%is,roc(iii)%i1,roc(iii)%i2,
     *roc(iii)%rho,iii=1,iroc-1)

          roc%rho=0.d0
          roc%ib=0
          roc%is=0
          roc%i1=0
          roc%i2=0

         continue


         enddo

       enddo

      enddo

      return
      end subroutine rop
************************************************************************
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

**********************************************************************

      subroutine readx(fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

c      allocate(xcc(ndimi,ndimj))
c      allocate(ndcc(ndimi))

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
  
       read(2)(xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%ig
     *,xcc(ilam+ilamp,i)%am,i=1,idphon)
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
************************************************************

      subroutine readc(fname,xcc,ndcc)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

c      allocate(xcc(ndimi,ndimj))
c      allocate(ndcc(ndimi))

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=20000
      ndiml=14000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(1:ndiml,ndimj))
      allocate(ndcc(1:ndiml))


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%ig
     *,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon
      enddo

      ilamps=ilamp
      ilamp=ilamp+no


c      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


c      deallocate(xcc,ndcc)
c      no=0
c      idphon=0

      enddo

  11  continue     

      close(2)

      return
      end subroutine readc
****************************************

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
  
       read(2)(xcc(1,i)%is,xcc(1,i)%ig
     *,xcc(1,i)%am,i=1,idphon)
       ndcc(1)=idphon

       do i=1,idphon
        ii=xcc(1,i)%is
        jj=xcc(1,i)%ig
        cc(ilam+ilamp,ii,jj)=xcc(1,i)%am
       enddo
      enddo


   

      ilamps=ilamp
      ilamp=ilamp+no


c      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


c      deallocate(xcc,ndcc)
c      no=0
c      idphon=0

      enddo

  11  continue     



      close(2)

      deallocate(xcc,ndcc)

      return
      end subroutine readc2
************************************************************







      end module rdens
      
