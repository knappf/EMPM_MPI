  module admatr

   contains

   subroutine admat(nf,ipar,jcal,phonbs,ndim,no,phon1,phon2,mxt)

   use anglib
   use input_sp

   implicit double precision (a-h,o-z)


   include 'types_eqm.inc'

   integer, dimension (:), allocatable :: phonus, phonmus
   type(phon_typ), dimension (:), allocatable :: phon1,phon2
   type(phonbase_typ), dimension (:), allocatable :: phonbs
   type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
   type(level_typ),dimension(:), allocatable :: levn,levp

   double precision, dimension(:,:,:,:),allocatable :: ronp1,ropp1,ronh1,roph1

   double precision, dimension(:,:,:,:), allocatable :: fnp,fpp,fnh,fph

   double precision, dimension(:,:,:,:,:), allocatable :: csixj,csixjd

   integer, dimension(:), allocatable :: nx,mxt,mxtr

   integer, dimension (:,:,:), allocatable :: ipozbr
   integer, dimension (:,:), allocatable :: ndbr
   double precision, dimension (:,:), allocatable :: rtdaovr
   double precision, dimension (:), allocatable :: a_m,d_m


   character*30 namernp,namerpp,namernh,namerph,namefnp,namefpp,namefnh,namefph,namernp1,namerpp1,namernh1,namerph1

    allocate(rtdaovr(5000,5000))
    rtdaovr=0.0d0
    open(881,file='tda_r_overl.dat',status='unknown',form='formatted')
     do while (.not.eof(881))
       read(881,*)iii,eee,rr
     enddo
      close(881)

      ilamcm=iii
      write(*,*)' CM phonon is number ',ilamcm

     open(881,file='tda_r_overl.dat',status='unknown',form='formatted')
      do while (.not.eof(881))
       read(881,*)iii,eee,rr
       rtdaovr(iii,ilamcm)=rr
       rtdaovr(ilamcm,iii)=rr
      enddo
     close(881)

    open(2,file='mxtr.dat',status='unknown',form='unformatted')
     read(2)idphontr
      allocate(mxtr(idphontr))
     read(2)(mxtr(i),i=1,idphontr)
    close(2)
 
     no=idphontr

      open(6,file='a_mat.dat',status='unknown',form='unformatted')
      open(7,file='d_mat.dat',status='unknown',form='unformatted')


      call loadsp(levn,levp)

      namefnp='V_phon_p_n.dat'
      namefpp='V_phon_p_p.dat'
      namefnh='V_phon_h_n.dat'
      namefph='V_phon_h_p.dat'

      if (nf.eq.2) then
      namerpp='1phonon/1f_rpp.dat'
      namernp='1phonon/1f_rnp.dat'
      namerph='1phonon/1f_rph.dat'
      namernh='1phonon/1f_rnh.dat'

      namernp1='1phonon/1f_rnp1.dat'
      namerpp1='1phonon/1f_rpp1.dat'
      namernh1='1phonon/1f_rnh1.dat'
      namerph1='1phonon/1f_rph1.dat'

      endif

      open(33,file=namernp,status='old',form='unformatted')
      open(34,file=namerpp,status='old',form='unformatted')
      open(43,file=namernh,status='old',form='unformatted')
      open(44,file=namerph,status='old',form='unformatted')

      open(331,file=namefnp,status='old',form='unformatted')
      open(341,file=namefpp,status='old',form='unformatted')
      open(431,file=namefnh,status='old',form='unformatted')
      open(441,file=namefph,status='old',form='unformatted')

      open(332,file=namernp1,status='old',form='unformatted')
      open(342,file=namerpp1,status='old',form='unformatted')
      open(432,file=namernh1,status='old',form='unformatted')
      open(442,file=namerph1,status='old',form='unformatted')

      call read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max)

!      write(*,*)ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max
      ndla=ndla_max
      nsi=isi_max

!      n1mn=1
!      n1mx=45
!      n2mn=1
!      n2mx=45

      iaaold=-10
      iaold=-10


      jmx=nsi !20

      allocate(csixj(0:jmx,0:jmx,0:jmx,0:jmx,0:jmx))
      csixj=0.0d0

      do i=0,jmx
       do j=0,jmx
        do k=0,jmx
         do l=0,jmx
          do m=0,jmx
           csixj(i,j,k,l,m)=sixj(2*jcal,2*i,2*j,2*k,2*l,2*m)
          enddo
         enddo
        enddo
       enddo
      enddo

     allocate(csixjd(0:jmx,0:jmx,0:jmx,0:jmx,0:jmx))
      csixjd=0.0d0

      do i=0,jmx
       do j=0,jmx
        do k=0,jmx
         do l=0,jmx
          do m=0,jmx
           csixjd(i,j,k,l,m)=sixj(2*i,2*j,2*k,2*jcal,2*l,2*m)
          enddo
         enddo
        enddo
       enddo
      enddo

      write(6)idphontr,ndim
      write(7)idphontr,ndim

      allocate(a_m(ndim),d_m(ndim))

      do i=1,idphontr

       iii=mxtr(i)
       iaa=phonbs(iii)%ilap  ! 1phonon index
       ia=phonbs(iii)%ila    ! n-1 phonon index
       write(912,*)' **** ',i,iaa,ia

      enddo



      do i=1,idphontr
       a_m=0.d0
       d_m=0.d0

       iii=mxtr(i)
       write(725,*)' **** ',i,idphontr
       iaa=phonbs(iii)%ilap  ! 1phonon index
       ia=phonbs(iii)%ila    ! n-1 phonon index


       jila=phon1(iaa)%j
       jial=phon2(ia)%j


       if (ia.ne.iaold) then
             call readro(33,ia,ronp,nronp)
             call readro(34,ia,ropp,nropp)
             call readro(43,ia,ronh,nronh)
             call readro(44,ia,roph,nroph)
             iaold=ia

            call redrsum(ndla,ronp,nronp,ropp,nropp,ronh,nronh,roph,nroph,ipozbr,ndbr)


        endif

       if (iaa.ne.iaaold) then

        call readfn(331,iaa,fnp,ndla,ipn_min,ipn_max,ipn_min,ipn_max,nsi)
        call readfn(341,iaa,fpp,ndla,ipp_min,ipp_max,ipp_min,ipp_max,nsi)
        call readfn(431,iaa,fnh,ndla,ihn_min,ihn_max,ihn_min,ihn_max,nsi)
        call readfn(441,iaa,fph,ndla,ihp_min,ihp_max,ihp_min,ihp_max,nsi)
!        iaaold=iaa
!       endif
!  
!       if (iaa.ne.iaaold) then

        call readro11(332,iaa,ronp1,ndla,ipn_min,ipn_max,ipn_min,ipn_max,nsi)
        call readro11(342,iaa,ropp1,ndla,ipp_min,ipp_max,ipp_min,ipp_max,nsi)
        call readro11(432,iaa,ronh1,ndla,ihn_min,ihn_max,ihn_min,ihn_max,nsi)
        call readro11(442,iaa,roph1,ndla,ihp_min,ihp_max,ihp_min,ihp_max,nsi)

        iaaold=iaa
       endif

       call OMP_SET_NUM_THREADS(5)
!$omp parallel default(shared) private(jjj,dd,aa, &
!$omp xsixj,ii,iii,ib,ibb,jlap,jilap,jialp,ibt,i1,i2, &
!$omp isi,faz,itid,ji1,ji2)

!$omp do schedule (dynamic,2000)


       do j=1,ndim

        jjj=j

        aa=0.d0
        dd=0.d0

        ib=phonbs(jjj)%ila
        ibb=phonbs(jjj)%ilap

        jilap=phon1(ibb)%j
        jialp=phon2(ib)%j

        if (ia.eq.ib.and.iaa.eq.ibb) aa=aa+phon1(iaa)%enf+phon2(ia)%enf
!        if (ipar.eq.-1.and.jcal.eq.1) then

!          if (ia.eq.ib.and.iaa.eq.ilamcm) then
!              dd=dd+rtdaovr(ibb)
!          endif

          if (ia.eq.ib) then
!              write(*,*)'***'
             aa=aa+rtdaovr(iaa,ibb)
          endif
          if (iaa.eq.ibb) then
!              write(*,*)'***'
             aa=aa+rtdaovr(ia,ib)
          endif

        if (nf.eq.2) then

          if (ia.eq.ib.and.iaa.eq.ibb) dd=dd+1.d0
          jlap=phon1(ibb)%j
          jla=phon1(iaa)%j
!          xsixj=sixj(2*jla,0,2*jla,2*jlap,2*jcal,2*jlap)
          xsixj=csixjd(0,jla,jla,jlap,jlap)
!c          faz=dfloat((-1)**((ji1-ji2)/2))
          if ((ia.eq.ibb).and.(ib.eq.iaa)) then
            dd=dd+xsixj*(dfloat(2*jla+1)*(2*jlap+1))**0.5d0
          endif
        endif

         faz1=-1.d0*dfloat((-1)**(jcal+phon1(ibb)%j+phon2(ia)%j))



        do iii=1,ndbr(1,ib)  !nronp          ! neutron particle
         ii=ipozbr(1,ib,iii)
         ibt=ronp(ii)%ilap
         if (ibt.eq.ib) then
          i1=ronp(ii)%i1
          i2=ronp(ii)%i2
          isi=ronp(ii)%j
          ji1=levn(i1)%j
          ji2=levn(i2)%j

!          xsixj=sixj(2*jcal,2*jilap,2*jialp,2*isi,2*jial,2*jila)
          xsixj=csixj(jilap,jialp,isi,jial,jila)
          faz=dfloat((-1)**(jcal+jial))
          aa=aa+xsixj*faz*ronp(ii)%ro*fnp(ibb,isi,i1,i2)

          xsixj=csixjd(isi,jilap,jila,jial,jialp)
          faz=dfloat((-1)**((ji1-ji2)/2))
          dd=dd+faz1*faz*xsixj*ronp(ii)%ro*ronp1(ibb,isi,i2,i1)

!          write(993,*)xsixj,faz,ronp(ii)%ro,fnp(ibb,isi,i1,i2)
         endif
        enddo

        do iii=1,ndbr(2,ib)        !nropp          ! proton particle
         ii=ipozbr(2,ib,iii)
         ibt=ropp(ii)%ilap
         if (ibt.eq.ib) then
          i1=ropp(ii)%i1
          i2=ropp(ii)%i2
          ji1=levp(i1)%j
          ji2=levp(i2)%j

          isi=ropp(ii)%j
!          xsixj=sixj(2*jcal,2*jilap,2*jialp,2*isi,2*jial,2*jila)
          xsixj=csixj(jilap,jialp,isi,jial,jila)
          faz=dfloat((-1)**(jcal+jial))
          aa=aa+xsixj*faz*ropp(ii)%ro*fpp(ibb,isi,i1,i2)
!          write(993,*)xsixj,faz,ropp(ii)%ro,fpp(ibb,isi,i1,i2)


          xsixj=csixjd(isi,jilap,jila,jial,jialp)
          faz=dfloat((-1)**((ji1-ji2)/2))
          dd=dd+faz1*faz*xsixj*ropp(ii)%ro*ropp1(ibb,isi,i2,i1)

         endif
        enddo

        do iii=1,ndbr(3,ib) ! nronh          ! neutron hole
         ii=ipozbr(3,ib,iii)
         ibt=ronh(ii)%ilap
         if (ibt.eq.ib) then
          i1=ronh(ii)%i1
          i2=ronh(ii)%i2

          ji1=levn(i1)%j
          ji2=levn(i2)%j

          isi=ronh(ii)%j
!          xsixj=sixj(2*jcal,2*jilap,2*jialp,2*isi,2*jial,2*jila)
          xsixj=csixj(jilap,jialp,isi,jial,jila)

          faz=dfloat((-1)**(jcal+jial))
          aa=aa+xsixj*faz*ronh(ii)%ro*fnh(ibb,isi,i1,i2)
!          write(993,*)xsixj,faz,ronh(ii)%ro,fnh(ibb,isi,i1,i2)

          xsixj=csixjd(isi,jilap,jila,jial,jialp)
          faz=dfloat((-1)**((ji2-ji1)/2))
          dd=dd+faz1*faz*xsixj*ronh(ii)%ro*ronh1(ibb,isi,i2,i1)

         endif
        enddo

!        write(993,*)'**************'

        do iii=1,ndbr(4,ib)  !nroph          ! proton hole
         ii=ipozbr(4,ib,iii)
         ibt=roph(ii)%ilap
         if (ibt.eq.ib) then
          i1=roph(ii)%i1
          i2=roph(ii)%i2
          ji1=levp(i1)%j
          ji2=levp(i2)%j

          isi=roph(ii)%j
!          xsixj=sixj(2*jcal,2*jilap,2*jialp,2*isi,2*jial,2*jila)
          xsixj=csixj(jilap,jialp,isi,jial,jila)

          faz=dfloat((-1)**(jcal+jial))
          aa=aa+xsixj*faz*roph(ii)%ro*fph(ibb,isi,i1,i2)
!          write(993,*)xsixj,faz,roph(ii)%ro,fph(ibb,isi,i1,i2)

          xsixj=csixjd(isi,jilap,jila,jial,jialp)

          faz=dfloat((-1)**((ji2-ji1)/2))
          dd=dd+faz1*faz*xsixj*roph(ii)%ro*roph1(ibb,isi,i2,i1)

         endif
        enddo

!        if (dabs(dd).gt.1.d-10) write(999,*)i,j,dd
!        if (dabs(aa).gt.1.d-10) write(6)i,j,aa
!        if (dabs(dd).gt.1.d-10) write(7)i,j,dd
       
        a_m(j)=aa
        d_m(j)=dd

       enddo  ! over j

!$omp end do

!$omp end parallel

       
       write(6)(a_m(jjj),jjj=1,ndim)
       write(7)(d_m(jjj),jjj=1,ndim)



      enddo  ! over i


      deallocate(a_m,d_m)


      close(33)
      close(34)
      close(43)
      close(44)
      close(331)
      close(341)
      close(431)
      close(441)
      close(332)
      close(342)
      close(432)
      close(442)

      close(6)
      close(7)



      return

      end subroutine admat







!*****************************************************************************
   subroutine readro(ifile,ig,ron,ndgg)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
      include 'types_eqm.inc'

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=30)fname

      rewind(ifile)

      ndro=5000000
      ndgg=0

      if (.not.allocated(ron)) allocate (ron(ndro))

       do igt=1,ig-1
         read(ifile)igggt,nggt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)xtt
       enddo


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


      return
      end subroutine readro

!***********************************************************************
      subroutine redrsum(ifonmx,ronp,nronp,ropp,nropp,ronh,nronh,roph,nroph,ipozbr,ndbr)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
      include 'types_eqm.inc'

      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      integer, dimension (:,:,:), allocatable :: ipozbr
      integer, dimension (:,:), allocatable :: ndbr

      character(len=30)fname

      ndipo=50000

      if (.not.allocated(ipozbr)) allocate(ipozbr(4,ifonmx,ndipo))
      ipozbr=0

      if (.not.allocated(ndbr)) allocate(ndbr(4,ifonmx))
      ndbr=0



      do i=1,nronp
       ibt=ronp(i)%ilap
       ndbr(1,ibt)=ndbr(1,ibt)+1

        if (ndbr(1,ibt).gt.ndipo) then
               write(*,*)' Increase dimension in redrsum'
               stop
            endif
       ipozbr(1,ibt,ndbr(1,ibt))=i

      enddo

      do i=1,nropp
       ibt=ropp(i)%ilap
       ndbr(2,ibt)=ndbr(2,ibt)+1

        if (ndbr(2,ibt).gt.ndipo) then
               write(*,*)' Increase dimension in redrsum'
               stop
            endif
       ipozbr(2,ibt,ndbr(2,ibt))=i
      enddo

      do i=1,nronh
       ibt=ronh(i)%ilap
       ndbr(3,ibt)=ndbr(3,ibt)+1

        if (ndbr(3,ibt).gt.ndipo) then
               write(*,*)' Increase dimension in redrsum'
               stop
            endif
       ipozbr(3,ibt,ndbr(3,ibt))=i

      enddo

     do i=1,nroph
       ibt=roph(i)%ilap
       ndbr(4,ibt)=ndbr(4,ibt)+1

        if (ndbr(4,ibt).gt.ndipo) then
               write(*,*)' Increase dimension in redrsum'
               stop
            endif
       ipozbr(4,ibt,ndbr(4,ibt))=i
      enddo



      return
      end subroutine redrsum
!***********************************************************************


      subroutine readro11(ifile,ig,ron,ndla,n1mn,n1mx,n2mn,n2mx,nsi)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      double precision, dimension(:,:,:,:), allocatable :: ron
      type(rho_typ), dimension(:), allocatable :: ronn



      character(len=30)fname

      rewind(ifile)

      ndro=5000000
      ndgg=0

      if (.not.allocated(ronn)) allocate (ronn(ndro))

       do igt=1,ig-1
         read(ifile)igggt,nggt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)xtt
       enddo


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



      if (.not.allocated(ron)) allocate(ron(ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))

      ron=0.d0

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

        ron(ibt,isit,i1t,i2t)=rot

       enddo



      return
      end subroutine readro11


!******************************************************************************
     subroutine readfn(ifile,ig,ron,ndla,n1mn,n1mx,n2mn,n2mx,nsi)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      double precision, dimension(:,:,:,:), allocatable :: ron
      type(rho_typ), dimension(:), allocatable :: ronn

      character(len=30)fname

      rewind(ifile)

      ndro=5000000
      ndgg=0

      if (.not.allocated(ronn)) allocate (ronn(ndro))

       do igt=1,ig-1
         read(ifile)igggt,nggt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)xtt
       enddo

       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readfn',ndgg,ndro
                stop
       endif

       read(ifile)(ronn(ii)%ilap,ii=1,ndgg)
       read(ifile)(ronn(ii)%j,ii=1,ndgg)
       read(ifile)(ronn(ii)%i1,ii=1,ndgg)
       read(ifile)(ronn(ii)%i2,ii=1,ndgg)
       read(ifile)(ronn(ii)%ro,ii=1,ndgg)

      if (.not.allocated(ron)) allocate(ron(ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))

       ron=0.d0

       do ii=1,ndgg

        ibt=ronn(ii)%ilap
        isit=ronn(ii)%j
        i1t=ronn(ii)%i1
        i2t=ronn(ii)%i2
        rot=ronn(ii)%ro

        if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
         write(*,*)' Small dimensions of ron in readfn',ibt,isit,i1t,i2t
         stop
       endif

        ron(ibt,isit,i1t,i2t)=rot

       enddo

       return

       end subroutine readfn

!***********************************************************************

      subroutine readf(ifile,ig,ron,ndla,n1mn,n1mx,n2mn,n2mx,nsi)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      double precision, dimension(:,:,:,:), allocatable :: ron

      character(len=30)fname


      if (.not.allocated(ron)) allocate(ron(ndla,0:nsi,n1mn:n1mx,n2mn:n2mx))


  111 continue
       ibt=1
       noig=0

       ron=0.d0

  112  read(ifile)igg

       if(igg.eq.0) goto 112


       do while (ibt.ne.0)
       read(ifile)ibt,isit,i1t,i2t,rot

       if (igg.eq.ig) then

       if (ibt.ne.0) then
       if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
         write(*,*)' Small dimensions of ron in read1',ibt,isit,i1t,i2t
         stop
       endif

        ron(ibt,isit,i1t,i2t)=rot
       endif

       endif

       enddo

        if (ig.ne.igg) then
!c            write(*,*)'WARNING: Loaded ig does not match '
            if (igg.gt.ig) rewind(ifile)
            goto 111
       endif

      return
      end subroutine readf
!***********************************************************************


      subroutine loadsp(levn,levp)

      use input_sp

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'
      include 'formats_eqm.inc'

      type(level_typ),dimension(:), allocatable :: levn,levp


!c      write(*,*)'Loading of input '


      open(1,file='input_eqm_coup.dat',status='old',form='formatted')

      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      read(1,26)alfa,beta
      read(1,*)
      read(1,15)iparmn,iparmx
      read(1,15)jminn,jmaxn
      close(1)

      allocate(levn(ipnmx+1000),levp(ippmx+1000))


      call inp_sp(levn,levp)

      end subroutine loadsp
!
      subroutine  read_input_par(ihp_min,ihp_max,ihn_min,ihn_max,ipp_min,ipp_max,ipn_min,ipn_max,isi_max,ndla_max)


      include 'formats_eqm.inc'

      character*30 name1f


      open(1,file='input_eqm_coup.dat',status='old',form='formatted')

      read(1,15)ia,iz
      read(1,15)ihn_min,ihn_max
      read(1,15)ihp_min,ihp_max
      read(1,15)ipn_min,ipn_max
      read(1,15)ipp_min,ipp_max
      read(1,26)alfa,beta
      read(1,*)
      read(1,15)iparmn,iparmx
      read(1,15)jminn,jmaxn
      close(1)

      isi_max=jmaxn

      name1f='1phonon/1f_states.dat'


      open (3,file=name1f,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
      enddo

      ndla_max=i

      close(3)


      end subroutine read_input_par

      end module admatr

