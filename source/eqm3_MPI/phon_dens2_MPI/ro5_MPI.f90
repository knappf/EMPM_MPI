!     last update 7.2 .2020
module rdens

include 'types_phon_dens.inc'


 contains 

!************************************************************************
      subroutine roo(jamax,jbmax,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax,iphous,iphous2,phonbs,nphon,ns1,ns2,myid,numprocs)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

!      include 'types_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (roc_typ), dimension(:), allocatable :: roc
      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,iphous2,ipozl
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:,:), allocatable :: rog

!     double precision, dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph

      real(kind=4), dimension(:,:,:,:,:), allocatable :: rpp,rnp,rph,rnh

      character(len=30) fnamex,fnamec
      character(len=30) fnamer1,fnamer2,fnamer,fnameo
      character(len=10) fnamern,name_new
      integer*8 ndimroc,iroc

      integer myid,numprocs
      integer, dimension(:), allocatable :: ig_resh


      open(23,file='2phonon/2f_rnh.dat',status='unknown',form='unformatted') 
      open(32,file='2phonon/2f_rnp.dat',status='unknown',form='unformatted')
      open(34,file='2phonon/2f_rph.dat',status='unknown',form='unformatted')
      open(43,file='2phonon/2f_rpp.dat',status='unknown',form='unformatted')

      jjmx=jamax

      allocate (csixj(0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx))

      csixj=0.d0

      if (myid.eq.0) write(*,*)'CG coef. for jmx =',jjmx

      do i1=0,jjmx
!       write(*,*)'i1 =',i1
       do i2=0,jjmx
        do i3=0,jjmx
         do i4=0,jjmx
          do i5=0,jjmx
           do i6=0,jjmx             
             csixj(i1,i2,i3,i4,i5,i6)=dfloat((-1)**(i1))*sixj(2*i1,2*i2,2*i3,2*i4,2*i5,2*i6)  ! redefined phase !!!!
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      ndimroc=2000000000
      allocate(roc(ndimroc))

      roc%rho=0.0
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

      allocate(rog(0:isimax,nphon(2),nphon(1)))
      rog=0.d0

      allocate(ipozl(nphon(2)))
      ipozl=0

      call readc2_part(fnamex,xx,ns2,nphon(1),nphon(1),iphous2,ipozl)

      write(*,*)'Proces ',myid, ' loading 1-phonon densities'

!      call readro1_all(ro1_all,nphon(1),ihpmn,ipnmx,ihpmn,ipnmx,isi_mx)
      call  readro1_all(rpp,rnp,rph,rnh,nphon(1),ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax)

!      do ipa=-1,1,2
!       do ja=0,jamax

!        write(*,*)'Ipar, J ',ipa,ja

!        call readx(fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)


! MPI paralelization 
      allocate(ig_resh(ns1))
      do i=1,ns1
       ig_resh(i)=i
      enddo

      call rperm(ns1, ig_resh)

!      do ig=1,ifmx
!       do ig=ig_min,ig_max

      if (mod(ns1,numprocs).eq.0) then
              n_seg=ns1/numprocs
      else
              n_seg=ns1/numprocs+1
      endif

      if (myid.eq.0) write(*,*) ' size of segment = ',n_seg


!        do ia_cal=1,ns1 !nphon(2)

       do ia_cal=myid*n_seg+1,min((myid+1)*n_seg,ns1)
         
          
          iaaa=iphous(ig_resh(ia_cal))

          ja=phonbs(nff,iaaa)%jj

          write(*,*) ' Process #',myid,'  calculating ia = ',iaaa, 'Phonon energy ',phonbs(nff,iaaa)%enf


!         iaaa=iphous(ia_cal)

          call readx_ia(iaaa,fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)


       enerf=phonbs(nff,iaaa)%enf

       
       iroc=1

!        if (iphous(iaaa).eq.1) then 


        ibb=ilac

      do il=1,nphon(1) ! lambda'

             rog=0.d0
             jilap=phonbs(1,il)%jj

             do i=1,nbb

            il1=cc(1,i)%is 
            ila=cc(1,i)%ig
            jila=phonbs(1,ila)%jj
            jil1=phonbs(1,il1)%jj

            ifaz=(-1)**(jilap+ja+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0 

            do ibbb=1,nphon(nff)

               jb=phonbs(nff,ibbb)%jj

                  isi_mn=iabs(ja-jb)
                  isi_mx=ja+jb
            
               rog(isi_mn:isi_mx,ibbb,ila)=rog(isi_mn:isi_mx,ibbb,ila)+fact*cc(1,i)%am*xx(ibbb,il1,il)*csixj(isi_mn:isi_mx,jilap,jila,jil1,ja,jb)

               ifaz2=(-1)**(jila+jb+jil1)
               fact2=dfloat(ifaz2)*(dfloat(2*ja+1))**0.5d0
                     
               rog(isi_mn:isi_mx,ibbb,il1)=rog(isi_mn:isi_mx,ibbb,il1)+fact2*cc(1,i)%am*xx(ibbb,il,ila)*csixj(isi_mn:isi_mx,jilap,jil1,jila,ja,jb)
        
            enddo   ! over ibbb
        
        !            enddo  ! isi
              
            enddo ! over i
        
        !            endif


            do ila=1,nphon(1)
                  do ibbb=1,nphon(nff)
                        do isi=0,isimax
           
                        if (dabs(rog(isi,ibbb,ila)).gt.1.d-8) then 
                          roc(iroc)%rho=rog(isi,ibbb,ila)       
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


!         endif
          
        if (iroc.ne.1) write(*,*)'iroc =', iroc
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

!       enddo

!      enddo

      return
      end subroutine roo

!************************************************************************
      subroutine ropar3(imin,imax,ityp,iph,jamax,jbmax,isimax,rom,irom,ifile,phonbs,nphon,fnamern,iaaa,name_new,rnp,rpp,rnh,rph)

            use anglib   ! angular momentum staff
      
            implicit double precision(a-h,o-z)
      
!            include 'types_phon_dens.inc'
            include 'input_phon_dens.inc'
            include 'formats_phon_dens.inc'
      
      
            type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
            type (roc_typ), dimension(:), allocatable :: roc,rom
            type (ro_real_typ), dimension(:), allocatable :: rh
            type(rho_typ), dimension(:), allocatable :: rop
      
            integer, dimension (:), allocatable :: ndx,ndc
      
            type (phon_typ), dimension (:,:), allocatable :: phonbs
            integer, dimension (:), allocatable :: nphon
            double precision, dimension (:,:,:), allocatable :: xx
            double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
            double precision, dimension (:,:,:,:,:), allocatable :: rho,rhon
            double precision, dimension (:,:,:,:), allocatable :: rr
      
      !      double precision, dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph
            real(kind=4), dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph,rxy
      
      
            character(len=30) fnamer1,fnamer2,fnamer,fnameo
      
            character(len=10) fnamern,name_new
            character(len=2) denst
            integer*8 irom,ii
      
      
            character*6 nlam
              if (ityp.eq.-1.and.iph.eq.1) denst='pp'
              if (ityp.eq.-1.and.iph.eq.-1) denst='ph'
              if (ityp.eq.1.and.iph.eq.1) denst='np'
              if (ityp.eq.1.and.iph.eq.-1) denst='nh'
      
             if (.not.allocated(rxy)) allocate(rxy(imin:imax,imin:imax,0:isimax,nphon(1),nphon(1))) !
      
             if (ityp.eq.-1.and.iph.eq.1) rxy=rpp
             if (ityp.eq.-1.and.iph.eq.-1) rxy=rph
             if (ityp.eq.1.and.iph.eq.1) rxy=rnp
             if (ityp.eq.1.and.iph.eq.-1) rxy=rnh
      
      
            iirg=0
      
            if (irom.ne.0) then
      
      !      isimx=0
            nff=2
            ndla=nphon(nff-1)
      
            if (.not.allocated(rr)) allocate(rr(imin:imax,imin:imax,0:isimax,nphon(2)))
            rr=0.d0
      
             ndrh=100000000
             if (.not.allocated(rh)) allocate(rh(ndrh))
      
              iirg=0
      
      !        do i1=imin,imax
               do ii=1,irom
      
               isir=rom(ii)%is
               iba=rom(ii)%ib
               ilam=rom(ii)%i1
               ilamp=rom(ii)%i2
      
      !          do i_2=imin,imax
      !          do i_1=imin,imax
      
      !         if (ityp.eq.-1.and.iph.eq.1) rrr=rpp(i_1,i_2,ilam,ilamp)
      !         if (ityp.eq.-1.and.iph.eq.-1) rrr=rph(i_1,i_2,ilam,ilamp)
      !         if (ityp.eq.1.and.iph.eq.1) rrr=rnp(i_1,i_2,ilam,ilamp)
      !         if (ityp.eq.1.and.iph.eq.-1) rrr=rnh(i_1,i_2,ilam,ilamp)
      
      !         if (i1.eq.23.and.i2.eq.23.and.ityp.eq.1.and.iph.eq.1) then
      !          write(771,'(3i5,3f15.10)')iba,ilam,ilamp,rr,rom(ii)%rho,rr(iba,isir)
      !         endif
      !         rr(iba,i_1,i_2)=rr(iba,i_1,i_2)+rrr*rom(ii)%rho
      
      !          rr(i_1,i_2,iba)=rr(i_1,i_2,iba)+rrr*rom(ii)%rho
!                write(*,*)imin,imax,isir,ilam,ilamp,iba
                rr(imin:imax,imin:imax,isir,iba)=rr(imin:imax,imin:imax,isir,iba)+rxy(imin:imax,imin:imax,isir,ilam,ilamp)*rom(ii)%rho
      !           rr(imin:imax,imin:imax,isir,iba)=rxy(imin:imax,imin:imax,isir,ilam,ilamp)*rom(ii)%rho
      
      
      !          enddo        ! i_1
      !          enddo        ! i_2
      
              enddo ! over ii
      
              iirg=0
               do iba=1,nphon(nff)
                do isi=0,isimax
                 do i_2=imin,imax
                  do i_1=imin,imax
      
                if (dabs(rr(i_1,i_2,isi,iba)).gt.1.d-10) then
      
                iirg=iirg+1
                if (iirg.gt.ndrh) then
                 write(*,*)'Increase ndrh in ropar2'
                 stop
                endif
                rh(iirg)%ib=iba
                rh(iirg)%is=int(isi,1)
                rh(iirg)%i1=int(i_1,2)
                rh(iirg)%i2=int(i_2,2)
                rh(iirg)%rho=real(rr(i_1,i_2,isi,iba))
      !          if (iba.eq.iaaa.and.i1.eq.i2.and.isir.eq.0) write(937,'(a5,5i5,f15.10)')denst,iaaa,iba,i_1,i_2,isir,rh(iirg)%rho
                endif
      
                   enddo ! i_1
                   enddo ! i_2
                enddo
               enddo  ! iba
      
            endif ! irom
      
      !      write(ifile)iaaa,iirg
            if (iirg.gt.0) then
      !      write(ifile)(rh(iii)%ib,iii=1,iirg)
      !      write(ifile)(rh(iii)%is,iii=1,iirg)
      !      write(ifile)(rh(iii)%i1,iii=1,iirg)
      !      write(ifile)(rh(iii)%i2,iii=1,iirg)
      !      write(ifile)(rh(iii)%rho,iii=1,iirg)
      !      do iii=1,iirg
      !       write(9991,'(5i8,f15.10)')iaaa,rh(iii)%ib,rh(iii)%is,rh(iii)%i1,rh(iii)%i2,rh(iii)%rho
      !      enddo
      
           write(nlam,'(i6.6)')iaaa
      !    write(nlam,'(i6.6)')iaaa
      
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
      
            deallocate(rxy)
            return
            end subroutine ropar3
      
      !******************************************************************************

     

      subroutine loadphon(phonbs,nphon)

      implicit double precision (a-h,o-z)

!      include 'types_phon_dens.inc'

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

!      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous


      allocate(iphous(nphon(nf)))
      iphous=0


      xthr_min=0.0d0
      xthr_max=10000000000000000.0d0

      ii=0
      do i=1,nphon(nf)
       xene=phonbs(nf,i)%enf
       jjf=phonbs(nf,i)%jj

       if (xene.le.xthr_max.and.xene.gt.xthr_min) then
        ii=ii+1
        iphous(ii)=i
       endif
      enddo

      write(*,*)'Energy threshold for 2 phon. dens.',xthr_min,xthr_max
      write(*,*)' Number of selected phonons a)',ii

      ns1=ii

      end subroutine selphon

!**********************************************************************
      subroutine selphon2(nf,phonbs,nphon,iphous,ns2)

      implicit double precision (a-h,o-z)

!      include 'types_phon_dens.inc'

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

!      include 'types_phon_dens.inc'

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
      subroutine readx_ia(ia,fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)

!      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0
      ndimj=100000
      ii=0

      if (allocated(xcc)) deallocate(xcc)

      do while (.not.eof(2))

        read(2)ipar,ijj,no,idphon
        ii=ii+no

        if (ia.gt.ii) then 

          do ilam=1,no
            read(2)
          enddo

        else


          allocate(xcc(1,no))

          do ilam=1,ia-ii+no-1
            read(2)
          enddo
          
          read(2)(xcc(1,i)%ig,xcc(1,i)%is,xcc(1,i)%am,i=1,idphon)
          close(2) 
          return 

        endif 
      enddo


      return
      end subroutine readx_ia


      subroutine readrho(fname,rho,nphon,isimx,imin,imax)

      implicit double precision (a-h,o-z)

!      include 'types_phon_dens.inc'

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
!      include 'types_phon_dens.inc'

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


!!******************************************************************************
      subroutine readro1_all(rpp,rnp,rph,rnh,ndla,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax)

            implicit double precision (a-h,o-z)
      
      !      include 'types_eqm.inc'
      !     include 'types_phon_dens.inc'
      
      !      double precision, dimension(:,:,:,:,:), allocatable :: rpp,rnp,rph,rnh
            real(kind=4), dimension(:,:,:,:,:), allocatable :: rpp,rnp,rph,rnh
            type(rho_typ), dimension(:), allocatable :: ronn
            logical je_tam_subor
      
      
            character(len=10)fname
            character(len=4)nlam
      
            if (.not.allocated(rpp)) allocate(rpp(ippmn:ippmx,ippmn:ippmx,0:isimax,ndla,ndla))
            rpp=0.0
            if (.not.allocated(rnp)) allocate(rnp(ipnmn:ipnmx,ipnmn:ipnmx,0:isimax,ndla,ndla))
            rnp=0.0
            if (.not.allocated(rph)) allocate(rph(ihpmn:ihpmx,ihpmn:ihpmx,0:isimax,ndla,ndla))
            rph=0.0
            if (.not.allocated(rnh)) allocate(rnh(ihnmn:ihnmx,ihnmn:ihnmx,0:isimax,ndla,ndla))
            rnh=0.0
      
            ndro=5000000
            if (.not.allocated(ronn)) allocate (ronn(ndro))
      
      !      return
      
            ifile=33
      
            do ityp=1,4
             if (ityp.eq.1) fname='1f_rnp.dat'
             if (ityp.eq.2) fname='1f_rpp.dat'
             if (ityp.eq.3) fname='1f_rnh.dat'
             if (ityp.eq.4) fname='1f_rph.dat'
      
            ifile=33
      
            do ig=1,ndla
      
      !       write(*,*)'Read ig',fname,ig
      
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
      !       if (isit.eq.0) then
              if (ityp.eq.1) rnp(i1t,i2t,isit,ig,ibt)=real(rot)
              if (ityp.eq.2) rpp(i1t,i2t,isit,ig,ibt)=real(rot)
              if (ityp.eq.3) rnh(i1t,i2t,isit,ig,ibt)=real(rot)
              if (ityp.eq.4) rph(i1t,i2t,isit,ig,ibt)=real(rot)
      
      !       endif
             enddo
             enddo
             enddo
      
            return
            end subroutine readro1_all
      
      
      !******************************************************************************

      subroutine readc2_part(fname,cc,n1,n2,n3,iphous2,ipozl)

      implicit double precision (a-h,o-z)

!      include 'types_phon_dens.inc'

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
!        xxr=real(xcc(1,i)%am)
!        cc(illl,ii,jj)=xxr
        cc(illl,ii,jj)=xcc(1,i)%am
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
!-----------------------------------------------------------------
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
!--------------------------------------------------------------


end module rdens
      
