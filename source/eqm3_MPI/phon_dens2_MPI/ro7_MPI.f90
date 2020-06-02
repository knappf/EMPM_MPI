!     last update 15.5 .2020
module rdens

include 'types_phon_dens.inc'


 contains 

!************************************************************************
      subroutine  roo(jamax,jbmax,ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax,iphous,iphous2,phonbs,nphon,ns1,ns2,myid,numprocs,levn,levp,idimnp,idimnh,idimpp,idimph)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

!      include 'types_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:), allocatable :: cc!,cc
      type (roc_typ), dimension(:), allocatable :: roc
      type (ro_real_typ), dimension(:), allocatable :: r2pp,r2np,r2nh,r2ph
      type (level_typ), dimension(:),allocatable :: levn,levp
      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,iphous2,ipozl
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:), allocatable :: rog
      integer,dimension(:),allocatable :: idimnp,idimnh,idimph,idimpp
!     double precision, dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph

      real(kind=4), dimension(:,:,:,:,:), allocatable :: rpp,rnp,rph,rnh

      character(len=30) fnamex,fnamec
      character(len=30) fnamer1,fnamer2,fnamer,fnameo
      character(len=10) fnamern,name_new
      character(len=10) myid_name
      integer*8 ndimroc,iroc


      integer myid,numprocs
      integer, dimension(:), allocatable :: ig_resh


      open(23,file='2phonon/2f_rnh.dat',status='unknown',form='unformatted') 
      open(32,file='2phonon/2f_rnp.dat',status='unknown',form='unformatted')
      open(34,file='2phonon/2f_rph.dat',status='unknown',form='unformatted')
      open(43,file='2phonon/2f_rpp.dat',status='unknown',form='unformatted')


      write(myid_name,'(i10.10)')myid
      open(63,file='2_phon_dens_calc_myid_'//myid_name,status='unknown',form='formatted')
      

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

      ndimroc=200000
      allocate(roc(ndimroc))

      roc%rho=0.0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0
    


      nff=2

!      call selphon(nff,phonbs,nphon,iphous,ns1)
!      call selphon2(nff,phonbs,nphon,iphous2,ns2)

      isi_mn=0
      isi_mx=isimax



      if (nff.eq.2) then 
       fnamex='2phonon/2f_x.dat'
       fnamec='2phonon/2f_c.dat'
      endif

      allocate(rog(0:isimax,nphon(1)))
      rog=0.d0

      allocate(ipozl(nphon(2)))
      ipozl=0

!      call readc2_part(fnamex,xx,ns2,nphon(1),nphon(1),iphous2,ipozl)
      call read_ampl(fnamex,xx,ns2,nphon(1),nphon(1),iphous2,ipozl)      

      write(*,*)'Proces ',myid, ': loading 1-ph densities'


      call  readro1_all(rpp,rnp,rph,rnh,nphon(1),ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax)

      ndrh=100000000
      allocate(r2nh(ndrh),r2ph(ndrh),r2pp(ndrh),r2np(ndrh))



! MPI paralelization 
      allocate(ig_resh(ns1))
      do i=1,ns1
       ig_resh(i)=i
      enddo

      call rperm(ns1, ig_resh)


      if (mod(ns1,numprocs).eq.0) then
              n_seg=ns1/numprocs
      else
              n_seg=ns1/numprocs+1
      endif

      if (myid.eq.0) write(*,*) ' size of segment = ',n_seg



    do ia_cal=myid*n_seg+1,min((myid+1)*n_seg,ns1)
         
          ndr2pp=0
          ndr2np=0
          ndr2ph=0
          ndr2nh=0
          
          iaaa=iphous(ig_resh(ia_cal))

          ja=phonbs(nff,iaaa)%jj

       write(*,*)'Process #',myid,': calc. ia = ',iaaa, ' E_2ph=',phonbs(nff,iaaa)%enf

!       call readx_ia(iaaa,fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)
       call read_amp_ia(iaaa,fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)

       enerf=phonbs(nff,iaaa)%enf

      do ibbb=1,nphon(nff)

         jb=phonbs(nff,ibbb)%jj
!         isi_mn=iabs(ja-jb)
!         isi_mx=ja+jb

         iroc=0

      do il=1,nphon(1) ! lambda'

             rog=0.d0
             jilap=phonbs(1,il)%jj

          do i=1,nbb

            il1=cc(i)%is 
            ila=cc(i)%ig
            jila=phonbs(1,ila)%jj
            jil1=phonbs(1,il1)%jj

            ifaz=(-1)**(jilap+ja+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0 
            
            rog(0:isimax,ila)=rog(0:isimax,ila)+fact*cc(i)%am*xx(ibbb,il1,il)*csixj(0:isimax,jilap,jila,jil1,ja,jb)

            ifaz2=(-1)**(jila+jb+jil1)
            fact2=dfloat(ifaz2)*(dfloat(2*ja+1))**0.5d0
                     
            rog(0:isimax,il1)=rog(0:isimax,il1)+fact2*cc(i)%am*xx(ibbb,il,ila)*csixj(0:isimax,jilap,jil1,jila,ja,jb)
               
          enddo ! over i

            do ila=1,nphon(1)
                do isi=0,isimax
           
                      if (dabs(rog(isi,ila)).gt.1.d-12) then 
                          iroc=iroc+1    
                          roc(iroc)%rho=rog(isi,ila)       
                          roc(iroc)%ib=ibbb
                          roc(iroc)%is=isi
                          roc(iroc)%i1=ila
                          roc(iroc)%i2=il
                      
                        if (iroc.gt.ndimroc) then 
                            write(*,*)' Increase dimension of array roc '
                           stop
                        endif
                       endif
           
                 enddo
            enddo
           
      enddo  !over il

       if (iroc.ne.0) then 
         call ropar3(iaaa,ibbb,ippmn,ippmx,isimax,roc,iroc,nphon,rpp,r2pp,ndr2pp,ndrh,ja,jb,levp)  
         call ropar3(iaaa,ibbb,ipnmn,ipnmx,isimax,roc,iroc,nphon,rnp,r2np,ndr2np,ndrh,ja,jb,levn)
         call ropar3(iaaa,ibbb,ihpmn,ihpmx,isimax,roc,iroc,nphon,rph,r2ph,ndr2ph,ndrh,ja,jb,levp) 
         call ropar3(iaaa,ibbb,ihnmn,ihnmx,isimax,roc,iroc,nphon,rnh,r2nh,ndr2nh,ndrh,ja,jb,levn)  
       endif
 
      enddo

!     writting densities
      if (ndr2pp.ne.0) then 
        fnamern='1f_rpp.dat'
        name_new='2f_rpp.dat'   
        if(ndr2pp.ne.idimpp(ja)) then 
         write(*,*) 'Wrong correspondence (pp) for',iaaa,'J=',ja,'expected=',idimpp(ja),'from dens=',ndr2pp
         stop
        endif 
        call  write_dens(iaaa,ndr2pp,r2pp,fnamern,name_new)
!write(833,*)iaaa,ja,ndr2pp,idimpp(ja),name_new
!do ii=1,ndr2pp
!write(833,*)ii,r2pp(ii)%ib,r2pp(ii)%is,r2pp(ii)%i2,r2pp(ii)%i1
!write(833,*)r2pp(ii)%rho
!end do
      endif 

      if (ndr2np.ne.0) then
        fnamern='1f_rnp.dat'
        name_new='2f_rnp.dat' 
        if(ndr2np.ne.idimnp(ja)) then 
            write(*,*) 'Wrong correspondence (np) for',iaaa,'J=',ja,'expected=',idimnp(ja),'from dens=',ndr2np
            stop 
        endif 
        call  write_dens(iaaa,ndr2np,r2np,fnamern,name_new)
!write(834,*)iaaa,ja,ndr2np,idimnp(ja),name_new
!do ii=1,ndr2np
!write(834,*)r2np(ii)%ib,r2np(ii)%is,r2np(ii)%i2,r2np(ii)%i1
!write(834,*)r2np(ii)%rho
!end do
      endif  

      if (ndr2ph.ne.0) then
        fnamern='1f_rph.dat'
        name_new='2f_rph.dat' 
        if(ndr2ph.ne.idimph(ja)) then 
            write(*,*) 'Wrong correspondence (ph) for',iaaa,'J=',ja,'expected=',idimph(ja),'from dens=',ndr2ph
            stop
        endif 
        call  write_dens(iaaa,ndr2ph,r2ph,fnamern,name_new)
!write(835,*)iaaa,ja,ndr2ph,idimph(ja),name_new
!do ii=1,ndr2ph
!write(835,*)r2ph(ii)%ib,r2ph(ii)%is,r2ph(ii)%i2,r2ph(ii)%i1
!write(835,*)r2ph(ii)%rho
!end do
      endif
      
      if (ndr2nh.ne.0) then
        fnamern='1f_rnh.dat'
        name_new='2f_rnh.dat' 
        if(ndr2nh.ne.idimnh(ja)) then 
            write(*,*) 'Wrong correspondence (nh) for',iaaa,'J=',ja,'expected=',idimnh(ja),'from dens=',ndr2nh
            stop
        endif
        call  write_dens(iaaa,ndr2nh,r2nh,fnamern,name_new)
!write(836,*)iaaa,ja,ndr2nhh,idimnh(ja),name_new
!do ii=1,ndr2nh
!write(836,*)r2nh(ii)%ib,r2nh(ii)%is,r2nh(ii)%i2,r2nh(ii)%i1
!write(836,*)r2nh(ii)%rho
!end do
      endif
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!      open(63,file='2_phon_dens_calc.dat',status='old',form='formatted',access='append')
       write(63,*)iaaa
!      close(63)
 enddo   ! ia_cal

return
end subroutine roo

!************************************************************************
      subroutine ropar3(iaaa,iba,imin,imax,isimax,rom,irom,nphon,rxy,rh,iirg,ndrh,ja,jb,lev)

            use anglib   ! angular momentum staff
      
            implicit double precision(a-h,o-z)
      
!            include 'types_phon_dens.inc'
            include 'input_phon_dens.inc'
            include 'formats_phon_dens.inc'
      
      
            type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
            type (roc_typ), dimension(:), allocatable :: roc,rom
            type (ro_real_typ), dimension(:), allocatable :: rh
            type(rho_typ), dimension(:), allocatable :: rop
            type (level_typ), dimension(:),allocatable :: lev
      
            integer, dimension (:), allocatable :: ndx,ndc
      
            type (phon_typ), dimension (:,:), allocatable :: phonbs
            integer, dimension (:), allocatable :: nphon
            double precision, dimension (:,:,:), allocatable :: xx
            double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
            double precision, dimension (:,:,:,:,:), allocatable :: rho,rhon
            double precision, dimension (:,:,:), allocatable :: rr
      
      !      double precision, dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph
            real(kind=4), dimension(:,:,:,:,:), allocatable :: rnp,rpp,rnh,rph,rxy
      
      
            character(len=30) fnamer1,fnamer2,fnamer,fnameo
      
            character(len=10) fnamern,name_new
            character(len=2) denst
            integer*8 irom,ii
      
      
            character*6 nlam
      
      
!            iirg=0
!            ndrh=100000000
      
            if (irom.ne.0) then
      
            nff=2
            ndla=nphon(nff-1)
      
            allocate(rr(imin:imax,imin:imax,0:isimax))
            rr=0.d0
            
      
             do ii=1,irom

                  
      
               isir=rom(ii)%is
               iba=rom(ii)%ib
               ilam=rom(ii)%i1
               ilamp=rom(ii)%i2

      
                rr(imin:imax,imin:imax,isir)=rr(imin:imax,imin:imax,isir)+rxy(imin:imax,imin:imax,isir,ilam,ilamp)*rom(ii)%rho
      
              enddo ! over ii

      
!              iirg=0
!               do iba=1,nphon(nff)
                do isi=0,isimax
                 if(abs(jb-ja).le.isi.and.isi.le.(jb+ja)) then
                 do i_2=imin,imax
                  do i_1=imin,imax
                   if(abs(lev(i_1)%j-lev(i_2)%j).le.2*isi.and.2*isi.le.(lev(i_1)%j+lev(i_2)%j)) then
                !if (dabs(rr(i_1,i_2,isi)).gt.1.d-10) then
      
                iirg=iirg+1
                if (iirg.gt.ndrh) then
                 write(*,*)'Increase ndrh in ropar3'
                 stop
                endif

                rh(iirg)%ib=iba
                rh(iirg)%is=int(isi,1)
                rh(iirg)%i1=int(i_1,2)
                rh(iirg)%i2=int(i_2,2)
                rh(iirg)%rho=real(rr(i_1,i_2,isi))
      !          if (iba.eq.iaaa.and.i1.eq.i2.and.isir.eq.0) write(937,'(a5,5i5,f15.10)')denst,iaaa,iba,i_1,i_2,isir,rh(iirg)%rho
                endif
      
                   enddo ! i_1
                   enddo ! i_2
                   end if
                   enddo
              !end if
!               enddo  ! iba
            endif ! irom
!            deallocate(rr)
            return
            end subroutine ropar3
!*********************************************************
!*********************************************************
!*********************************************************
      subroutine write_dens(iaaa,iirg,rh,fnamern,name_new)

      implicit double precision (a-h,o-z)

      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'
      type (ro_real_typ), dimension(:), allocatable :: rh
      character(len=10) fnamern,name_new
      character(len=6) nlam

            if (iirg.gt.0) then
                       write(nlam,'(i6.6)')iaaa                  
                        open(73,file='scratch/'//name_new//'_'//nlam,status='unknown',form='unformatted')
                        write(73)iirg
                        !write(73)(rh(iii)%ib,iii=1,iirg)
                        !write(73)(rh(iii)%is,iii=1,iirg)
                        !write(73)(rh(iii)%i1,iii=1,iirg)
                        !write(73)(rh(iii)%i2,iii=1,iirg)
                        write(73)(rh(iii)%rho,iii=1,iirg)
                        close(73)
            endif
     end subroutine write_dens
!*************************************
!*************************************
!*************************************
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
      subroutine selphon(nf,phonbs,nphon,iphous,ns1,myid)

      implicit double precision (a-h,o-z)

!      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous,list_2ph


      allocate(iphous(nphon(nf)))
      iphous=0


!      xthr_min=0.0d0
!      xthr_max=10000000000000000.0d0

      allocate(list_2ph(nphon(nf)))
      list_2ph=0

      open(881,file='2_phon_dens_list.dat',status='old',form='formatted')
      do while(.not.eof(881))
         read(881,*)ii
         list_2ph(ii)=1
      enddo
      close(881)
      

      ii=0
      do i=1,nphon(nf)

!      xene=phonbs(nf,i)%enf
!       jjf=phonbs(nf,i)%jj
!        do j=1,nlist
!       if (xene.le.xthr_max.and.xene.gt.xthr_min) then
        if (list_2ph(i).eq.1) then     
          ii=ii+1
          iphous(ii)=i
        endif
!          enddo
      enddo

!      write(*,*)'Energy threshold for 2 phon. dens.',xthr_min,xthr_max
      if(myid.eq.0) write(*,*)' Number of selected phonons a)',ii

      ns1=ii

      end subroutine selphon

!**********************************************************************
      subroutine selphon2(nf,phonbs,nphon,iphous,ns2,myid)

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

      if (myid.eq.0) write(*,*)' Number of selected phonons b) ',ii
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
!************************************************************
     subroutine read_amp_ia(ia,fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)
      
      !      include 'types_phon_dens.inc'
      
            type (amp2_typ), dimension(:), allocatable :: xcc
            integer, dimension (:), allocatable :: ndcc
      
            character(len=30)fname
      
            open(2,file=fname,status='old',form='unformatted')
      
            ilamp=0
            ndimj=100000
            ii=0
      
            
      
            do while (.not.eof(2))
      
              read(2)ipar,ijj,no,idphon
      !        if (allocated(xcc)) deallocate(xcc)
      !        allocate(xcc(idphon))
      !        read(2)(xcc(i)%ig,xcc(i)%is,i=1,idphon)
      
              ii=ii+no
      
              if (ia.gt.ii) then 
                do ilam=1,no+1
                  read(2)
                enddo
      
              else
              if (allocated(xcc)) deallocate(xcc)
                allocate(xcc(idphon))
                read(2)(xcc(i)%ig,xcc(i)%is,i=1,idphon)
              do ilam=1,ia-ii+no-1
                 read(2)
              enddo
                
                read(2)(xcc(i)%am,i=1,idphon)
                close(2) 
               return 
              endif 
      
            enddo
      
      
            return
            end subroutine read_amp_ia

      !******************************************************************************

      subroutine read_ampl(fname,cc,n1,n2,n3,iphous2,ipozl)

            implicit double precision (a-h,o-z)
      
      !      include 'types_phon_dens.inc'
      
            type (amp2_typ), dimension(:), allocatable :: xcc
            double precision, dimension (:,:,:), allocatable :: cc
            integer, dimension (:), allocatable :: ndcc,iphous2,ipozl
            real :: xxr
            character(len=30)fname
      
            allocate (cc(n1,n2,n3))
            cc=0.0d0
            illl=1
      
            open(2,file=fname,status='old',form='unformatted')
      
            ilamp=0
      
            do while (.not.eof(2))
      
             read(2)ipar,ijj,no,idphon
             if (.not.allocated(xcc)) allocate(xcc(idphon))
             read(2)(xcc(i)%ig,xcc(i)%is,i=1,idphon)
              
            do ilam=1,no
      
             if (iphous2(ilam+ilamp).eq.1) then
      
             read(2)(xcc(i)%am,i=1,idphon)
      
             do i=1,idphon
              cc(illl,xcc(i)%is,xcc(i)%ig)=xcc(i)%am
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
      
            deallocate(xcc)
      
            enddo
      
            close(2)
      
            
            return
            end subroutine read_ampl
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
      
