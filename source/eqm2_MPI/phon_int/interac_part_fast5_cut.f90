!*     Phoninterac contains routines for calculation of redefined phonon
!*     interaction in p-h J-scheme
!*     last update 10.2.2015

      module phoninterac

      use eofmod

      contains
      
      subroutine vinth(ipart,ifmx,isimax,lev,xthrun_min,xthrun_max,myid,numprocs)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_int.inc'
      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'


!      double precision, dimension(:,:,:,:,:),allocatable ::fp,fpn
      real(kind=4), dimension(:,:,:,:,:),allocatable ::fp,fpn

      double precision, dimension(:,:,:),allocatable :: vint

      integer, dimension(:), allocatable :: jphon,ironp,iropp,ironh,iroph,ndcamn,ndcamp,ius,ig_resh

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh

      integer ierr,myid,numprocs,irc
      integer :: n_seg


      character*10 namer,namerpp,namernp,namerph,namernh,namev
      character*30 namefp,namefpn,namecp,namecn
      character*4 nlam

   

      
      ndrho=10000000
      allocate(rh(ndrho))

      allocate(ius(ifmx))
      ius=0

     
      jmin=0
      jmax=isimax

      ndlam=10000      
      ndamp=10000
      ifmxx=10000
      allocate (jphon(ifmxx))
      jphon=0


      ius=0
      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
       if (en.le.xthrun_max.and.en.ge.xthrun_min) ius(i)=1
      enddo

      close(3)

      ifmx=i

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo

      write(*,*)' Number of 0hom+1hom+2hom+3hom phonons ',icount


      

      if (ipart.eq.1) then 
              namecp='1phonon/1f_cn.dat'
              namecn='1phonon/1f_cp.dat'
              namefp='fmat_n.dat'
              namefpn='fmat_pn.dat'
              namerpp='1f_rnp.dat'
              namernp='1f_rpp.dat'
              namerph='1f_rnh.dat'
              namernh='1f_rph.dat'
              namev='V_phon_h_n.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
              ihminp=ihpmn
              ipmaxp=ippmx
      endif

      if (ipart.eq.-1) then 
              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
              namefp='fmat_p.dat'
              namefpn='fmat_np.dat'
              namerpp='1f_rpp.dat'
              namernp='1f_rnp.dat'
              namerph='1f_rph.dat'
              namernh='1f_rnh.dat'
              namev='V_phon_h_p.dat'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
              ihminp=ihnmn
              ipmaxp=ipnmx

      endif

!     read(*,*)ig_min,ig_max
     ig_min=1
     ig_max=icount

  if (myid.eq.0) then
      write(*,*)'Phonons for which interaction is calculated?        <',1,ifmx,'       > '
      write(*,*)ig_min,ig_max
  endif



!c     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihmin,ipmax,ihmin,ihmax,fp)
!c     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihminp,ipmaxp,ihmin,ipmax,fpn)
!c     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
!c     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      

     if (myid.eq.0) then
      if (ipart.eq.1) write(*,*)'Calculation of phonon neutron hole interaction'
      if (ipart.eq.-1) write(*,*)'Calculation of phonon proton hole interaction'
     endif 



      allocate(vint(0:isimax,ihmin:ihmax,ihmin:ihmax))

      vint=0.0d0

!     reshufling
      allocate(ig_resh(ifmx))
      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) then
                 icount=icount+1
                 ig_resh(icount)=i
       endif
      enddo

      call rperm(icount, ig_resh)

!      do ig=1,ifmx
!       do ig=ig_min,ig_max

      if (mod(icount,numprocs).eq.0) then
              n_seg=icount/numprocs
      else
              n_seg=icount/numprocs+1
      endif
!      n_rem=mod(icont,numprocs)

if (myid.eq.0) write(*,*) ' size of segment = ',n_seg

!      do ig=1,ifmx
!      do ig=ig_min,ig_max
  do igg=myid*n_seg+1,min((myid+1)*n_seg,icount)

       ig=ig_resh(igg)
       write(*,*) ' Process #',myid,'  calculating ig = ',ig


      iirg=0

!      write(5)ig

      if (ius(ig).ne.0) then


       write(nlam,'(i4.4)')ig

       open(5,file='scratch/'//namev//'_'//nlam,status='unknown',form='unformatted')



!       write(997,*)ig,ifmx


       jig=jphon(ig)


!      call readro(33,ig,ronp,nronp)
!      call readro(34,ig,ropp,nropp)
!      call readro(43,ig,ronh,nronh)
!      call readro(44,ig,roph,nroph)

!      call readro('1f_rnp.dat',ig,ronp,nronp)
!      call readro('1f_rpp.dat',ig,ropp,nropp)
!      call readro('1f_rnh.dat',ig,ronh,nronh)
!
!      call readro('1f_rph.dat',ig,roph,nroph)

       call readro(namernp,ig,ronp,nronp)
       call readro(namerpp,ig,ropp,nropp)
       call readro(namernh,ig,ronh,nronh)
       call readro(namerph,ig,roph,nroph)


    
      do ib=1,ifmx
       jib=jphon(ib)

!        isi_min=iabs(jib-jig)
!        isi_max=jib+jig
        isi_min=0
        isi_max=isimax

      
!      do isi=0,isimax 
!      ifaz=(-1)**(jib+isi)
       
      allocate(ironp(nronp))
      ironp=0 
      allocate(iropp(nropp))
      iropp=0 
      allocate(ironh(nronh))
      ironh=0 
      allocate(iroph(nroph))
      iroph=0 

      call rosub2(ib,ropp,nropp,iropp,nropps)
      call rosub2(ib,ronp,nronp,ironp,nronps)     
      call rosub2(ib,roph,nroph,iroph,nrophs)
      call rosub2(ib,ronh,nronh,ironh,nronhs)     
    
     
!      do ih1=ihmin,ihmax
!       jh1=lev(ih1)%j 
!      do ih2=ihmin,ihmax
!       jh2=lev(ih2)%j
             
        vint=0.d0

        do ii=1,nropps
         i1=ropp(iropp(ii))%i1
         i2=ropp(iropp(ii))%i2
         isi=ropp(iropp(ii))%j
         ifaz=(-1)**(jib+isi)
        do ih1=ihmin,ihmax
          do ih2=ihmin,ihmax
         vint(isi,ih1,ih2)=vint(isi,ih1,ih2)+dfloat(ifaz)*ropp(iropp(ii))%ro*fp(isi,i1,i2,ih1,ih2)
          enddo
        enddo
        enddo

        do ii=1,nronps
         i1=ronp(ironp(ii))%i1
         i2=ronp(ironp(ii))%i2
         isi=ronp(ironp(ii))%j
         ifaz=(-1)**(jib+isi)

        do ih1=ihmin,ihmax
          do ih2=ihmin,ihmax
         vint(isi,ih1,ih2)=vint(isi,ih1,ih2)+dfloat(ifaz)*ronp(ironp(ii))%ro*fpn(isi,i1,i2,ih1,ih2)
          enddo
        enddo
 
        enddo

        do ii=1,nrophs
         i1=roph(iroph(ii))%i1
         i2=roph(iroph(ii))%i2
         isi=roph(iroph(ii))%j
         ifaz=(-1)**(jib+isi)
        do ih1=ihmin,ihmax
          do ih2=ihmin,ihmax
         vint(isi,ih1,ih2)=vint(isi,ih1,ih2)+0.5d0*dfloat(ifaz)*roph(iroph(ii))%ro*fp(isi,i1,i2,ih1,ih2)
          enddo
        enddo
 
        enddo

        do ii=1,nronhs
         i1=ronh(ironh(ii))%i1
         i2=ronh(ironh(ii))%i2
         isi=ronh(ironh(ii))%j
         ifaz=(-1)**(jib+isi)
        do ih1=ihmin,ihmax
          do ih2=ihmin,ihmax
         vint(isi,ih1,ih2)=vint(isi,ih1,ih2)+dfloat(ifaz)*ronh(ironh(ii))%ro*fpn(isi,i1,i2,ih1,ih2)
          enddo
        enddo
        enddo


!        factor=dfloat((-1)**(jib+(jh1+jh2)/2))*
!     *((2*jib+1)*(2*jig+1)*(2*isi+1))**0.5d0

        do ii=1,ndcamn(ib)
          ipp=camn(ib,ii)%par
          ihp=camn(ib,ii)%hol
          campp=camn(ib,ii)%am
         do jj=1,ndcamp(ig)
           ip=camp(ig,jj)%par
           ih=camp(ig,jj)%hol
!          if (ih.eq.ih2) then 
        
           jp=lev(ip)%j
           jh2=lev(ih)%j
           campn=camp(ig,jj)%am

  
          do ih1=ihmin,ihmax
          jh1=lev(ih1)%j

          do isi=isi_min,isi_max
 

           factor=dfloat((-1)**(jib+(jh1+jh2)/2))*((2*jib+1)*(2*jig+1)*(2*isi+1))**0.5d0

           xsixj=sixj(2*jib,2*isi,2*jig,jh2,jp,jh1)          
           vint(isi,ih1,ih)=vint(isi,ih1,ih)+factor*campp*campn*xsixj*fpn(jib,ipp,ihp,ih1,ip)
           enddo
           enddo

!          endif
         enddo
        enddo
             
!c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.1) 
!c     *write(997,102)isi,ig,ib,ih1,ih2,vint

!c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.-1) 
!c     *write(998,102)isi,ig,ib,ih1,ih2,vint


!        if (dabs(vint).gt.xrotrunc) 
!     *write(5)ib,isi,ih1,ih2,vint

!       do isi=0,isimax
       do isi=isi_min,isi_max

       do ih1=ihmin,ihmax
       do ih2=ihmin,ihmax

       if (dabs(vint(isi,ih1,ih2)).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in rop!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ih1
               rh(iirg)%i2=ih2
               rh(iirg)%rho=real(vint(isi,ih1,ih2))
!c        write(33)ib,isi,ip1,ip2,ronp
       endif
       enddo
       enddo
       enddo




!        enddo ! loop ih2
!      enddo ! loop ih1
      deallocate (ironp,ironh,iropp,iroph)
!      enddo ! loop isi
      enddo ! loop ib

!      endif

      write(5)ig,iirg
      if (iirg.gt.0) then
      write(5)(rh(iii)%ib,iii=1,iirg)
      write(5)(rh(iii)%isi,iii=1,iirg)
      write(5)(rh(iii)%i1,iii=1,iirg)
      write(5)(rh(iii)%i2,iii=1,iirg)
      write(5)(rh(iii)%rho,iii=1,iirg)
      endif

!      do iii=1,iirg
!       write(332,'(5i5,f15.10)')ig,rh(iii)%ib,rh(iii)%isi,rh(iii)%i1,rh(iii)%i2,rh(iii)%rho
!      enddo
!      write(332,*)ig,iirg


      close(5)

      endif

!      if (iirg.eq.0) then
!      write(5)iirg
!      write(5)iirg
!      write(5)iirg
 !     write(5)iirg
 !     write(5)dfloat(iirg)
 !     endif



!      write(5)0,0,0,0,0.d0      
      enddo ! loop ig      

      deallocate(rh,camp,camn,jphon,fp,fpn,ronp,ronh,ropp,roph)

!      write(5)10000000
!      write(5)0,0,0,0,0.d0 

!      close(33)
!      close(34)
!      close(43)
!      close(44)
 !     close(5)
      return
      end subroutine vinth
!*************************************************************************
      subroutine vintp(ipart,ifmx,isimax,lev,xthrun_min,xthrun_max,myid,numprocs)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_int.inc'
      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'


!      double precision, dimension(:,:,:,:,:),allocatable ::fp,fpn
      real(kind=4), dimension(:,:,:,:,:),allocatable ::fp,fpn

      double precision, dimension(:,:,:),allocatable :: vint

      integer, dimension(:), allocatable :: jphon,ironp,iropp,ironh,iroph,ndcamn,ndcamp,ius,ig_resh

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      type(level_typ),dimension(*) :: lev
      type(ro_typ),dimension(:), allocatable :: rh

      integer ierr,myid,numprocs,irc
      integer :: n_seg



      character*10 namer,namerpp,namernp,namerph,namernh
      character*30 namefp,namefpn,namecp,namecn
      character*10 namev
      character*4 nlam

      ndrho=100000000
      allocate(rh(ndrho))


      allocate(ius(ifmx))
      ius=0
      
     
      jmin=0
      jmax=isimax

      ndlam=10000      
      ndamp=10000
      ifmxx=10000
      allocate (jphon(ifmxx))
      jphon=0

      ius=0

      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
!c       write(991,*)i,ipar,ijj,en
       jphon(i)=ijj
!       if (en.le.xthrun) ius(i)=1
       if (en.le.xthrun_max.and.en.ge.xthrun_min) ius(i)=1

      enddo

      close(3)

      ifmx=i

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo

      if (myid.eq.0) write(*,*)' Number of 0hom+1hom+2hom+3hom phonons ',icount


      if (ipart.eq.1) then 
              namecp='1phonon/1f_cn.dat'
              namecn='1phonon/1f_cp.dat'
              namefp='fmat_n.dat'
              namefpn='fmat_pn.dat'
              namerpp='1f_rnp.dat'
              namernp='1f_rpp.dat'
              namerph='1f_rnh.dat'
              namernh='1f_rph.dat'
              namev='V_phon_p_n'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
              ihminp=ihpmn
              ipmaxp=ippmx
      endif

      if (ipart.eq.-1) then 
              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
              namefp='fmat_p.dat'
              namefpn='fmat_np.dat'
              namerpp='1f_rpp.dat'
              namernp='1f_rnp.dat'
              namerph='1f_rph.dat'
              namernh='1f_rnh.dat'
              namev='V_phon_p_p'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
              ihminp=ihnmn
              ipmaxp=ipnmx

      endif


!     read(*,*)ig_min,ig_max
     ig_min=1
     ig_max=icount

  if (myid.eq.0) then 
      write(*,*)'Phonons for which interaction is calculated?        <',1,ifmx,'       > '
      write(*,*)ig_min,ig_max
  endif 


!     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihmin,ipmax,ipmin,ipmax,fp)
!     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihminp,ipmaxp,ihmin,ipmax,fpn)
!     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
!     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      
      if (myid.eq.0) then 
       if (ipart.eq.1) write(*,*)'Calculation of phonon neutron particle interaction'
       if (ipart.eq.-1) write(*,*)'Calculation of phonon proton particle interaction'
      end if

      allocate(vint(0:isimax,ipmin:ipmax,ipmin:ipmax))
      vint=0.0d0

!     reshufling 
      allocate(ig_resh(ifmx))
      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) then
                 icount=icount+1
                 ig_resh(icount)=i
       endif
      enddo

      call rperm(icount, ig_resh)


!      do ig=1,ifmx
!       do ig=ig_min,ig_max 

      if (mod(icount,numprocs).eq.0) then
              n_seg=icount/numprocs
      else
              n_seg=icount/numprocs+1
      endif
!      n_rem=mod(icont,numprocs)

if (myid.eq.0) write(*,*) ' size of segment = ',n_seg

       
      do igg=myid*n_seg+1,min((myid+1)*n_seg,icount)

       ig=ig_resh(igg)
       write(*,*) ' Process #',myid,'  calculating ig = ',ig

   
       iirg=0 

       jig=jphon(ig)

!      write(5)ig
      if (ius(ig).ne.0) then

!       open(5,file=namev,status='unknown'
!     *,form='unformatted')



       write(nlam,'(i4.4)')ig
       
       open(5,file='scratch/'//namev//'_'//nlam,status='unknown',form='unformatted')


!       write(996,*)ig,ifmx

      call readro(namernp,ig,ronp,nronp)
      call readro(namerpp,ig,ropp,nropp)
      call readro(namernh,ig,ronh,nronh)
      call readro(namerph,ig,roph,nroph)


    
      do ib=1,ifmx
       jib=jphon(ib)
       isi_min=0 !iabs(jib-jig)
       isi_max=isimax !jib+jig

      
!      do isi=0,isimax 
      ifaz=(-1)**(jib+isi)
       
      allocate(ironp(nronp))
      ironp=0 
      allocate(iropp(nropp))
      iropp=0 
      allocate(ironh(nronh))
      ironh=0 
      allocate(iroph(nroph))
      iroph=0 

      call rosub2(ib,ropp,nropp,iropp,nropps)
      call rosub2(ib,ronp,nronp,ironp,nronps)     
      call rosub2(ib,roph,nroph,iroph,nrophs)
      call rosub2(ib,ronh,nronh,ironh,nronhs)     
    
     
!      do ip1=ipmin,ipmax
!       jp1=lev(ip1)%j 
!      do ip2=ipmin,ipmax
!       jp2=lev(ip2)%j
             
        vint=0.d0

        do ii=1,nropps
         i1=ropp(iropp(ii))%i1
         i2=ropp(iropp(ii))%i2
         isi=ropp(iropp(ii))%j
         ifaz=(-1)**(jib+isi)

        do ip1=ipmin,ipmax
         do ip2=ipmin,ipmax
         vint(isi,ip1,ip2)=vint(isi,ip1,ip2)+0.5d0*dfloat(ifaz)*ropp(iropp(ii))%ro*fp(isi,i1,i2,ip1,ip2)
         enddo
        enddo

        enddo

        do ii=1,nronps
         i1=ronp(ironp(ii))%i1
         i2=ronp(ironp(ii))%i2
         isi=ronp(ironp(ii))%j
         ifaz=(-1)**(jib+isi)

        do ip1=ipmin,ipmax
         do ip2=ipmin,ipmax
         vint(isi,ip1,ip2)=vint(isi,ip1,ip2)+dfloat(ifaz)*ronp(ironp(ii))%ro*fpn(isi,i1,i2,ip1,ip2)
         enddo
        enddo

        enddo

        do ii=1,nrophs
         i1=roph(iroph(ii))%i1
         i2=roph(iroph(ii))%i2
         isi=roph(iroph(ii))%j
         ifaz=(-1)**(jib+isi)

        do ip1=ipmin,ipmax
         do ip2=ipmin,ipmax
          vint(isi,ip1,ip2)=vint(isi,ip1,ip2)+dfloat(ifaz)*roph(iroph(ii))%ro*fp(isi,i1,i2,ip1,ip2)
         enddo
        enddo
        enddo

        do ii=1,nronhs
         i1=ronh(ironh(ii))%i1
         i2=ronh(ironh(ii))%i2
         isi=ronh(ironh(ii))%j
         ifaz=(-1)**(jib+isi)
        do ip1=ipmin,ipmax
         do ip2=ipmin,ipmax
         vint(isi,ip1,ip2)=vint(isi,ip1,ip2)+dfloat(ifaz)*ronh(ironh(ii))%ro*fpn(isi,i1,i2,ip1,ip2)
         enddo
        enddo
        enddo


        do ii=1,ndcamn(ib)
          ipp=camn(ib,ii)%par
          ihp=camn(ib,ii)%hol
          campp=camn(ib,ii)%am
         do jj=1,ndcamp(ig)
           ip=camp(ig,jj)%par
           ih=camp(ig,jj)%hol 

!          if (ip.eq.ip1) then 
           jh=lev(ih)%j
           jp1=lev(ip)%j

           campn=camp(ig,jj)%am

!        do isi=isi_min,isi_max

!        do isi=0,isimax

       do ip2=ipmin,ipmax
        jp2=lev(ip2)%j

        do isi=isi_min,isi_max

           factor=dfloat((-1)**(isi-jig+(jp1+jp2)/2))*((2*jib+1)*(2*jig+1)*(2*isi+1))**0.5d0

           xsixj=sixj(2*jib,2*isi,2*jig,jp1,jh,jp2)          
           vint(isi,ip,ip2)=vint(isi,ip,ip2)-factor*campp*campn*xsixj*fpn(jib,ipp,ihp,ih,ip2)

       enddo

        
        enddo

!          endif
         enddo
        enddo



       do isi=isi_min,isi_max

        do ip1=ipmin,ipmax
         do ip2=ipmin,ipmax


       if (dabs(vint(isi,ip1,ip2)).gt.xrotrunc) then
               iirg=iirg+1
               if (iirg.gt.ndrho) then
                  write(*,*)' Increase dimension of ndrho in rop!!'
                  stop
               endif
               rh(iirg)%ib=ib
               rh(iirg)%isi=isi
               rh(iirg)%i1=ip1
               rh(iirg)%i2=ip2
               rh(iirg)%rho=real(vint(isi,ip1,ip2))
!c        write(33)ib,isi,ip1,ip2,ronp
       endif

         enddo
        enddo

       enddo




!        enddo ! loop ih2
!      enddo ! loop ih1
      deallocate (ironp,ironh,iropp,iroph)
!      enddo ! loop isi
      enddo ! loop ib
  
!      endif

      write(5)ig,iirg
      if (iirg.gt.0) then
      write(5)(rh(iii)%ib,iii=1,iirg)
      write(5)(rh(iii)%isi,iii=1,iirg)
      write(5)(rh(iii)%i1,iii=1,iirg)
      write(5)(rh(iii)%i2,iii=1,iirg)
      write(5)(rh(iii)%rho,iii=1,iirg)
      endif

!      do iii=1,iirg
!       write(331,'(5i5,f15.10)')ig,rh(iii)%ib,rh(iii)%isi,rh(iii)%i1,rh(iii)%i2,rh(iii)%rho
!      enddo
!      write(331,*)ig,iirg


      close(5)

      endif

!      if (iirg.eq.0) then
!      write(5)iirg
!      write(5)iirg
!      write(5)iirg
!      write(5)iirg
!      write(5)dfloat(iirg)
!      endif



!      write(5)0,0,0,0,0.d0 

      enddo ! loop ig      

      deallocate(rh,camp,camn,jphon,fp,fpn,ronp,ronh,ropp,roph)

!      write(5)10000000
!      write(5)0,0,0,0,0.d0 

!      close(33)
!      close(34)
!      close(43)
!      close(44)
!      close(5)
      return
      end subroutine vintp


!************************************************************************
      subroutine readcam(fname,ndimi,ndimj,cam,ndcc)

      implicit double precision (a-h,o-z)

      include 'formats_phon_int.inc'
      include 'types_phon_int.inc'

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

!
!*****************************************************************************
      subroutine readro(fname,ig,ron,ndgg)

      implicit double precision (a-h,o-z)

!      include 'formats_eqm.inc'
      include 'types_phon_int.inc'

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=10)fname
      character(len=4)nlam
      logical je_tam

      ifile=33

      write(nlam,'(i4.4)')ig

      inquire(file='scratch/'//fname//'_'//nlam,exist=je_tam)

      if (je_tam.eqv..FALSE.) then
        ndgg=0
        return
      endif


      open(ifile,file='scratch/'//fname//'_'//nlam,status='unknown',form='unformatted')

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




      subroutine rosub(j,ilam,rop,nrop,irop,nrops)

      implicit double precision (a-h,o-z)

      include 'types_phon_int.inc'

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
      subroutine rosub2(ilam,rop,nrop,irop,nrops)

      implicit double precision (a-h,o-z)

      include 'types_phon_int.inc'

      type(rho_typ), dimension(:), allocatable :: rop
      integer, dimension(:), allocatable :: irop

      ii=0
      do i=1,nrop
      if (rop(i)%ilap.eq.ilam) then
              ii=ii+1
              irop(ii)=i
             endif

      enddo

      nrops=ii

      end subroutine rosub2
!***********************************************************************


      subroutine readfin(fname,jmin,jmax,imin,imax,kmin,kmax,fpp)

      implicit double precision (a-h,o-z)

      include 'formats_phon_int.inc'

!      double precision, dimension(:,:,:,:,:), allocatable ::fpp
      real(kind=4), dimension(:,:,:,:,:), allocatable ::fpp
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

!      read(2)itt,ipt,ijt,i,j,k,l,vint


!c        if (ipt.ne.ipar) goto 11
      if (ijt.gt.jmax.or.ijt.lt.jmin) goto 11
      
      if (i.gt.imax.or.i.lt.imin.or.j.gt.imax.or.j.lt.imin.or.k.gt.kmax.or.k.lt.kmin.or.l.gt.kmax.or.l.lt.kmin) goto 11


      fpp(ijt,i,j,k,l)=real(vint)
  
 11   enddo
      close(2)


      return
      end subroutine readfin
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


     
end module phoninterac 
      
      
