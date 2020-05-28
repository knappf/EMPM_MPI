!*     Phoninterac contains routines for calculation of redefined phonon
!*     interaction in p-h J-scheme

!*     last update 6.2.2020

!      modification used natural orbitals calculations  

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

subroutine vintn01(nf,ipcal,jcal,levn,levp)

      use anglib   ! angular momentum staff
      use input_sp
     
      implicit double precision(a-h,o-z)
     
      include 'types_ndgi_int.inc'
      include 'input_ndgi_int.inc'
      include 'formats_ndgi_int.inc'
     
     
      double precision, dimension(:,:,:,:,:),allocatable ::fp,fn,fpn
     
      integer, dimension(:), allocatable :: jphon,jphonm,ironp,iropp,ironh,iroph,ndcamn,ndcamp,iphon,iphonm,phonus
     
     double precision, allocatable, dimension(:,:) :: kin_p, kin_n
     
      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      type(level_typ),dimension(:), allocatable :: levn,levp
     
      character*10 namer
      character*30 namefp,namefn,namefpn,namecp,namecn,namerpp,namernp,namerph,namernh,namev,nameff
          
     
                 
     ndamp=1600
     ifmx=1600
     ifmmx=24000
     allocate (jphon(ifmx))  ! 1phonon 
     jphon=0
     allocate (iphon(ifmx))  ! 1phonon
     iphon=0
          
     
     open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')
     
     do while (.not.eof(3))
      read(3)i,ipar,ijj,en
      jphon(i)=ijj
      iphon(i)=ipar
     enddo
     
     close(3)
     
     ifmx=i
     ndlam=ifmx
        
     allocate(phonus(ifmx))
     phonus=1
     namecp='1phonon/1f_cp.dat'
     namecn='1phonon/1f_cn.dat'
     namefp='fmat_p.dat'
     namefn='fmat_n.dat'
     namefpn='fmat_pn.dat'
     namev='Vint_phon01.dat'  
       
       
           write(*,*)'jmin,jmax',jmin,jmax


           call input_kin(kin_p,kin_n)
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
     
     
           open(5,file=namev,status='unknown',form='unformatted') !,access='append')
     
     

           
 
     do isi=1,ifmx   ! cycle over 1-phonon states
     !
            ips=iphon(isi)  ! parity               
            jisi=jphon(isi)   ! j 

            vint=0.d0
            
           if (jisi.eq.0.and.ips.eq.1) then 

           do jj=1,ndcamp(isi)
              ip=camp(isi,jj)%par
              ih=camp(isi,jj)%hol
              jp=levp(ip)%j
              jh=levp(ih)%j
              if (jp.eq.jh) then 
                campn=camp(isi,jj)%am
                vint=vint+(dfloat(jh+1))**0.5d0*campn*kin_p(ih,ip)
              endif

              do i1=ihpmn,ihpmx
                 ji1=levp(i1)%j
                campn=camp(isi,jj)%am
                vint=vint+0.5d0*(dfloat(ji1+1))**0.5d0*campn*fp(0,ih,ip,i1,i1)
              enddo
              do i1=ihnmn,ihnmx
                 ji1=levn(i1)%j
                campn=camp(isi,jj)%am
                vint=vint+campn*(dfloat(ji1+1))**0.5d0*fpn(0,ih,ip,i1,i1)
               enddo
           enddo    

           do jj=1,ndcamn(isi)
                 ip=camn(isi,jj)%par
                 ih=camn(isi,jj)%hol
                 jp=levn(ip)%j
                 jh=levn(ih)%j
                 if (jp.eq.jh) then 
                   campn=camn(isi,jj)%am
                   vint=vint+(dfloat(jh+1))**0.5d0*campn*kin_n(ih,ip)
                 endif

                 do i1=ihnmn,ihnmx
                      ji1=levn(i1)%j
                       campn=camn(isi,jj)%am
                       vint=vint+0.5d0*(dfloat(ji1+1))**0.5d0*campn*fn(0,ih,ip,i1,i1)
                 enddo

                 do i1=ihpmn,ihpmx
                       ji1=levp(i1)%j
                       campn=camn(isi,jj)%am
                     vint=vint+campn*(dfloat(ji1+1))**0.5d0*fpn(0,i1,i1,ih,ip)
                 enddo
           enddo    
                                      

             if (dabs(vint).gt.xrotrunc) write(5)isi,vint
     
             if (dabs(vint).gt.xrotrunc) write(997,*)isi,vint
     
          
           endif 
     enddo ! loop isi
     
                
     deallocate(camp,camn,fp,fpn,kin_n,kin_p)


           close(5)
           return
           end subroutine vintn01
     
!***************************************************************************
  
 !     
subroutine vintn12(nf,ipcal,jcal,levn,levp)

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
                namev='Vint_phon12.dat'
!              endif

!              if (nf.eq.3) then 
!                namerpp='2f_rpp.dat'
!                namernp='2f_rnp.dat'
!                namerph='2f_rph.dat'
!                namernh='2f_rnh.dat'
!                namev='Vint_phon23.dat'
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
      
      write(*,*)'Calculation of redefined interaction'

      open(5,file=namev,status='unknown',form='unformatted') !,access='append')


      do ia=1,ifmmx

      write(*,*)ia,ifmmx

       jia=jphonm(ia)
       ipa=iphonm(ia)

       xfact=(dfloat(2*jia+1))**(-1.d0)



       if (jia.eq.jcal.and.ipa.eq.ipcal) then


!      if (jia.eq.jcal.and.ipa.eq.ipcal) then 
!        if (jia.eq.jia) then 

!c       xfact=(dfloat(2*jia+1))**(-1.d0)

!c      write(5)ia
           call readro(namernp,ia,ronp,nronp)
           call readro(namerpp,ia,ropp,nropp)
           call readro(namernh,ia,ronh,nronh)
           call readro(namerph,ia,roph,nroph)

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



!      if (jisi.ne.jisiprev) then 

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

        if (dabs(vint).gt.xrotrunc) write(5)ia,ig,isi,xfact*dfloat(ifaz)*vint

!        if (dabs(vint).gt.xrotrunc) write(991,*)ia,ig,isi,xfact*dfloat(ifaz)*vint


!        enddo ! loop ih2
!      enddo ! loop ih1


      enddo ! loop isi


      enddo ! loop ib


      endif
!c      write(5)0,0,0,0,0.d0      
      enddo ! loop ig      

      
      deallocate(camp,camn,jphon,jphonm,fp,fpn,ronp,ronh,ropp,roph)
      deallocate (ironp,iropp,ironh,iroph)

!      write(5)10000000
!c      write(5)0,0,0,0,0.d0 

!      close(33)
!      close(34)
!      close(43)
!      close(44)
      close(5)
      return
      end subroutine vintn12



!***************************************************************************
      subroutine vintn23(nf,ipcal,jcal,levn,levp)

            use anglib   ! angular momentum staff
            
            implicit double precision(a-h,o-z)
            
            include 'types_ndgi_int.inc'
            include 'input_ndgi_int.inc'
            include 'formats_ndgi_int.inc'
            
            
            double precision, dimension(:,:,:,:,:),allocatable ::fp,fn,fpn
            double precision, dimension(:),allocatable ::vint
            integer, dimension(:), allocatable :: jphon,jphonm,ironp,iropp,ironh,iroph,ndcamn,ndcamp,iphon,iphonm,phonus
            type(amp_typ), dimension(:,:), allocatable :: camp,camn
            type(rho2_typ), dimension(:,:), allocatable :: c_pp,c_ph,c_np,c_nh
            Real(Kind=4), dimension(:),allocatable:: ronp,ropp,ronh,roph
            integer, allocatable, dimension(:):: idimnp,idimpp,idimph,idimnh
            type(level_typ),dimension(:), allocatable :: levn,levp
            character*10 namer
            character*30 namefp,namefn,namefpn,namecp,namecn,namerpp,namernp,namerph,namernh,namev,nameff
                 
                        
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
            
                  jmax2=0
                  open (3,file=nameff,status='old',form='unformatted')
            
                  do while (.not.eof(3))
                   read(3)i,ipar,ijj,en
                   jphonm(i)=ijj
                   if(ijj.gt.jmax2) jmax2=ijj
                   iphonm(i)=ipar
                  enddo
            
                  close(3)
            
                  ifmmx=i
                  idim2=i
               
                  allocate(phonus(ifmx))
                  phonus=1
            
            
                          namecp='1phonon/1f_cp.dat'
                          namecn='1phonon/1f_cn.dat'
                          namefp='fmat_p.dat'
                          namefn='fmat_n.dat'
                          namefpn='fmat_pn.dat'
            
                          if (nf.eq.2) then 
                            namerpp='1f_rpp.dat'
                            namernp='1f_rnp.dat'
                            namerph='1f_rph.dat'
                            namernh='1f_rnh.dat'
                            namev='Vint_phon12.dat'
                          endif
            
                          if (nf.eq.3) then 
                            namerpp='2f_rpp.dat'
                            namernp='2f_rnp.dat'
                            namerph='2f_rph.dat'
                            namernh='2f_rnh.dat'
                            namev='Vint_phon23.dat'
                          endif
            
            
                   write(*,*)'jmin,jmax',jmin,jmax


!**************
allocate(idimnp(0:jmax2),idimnh(0:jmax2),idimpp(0:jmax2),idimph(0:jmax2))
    idimnp=0
    idimnh=0
    idimph=0
    idimpp=0
    ndro=1500000

      open(1,file='input_tda_coup.dat',status='old',form='formatted')

      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      close(1)



      jmaxn=0
         do i=1,ippmx
          if (levp(i)%j.gt.jmaxn) jmaxn=levp(i)%j
         enddo

        do i=1,ipnmx
         if (levn(i)%j.gt.jmaxn) jmaxn=levn(i)%j
         enddo

      write(*,*)'jmax =',jmaxn
isimax=jmaxn
allocate(c_pp(0:jmax2,ndro),c_nh(0:jmax2,ndro),c_np(0:jmax2,ndro),c_ph(0:jmax2,ndro))

   do jj=0,jmax2

      do ib=1,idim2
      
      jb=jphonm(ib)!%j
       do isi=0,isimax
        if(abs(jb-jj).le.isi.and.isi.le.(jb+jj))then
         !write(19,*) !iaaa,ib,isi,ja,jb!ib,iaaa,isi,ja,jb
         do i_2=ipnmn,ipnmx
          do i_1=ipnmn,ipnmx
           if(abs(levn(i_1)%j-levn(i_2)%j).le.2*isi.and.2*isi.le.(levn(i_1)%j+levn(i_2)%j)) then
           idimnp(jj)=idimnp(jj)+1
           if(idimnp(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
           if(idimnp(jj).gt.ndro) stop
           c_np(jj,idimnp(jj))%ilap=ib
           c_np(jj,idimnp(jj))%j=isi
           c_np(jj,idimnp(jj))%i1=i_1
           c_np(jj,idimnp(jj))%i2=i_2
           !write(58)ib,isi,i_2,i_1
           !write(588,*) idimnp(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ihnmn,ihnmx
          do i_1=ihnmn,ihnmx
           if(abs(levn(i_1)%j-levn(i_2)%j).le.2*isi.and.2*isi.le.(levn(i_1)%j+levn(i_2)%j)) then         
           idimnh(jj)=idimnh(jj)+1
           if(idimnh(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
           if(idimnh(jj).gt.ndro) stop
           c_nh(jj,idimnh(jj))%ilap=ib
           c_nh(jj,idimnh(jj))%j=isi
           c_nh(jj,idimnh(jj))%i1=i_1
           c_nh(jj,idimnh(jj))%i2=i_2
!           write(59)ib,isi,i_2,i_1
!           write(590,*) idimnh(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ippmn,ippmx
          do i_1=ippmn,ippmx
           if(abs(levp(i_1)%j-levp(i_2)%j).le.2*isi.and.2*isi.le.(levp(i_1)%j+levp(i_2)%j)) then
           idimpp(jj)=idimpp(jj)+1
          if(idimpp(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
          if(idimpp(jj).gt.ndro) stop
          c_pp(jj,idimpp(jj))%ilap=ib
          c_pp(jj,idimpp(jj))%j=isi
          c_pp(jj,idimpp(jj))%i1=i_1
          c_pp(jj,idimpp(jj))%i2=i_2
!           write(60)ib,isi,i_2,i_1
!           write(600,*) idimpp(jj),ib,isi,i_2,i_1
           end if
          end do
         end do

         do i_2=ihpmn,ihpmx
          do i_1=ihpmn,ihpmx
           if(abs(levp(i_1)%j-levp(i_2)%j).le.2*isi.and.2*isi.le.(levp(i_1)%j+levp(i_2)%j)) then         
           idimph(jj)=idimph(jj)+1
           if(idimph(jj).gt.ndro) write(*,*) 'in eqm.f90 increase ndro'
           if(idimph(jj).gt.ndro) stop
           c_ph(jj,idimph(jj))%ilap=ib
           c_ph(jj,idimph(jj))%j=isi
           c_ph(jj,idimph(jj))%i1=i_1
           c_ph(jj,idimph(jj))%i2=i_2
           !write(61)ib,isi,i_2,i_1
           !write(610,*) idimph(jj),ib,isi,i_2,i_1
           end if
          end do
         end do
         
        end if!jj
       end do!isi
      end do !ib
     end do
!****************


            
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

                  ndro=10000000
            
                  allocate (ronp(ndro))
                  allocate (ropp(ndro))
                  allocate (ronh(ndro))
                  allocate (roph(ndro))
                  allocate (vint(ifmmx))

      
                  open(5,file=namev,status='unknown',form='unformatted') !,access='append')
            
            
            do ia=1,ifmmx
            
                  write(*,*)ia,ifmmx
            
                   jia=jphonm(ia)
                   ipa=iphonm(ia)
            
                   xfact=(dfloat(2*jia+1))**(-1.d0)
            
            
            
               if (jia.eq.jcal.and.ipa.eq.ipcal) then
                  
                   call readro2(namernp,ia,ronp,nronp)
if(idimnp(jia).ne.nronp) write(*,*) 'Warning idim ne nro for np',idimnp(jia),nronp
if(idimnp(jia).ne.nronp) stop
                   call readro2(namerpp,ia,ropp,nropp)
if(idimpp(jia).ne.nropp) write(*,*) 'Warning idim ne nro for np',idimpp(jia),nropp
if(idimnp(jia).ne.nronp) stop
                   call readro2(namernh,ia,ronh,nronh)
if(idimnh(jia).ne.nronh) write(*,*) 'Warning idim ne nro for np',idimnh(jia),nronh
if(idimnh(jia).ne.nronh) stop
                   call readro2(namerph,ia,roph,nroph)
if(idimph(jia).ne.nroph) write(*,*) 'Warning idim ne nro for np',idimph(jia),nroph
if(idimph(jia).ne.nroph) stop


              do isi=1,ifmx
                  jisi=jphon(isi)
                  ips=iphon(isi)
                  vint=0.d0

                  do ii=1,nropp
                        
                      if (c_pp(jia,ii)%j.eq.jisi) then 
                        ig=c_pp(jia,ii)%ilap
                        ipg=iphonm(ig)
                        if ((ips*ipg).eq.ipa) then 
                        i1=c_pp(jia,ii)%i1
                        i2=c_pp(jia,ii)%i2
                        ji1=levp(i1)%j
                        ji2=levp(i2)%j
                        iffz=(-1)**(jig-jia+(ji1-ji2)/2)
                         do jj=1,ndcamp(isi)
                          ip=camp(isi,jj)%par
                          ih=camp(isi,jj)%hol
               !c           jp=lev(ip)%j
                          campn=camp(isi,jj)%am
                          vint(ig)=vint(ig)+0.5d0*campn*fp(jisi,ip,ih,i1,i2)*ropp(ii)
                        enddo
               
                        do jj=1,ndcamn(isi)
                          ip=camn(isi,jj)%par
                          ih=camn(isi,jj)%hol
               !c           jp=lev(ip)%j
                          campn=camn(isi,jj)%am
                          vint(ig)=vint(ig)+campn*fpn(jisi,i1,i2,ip,ih)*ropp(ii)
                        enddo
                       endif
                      endif 
                  enddo


                  do ii=1,nronp

                  if (c_np(jia,ii)%j.eq.jisi) then 
                        ig=c_np(jia,ii)%ilap
                        ipg=iphonm(ig)
                        if ((ips*ipg).eq.ipa) then
                        i1=c_np(jia,ii)%i1
                        i2=c_np(jia,ii)%i2
                        ji1=levn(i1)%j
                        ji2=levn(i2)%j
                        iffz=(-1)**(jig-jia+(ji1-ji2)/2)
                        
                         do jj=1,ndcamn(isi)
                          ip=camn(isi,jj)%par
                          ih=camn(isi,jj)%hol
               !           jp=lev(ip)%j
                          campn=camn(isi,jj)%am
                          vint(ig)=vint(ig)+0.5d0*campn*fn(jisi,ip,ih,i1,i2)*ronp(ii)
               
                        enddo
               
                        do jj=1,ndcamp(isi)
                          ip=camp(isi,jj)%par
                          ih=camp(isi,jj)%hol
               !           jp=lev(ip)%j
                          campn=camp(isi,jj)%am
                          vint(ig)=vint(ig)+campn*fpn(jisi,ip,ih,i1,i2)*ronp(ii)
               
                        enddo
                      endif
                    endif
                  enddo


                  do ii=1,nroph
                   if (c_ph(jia,ii)%j.eq.jisi) then
                        ig=c_ph(jia,ii)%ilap
                        ipg=iphonm(ig) 
                        if ((ips*ipg).eq.ipa) then      
                        i1=c_ph(jia,ii)%i1
                        i2=c_ph(jia,ii)%i2
                        ji1=levp(i1)%j
                        ji2=levp(i2)%j
                        iffz=(-1)**(jig-jia+(ji1-ji2)/2)
                         do jj=1,ndcamp(isi)
                          ip=camp(isi,jj)%par
                          ih=camp(isi,jj)%hol
               !           jp=lev(ip)%j
                          campn=camp(isi,jj)%am
                          vint(ig)=vint(ig)+0.5d0*campn*fp(jisi,ip,ih,i1,i2)*roph(ii)
                   enddo
               
                        do jj=1,ndcamn(isi)
                          ip=camn(isi,jj)%par
                          ih=camn(isi,jj)%hol
               !           jp=lev(ip)%j
                          campn=camn(isi,jj)%am
                          vint(ig)=vint(ig)+campn*fpn(jisi,i1,i2,ip,ih)*roph(ii)
                        enddo
                    endif
                    endif
                  enddo

                  do ii=1,nronh
                    if (c_nh(jia,ii)%j.eq.jisi) then   
                        ig=c_nh(jia,ii)%ilap
                        ipg=iphonm(ig)
                        if ((ips*ipg).eq.ipa) then 
                        i1=c_nh(jia,ii)%i1
                        i2=c_nh(jia,ii)%i2
                        ji1=levn(i1)%j
                        ji2=levn(i2)%j
                        iffz=(-1)**(jig-jia+(ji1-ji2)/2)
                        
                        do jj=1,ndcamn(isi)
                          ip=camn(isi,jj)%par
                          ih=camn(isi,jj)%hol
               !           jp=lev(ip)%j
                          campn=camn(isi,jj)%am
                          vint(ig)=vint(ig)+0.5d0*campn*fn(jisi,ip,ih,i1,i2)*ronh(ii)
                        enddo
               
                        do jj=1,ndcamp(isi)
                          ip=camp(isi,jj)%par
                          ih=camp(isi,jj)%hol
               !           jp=lev(ip)%j
                          campn=camp(isi,jj)%am
                          vint(ig)=vint(ig)+campn*fpn(jisi,ip,ih,i1,i2)*ronh(ii)
                        enddo
                    endif 
                    endif
                  enddo

                  do ig=1,ifmmx
                      jig=jphonm(ig)  
                      ifaz=(-1)**(jia+jig+jisi)  
                      if (dabs(vint(ig)).gt.1.d-10) then 
                        write(5)ia,ig,isi,xfact*dfloat(ifaz)*vint(ig)
                        write(9991,*)ia,ig,isi,xfact*dfloat(ifaz)*vint(ig)
                      endif 
                  enddo
            
              enddo ! loop isi
                        
             endif
                  
          enddo ! loop ia
            
                  
            deallocate(camp,camn,jphon,jphonm,fp,fpn,ronp,ronh,ropp,roph)
            
            close(5)
            return
            end subroutine vintn23
            
            
            
!***************************************************************************
            
      subroutine vintn2

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_ndgi_int.inc'
      include 'input_ndgi_int.inc'
      include 'formats_ndgi_int.inc'


      double precision, dimension(:,:,:,:,:),allocatable :: fp,fn,fpn

      integer, dimension(:), allocatable :: jphon,jphonm,ironp,iropp,ironh,iroph,ndcamn,ndcamp

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
!      type(level_typ),dimension(*) :: lev

      character*10 namer
      character*30 namefp,namefn,namefpn,namecp,namecn,namerpp,namernp,namerph,namernh,namev,nameff
     

            
      ndamp=1600

      ifmx=1600
      ifmmx=24000
      allocate (jphon(ifmx))  ! 1phonon 
      jphon=0

      allocate (jphonm(ifmmx)) ! n-1 phonon
      jphonm=0



      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
      enddo

      close(3)

      ifmx=i
      ndlam=ifmx

     

              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
              namefp='fmat_p.dat'
              namefn='fmat_n.dat'
              namefpn='fmat_pn.dat'
              namev='Vint_phon2.dat'
 


!     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihpmn,ippmx,ihpmn,ippmx,fp)
      call readfin(namefn,jmin,jmax,ihnmn,ipnmx,ihnmn,ipnmx,fn)
!     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihpmn,ippmx,ihnmn,ipnmx,fpn)
!     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
!     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      
      write(*,*)'Calculation of redefined interaction'



      open(5,file=namev,status='unknown',form='unformatted')

      do ia=1,ifmx

      write(*,*)ia,ifmx

       jig=jphon(ia)

!      write(5)ia


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
!           jp=lev(ip)%j
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
!           jp=lev(ip)%j
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
!           jp=lev(ip)%j
           campn=camn(isi,jj)%am

         vint=vint+campn*campp*fpn(jisi,ipa,iha,ips,ihs)

         enddo

         enddo

          endif




        if (dabs(vint).gt.xrotrunc) write(997,*)ia,isi,vint

        if (dabs(vint).gt.xrotrunc) write(5)ia,isi,vint


      enddo ! loop isi
      enddo ! loop ib

      
      deallocate(camp,camn,jphon,fp,fpn)

!      write(5)10000000
!      write(5)0,0,0,0,0.d0 

      close(33)
      close(34)
      close(43)
      close(44)
      close(5)
      return
      end subroutine vintn2


!************************************************************************


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
      


      subroutine loadsp(levn,levp)

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


      call inp_sp(levn,levp)

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
      
      
