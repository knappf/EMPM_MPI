!     update 11.8. 2011
      program runsum
      parameter (nf=100000)    

      implicit double precision (a-h,o-z)

      dimension be(0:nf),en(0:nf),rsum(0:nf),rsumnew(0:nf)
c      dimension sil(1:nf)
      
      character(len=25) filein     


101   format(i5,2f20.10)
102   format(i5,3f20.10)

      write(*,*) 'Name of file with BE values:'
      read(*,'(A )') filein

      write(*,*)'Energy of ground state in MeV :'
      read(*,*)egs

      ir=0
      open(1,file=filein,form='formatted',status='old')
  
      read(1,*)
 
      do while (.not.eof(1))
      read(1,*)ent,bet
      ir=ir+1
      be(ir)=bet
      en(ir)=ent
c      read(1,101)ir,be(ir),en(ir)
      enddo

      close(1)
      write(*,*)'ir =',ir
        
c      do i=1,ir
c      ipoz(i)=i
c      enddo

      
      do i=1,ir
      do j=1,ir-i
      if (en(j).gt.en(j+1)) then
         xe=en(j)
         xbe=be(j)
         en(j)=en(j+1)
         be(j)=be(j+1)
         en(j+1)=xe
         be(j+1)=xbe
c         ipozz=ipoz(j)
c         ipoz(j)=ipoz(j+1)
c         ipoz(j+1)=ipozz
      endif
      enddo
      enddo
     

      write(*,*)' sum rule value'
      read(*,*)sumrule

      rsum(0)=0.d0
      en(0)=egs
      do i=1,ir
      rsum(i)=rsum(i-1)+be(i)*(en(i)-egs)
      rsumnew(i)=rsumnew(i-1)+be(i)
      enddo
    
      open(1,file='rsum_'//filein,form='formatted',status='unknown')
      open(11,file='rsumnew_'//filein,form='formatted',status='unknown')

    
      do i=0,ir
      write(1,102)i,en(i)-egs,be(i),rsum(i)/sumrule
      write(11,102)i,en(i)-egs,be(i),rsumnew(i)
      enddo
    
      close(1)

      write(*,*)'Centroid ',rsum(ir-1)/rsumnew(ir-1)



 998  write(*,*)'width of lorentzian ?'
      read(*,*)xdeltar
!      write(*,*)'ilambda, ilambda2f ?'
!      read(*,*)ilambda,ilambda2f
!     for E1
!       xfact=140530.d0
       xfact=70265.d0
       ilambda=1
       ilambda2f=9      ! [(2*ilambda+1)!!]^2

      open(1,file='silf_'//filein,status='unknown',form='formatted')
      open(2,file='crsec'//filein,status='unknown',form='formatted')
      emin=egs
      emax=100.d0
      ibod=100000
      xe=emin
      xthres=20.0d0
      write(*,*)'Energy threshold for wider width'
      read(*,*)xthres
      do j=1,ibod
      estep=(emax-emin)/ibod
      xe=xe+estep
 
      xsil=0.d0
      do ii=1,ir
      xdelta=xdeltar
!      xdelta=xdeltar*((en(ii)-emin))**3.0d0/5000.0d0
!      write(*,*)xdelta
      if ((en(ii)-egs).lt.xthres) xdelta=0.02d0     
      xsil=xsil+be(ii)*xlor(xe-en(ii),xdelta)
      enddo
      write(1,*)xe-emin,xsil
!c      xcr=xsil*xfact*(xe-emin)**dfloat(2*ilambda-1)*(ilambda+1)
!c     */(197.33d0**(dfloat(2*ilambda))*ilambda*ilambda2f)
     
      crfactor=0.28d0*10.d0 ! [in milibarns]
      xcr=xsil*crfactor*(xe-emin)

      write(2,*)xe-emin,xcr  ! in mb

      enddo
      close(1)


      open(1,file='step_'//filein,status='unknown',form='formatted')

      sil=0.0d0
      emin=egs
      emax=100.d0
      ebin=0.02d0
      ibod=int((emax-emin)/ebin)
      xe=emin

      do j=1,ibod
c       estep=(emax-emin)/ibod
       

      sil=0.0d0

      do ii=1,ir

      if ((en(ii)-egs).ge.xe.and.(en(ii)-egs).lt.(xe+ebin)) then    
        sil=sil+be(ii)
c        xsil=xsil+be(ii)*xlor(xe-en(ii),xdelta)
      endif


      enddo

      xe=xe+ebin

      write(1,*)xe-emin+ebin/2.0d0,sil
     


      enddo

      close(1)


c      goto 998


      
      end

********************************************************************
      double precision function xlor(en,xdelt)

      implicit double precision (a-h,o-z)

      xPi=3.141592653589793D0
      xpom=xdelt/(2.d0*xPi)
      xpom=xpom/(en**2.d0+xdelt**2.d0/4.d0)

      xlor=xpom

      return
      end

*******************************************************************

