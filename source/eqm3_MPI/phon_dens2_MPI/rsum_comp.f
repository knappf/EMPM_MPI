      program runsum
      parameter (nf=100000)    

      implicit double precision (a-h,o-z)

      dimension be(0:nf),en(0:nf),rsum(0:nf)
      
      character(len=20) filein     


101   format(i5,2f20.10)
102   format(i5,3f20.10)

      write(*,*) 'Name of file with BE values:'
      read(*,'(A )') filein

      write(*,*)'Energy of ground state in MeV :'
      read(*,*)egs

      ir=0
      open(1,file=filein,form='formatted',status='old')
   
      do while (.not.eof(1))
      read(1,101)irt,bet,ent
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
      enddo
    
      open(1,file='rsum_'//filein,form='formatted',status='unknown')
    
      do i=0,ir
      write(1,102)i,en(i)-egs,be(i),rsum(i)/sumrule
      enddo
    
      close(1)


 998  write(*,*)'width of lorentzian ?'
      read(*,*)delta
      write(*,*)'ilambda, ilambda2f ?'
      read(*,*)ilambda,ilambda2f
      open(1,file='silf_'//filein,status='unknown',form='formatted')
      open(2,file='crsec'//filein,status='unknown',form='formatted')
      emin=egs
      emax=0.d0
      ibod=100000
      xe=emin
      do j=1,ibod
      estep=(emax-emin)/ibod
      xe=xe+estep
      xsil=0.d0
      do ii=1,ir
      xsil=xsil+be(ii)*xlor(xe-en(ii),delta)
      enddo
      write(1,*)xe-emin,xsil
      xcr=xsil*(xe-emin)**dfloat(2*ilambda-1)*(ilambda+1)
     */(197.33d0**(dfloat(2*ilambda))*ilambda*ilambda2f)
      write(2,*)xe-emin,xcr

      enddo
      close(1)

      goto 998


      
      end

********************************************************************
      double precision function xlor(en,delt)

      implicit double precision (a-h,o-z)

      Pi=3.141592653589793D0
      xpom=delt/2.d0*Pi
      xpom=xpom/(en**2.d0+delt**2.d0/4.d0)

      xlor=xpom

      return
      end

*******************************************************************

