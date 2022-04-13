c       Cholesky    procedure         
c     last modification 13.5.2010

      module choleski

      use read_admat      

      contains
       
      subroutine cholesk(ndim,no,noo,dd,nx,mxt)

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'
     
      
      double precision, dimension(:,:), allocatable :: dd,cq,d1,cdu
    
      double precision, dimension(:), allocatable :: r,rcp 

      integer, dimension(:), allocatable :: nx,mxt


!      write(*,*)' Cholesky analysis of D matrix '

c      if (ndim.gt.nda) then
c              write(*,*)'Small dimensions of arrays in cholesk'
c              stop
c                        endif

 
      allocate(dd(ndim,ndim))
      dd=0.d0
      allocate(d1(noo,ndim))
      d1=0.d0
      allocate(cdu(ndim,noo))
      cdu=0.d0
      allocate(rcp(ndim+1),r(ndim),nx(ndim),mxt(ndim))
      rcp=0.d0
      r=0.d0
      nx=0
      mxt=0
      
      norder=ndim

      call read_dmat(dd)

c      do i=1,10
c      write(981,102)(dd(i,j),j=1,10)
c      enddo

 
      do im=1,ndim
 
      r(im)=dd(im,im)

c      write(109,*)im,im,r(im)
c      r(im)=dmat(nf,ib,ip,ih,ib,ip,ih) 
c       write(*,*)im,ib,ifon(nf,im),ip,ih           
      end do
      
c       write(87,*)(r(i),i=1,ndim)

c       xsum=0.d0
c       do i=1,ndim
c       xsum=xsum+r(i)
c       end do
c       write(*,*)'sum',xsum      
                      
      do kk=1,norder
        mxt(kk)=kk
        rcp(kk)=r(kk)
      enddo
c

      rtr=0.0001d0
      do 200 i=1,norder
      tst=r(i)
      tcp=rcp(i)
      numb=i
      if(i.eq.norder) go to 201
      do 210 j=i+1,norder
      if(r(j).gt.tst) then
      tst=r(j)
      tcp=rcp(j)
      numb=j
      endif
  210 continue
      its=mxt(numb)
      r(numb)=r(i)
      rcp(numb)=rcp(i)
      mxt(numb)=mxt(i)
      r(i)=tst
      rcp(i)=tcp
      mxt(i)=its
  201 continue
c     
      if(r(i).lt.rtr) go to 299
      r(i)=dsqrt(r(i))
c
      if(i.eq.1) then
      cdu(1,1)=r(1)
      mr=mxt(1)
c
c
      d1(1,mr)=rcp(1)
      do 220 j=2,norder
      ml=mxt(j)
c
c
c ****** ovu calculates elements D1(mr,ml)
c       ib=ifon(nf,ml+imin-1)
c       ibp=ifon(nf,mr+imin-1)
c       ip=int(iphbs(iph(nf,ml+imin-1))/10000)
c       ih=iphbs(iph(nf,ml+imin-1))-10000*ip
c       ipp=int(iphbs(iph(nf,mr+imin-1))/10000)
c       ihp=iphbs(iph(nf,mr+imin-1))-10000*ipp
c              write(*,*)'***',im,jm,ibp,ipp,ihp,ib,ip,ih  
  
c      d1(1,ml)=dmat(nf,ibp,ipp,ihp,ib,ip,ih)

c      tst1=dmat(tt,mr,mr)
c      if (tst1.lt.rtr) then  
c      d1(1,ml)=0.d0
c      else 
c      tst2=dmat(tt,ml,ml)
c      if (tst2.lt.rtr) then
c      d1(1,ml)=0.d0
c      else            
      d1(1,ml)=dd(mr,ml)
c      endif
c      endif
      
c      write(109,*)mr,ml,d1(1,ml) 
c 
      cdu(j,1)=d1(1,ml)/cdu(1,1)
  220 r(j)=r(j)-cdu(j,1)**2
c
      else
c
      cdu(i,i)=r(i)
      if(i.eq.norder) go to 200
      mr=mxt(i)
c
c
      d1(i,mr)=rcp(i)
      do 230 l=1,i-1
      tsl=cdu(numb,l)
      cdu(numb,l)=cdu(i,l)
      cdu(i,l)=tsl
      do 240 j=i+1,norder
      if(l.eq.1) then
      ml=mxt(j)
c
c
c ****** ovu calculates elements D1(mr,ml)
c       ib=ifon(nf,ml+imin-1)
c       ibp=ifon(nf,mr+imin-1)
c       ip=int(iphbs(iph(nf,ml+imin-1))/10000)
c       ih=iphbs(iph(nf,ml+imin-1))-10000*ip
c       ipp=int(iphbs(iph(nf,mr+imin-1))/10000)
c       ihp=iphbs(iph(nf,mr+imin-1))-10000*ipp

c      d1(i,ml)=dmat(nf,ibp,ipp,ihp,ib,ip,ih)
c      d1(i,ml)=dmat(tt,mr,ml)
c
c      tst1=dmat(tt,mr,mr)
c      if (tst1.lt.rtr) then  
c      d1(i,ml)=0.d0
c      else 
c      tst2=dmat(tt,ml,ml)
c      if (tst2.lt.rtr) then
c      d1(i,ml)=0.d0
c      else            
      d1(i,ml)=dd(mr,ml)
c      endif
c      endif






      
c      write(109,*)mr,ml,d1(i,ml)
c
     
c
      cdu(j,i)=d1(i,ml)
      endif
  240 cdu(j,i)=cdu(j,i)-cdu(i,l)*cdu(j,l)
  230 continue
      do 250 j=i+1,norder
      cdu(j,i)=cdu(j,i)/cdu(i,i)
  250 r(j)=r(j)-cdu(j,i)**2
c
      endif
c
  200 continue
c
  299 if (i-1.ne.ndim) then  ! toto som tu upravil pre pripad ndim=no
      its=mxt(i)
      mxt(i)=mxt(numb)
      mxt(numb)=its
      endif
      NO=I-1
c
      write(*,*)' dimension = ',norder
c      write(6,2500) no     
c 2500 format(1x,' no =',i6)
c -----------------------
c -----------------------
c      write(21) itb,jb,ipb,no
c      write(22) itb,jb,ipb,no,norder
      if(no.eq.0) go to 3000
c      write(22) ((ivec(l,m),l=1,norder),m=1,2)
c
c *** calcolates the inverse of the metric matrixc
c



      goto 6666
      write(*,*)'Calculation of D^-1 '
 
      DO 3270 I=1,NO
      cdu(I,I)=1.d0/cdu(I,I)
      IF(I.EQ.NO) GO TO 3270
      IP1=I+1
      DO 3260 J=IP1,NO
      TST=0.d0
      JM1=J-1
      DO 3250 K=I,JM1
 3250 TST=TST-cdu(J,K)*cdu(K,I)/cdu(J,J)
      cdu(J,I)=TST
 3260 CONTINUE
 3270 CONTINUE
C
c      rprec=10.d0**(-5)
      rprec=1.d0
      DO 3290 I=1,NO
      DO 3290 J=1,NO
      TEMP=0.d0
      DO 3280 K=J,NO
 3280 TEMP=TEMP+rprec*cdu(K,I)*cdu(K,J)
 3290 cdu(I,J)=TEMP
      DO 3330 I=1,NO
      DO 3333 J=I,NO
      cdu(I,J)=cdu(I,J)
 3333 cdu(J,I)=cdu(I,J)
 3330 CONTINUE


 6666 continue     

      do il=1,no
      nx(mxt(il))=1
      do jl=il,norder
      r(jl)=d1(il,mxt(jl))
      enddo
      do jl=il,norder
      d1(il,jl)=r(jl)
      enddo
      enddo
c
      do il=1,no
      do jl=il,no
      d1(jl,il)=d1(il,jl)
c      cq(il,jl)=cdu(il,jl)
c      cq(jl,il)=cdu(jl,il)
      enddo
      enddo
c
c      open(77,file='chol.o',status='unknown',form='formatted')
      
 3000 if(no.eq.0) then
c      write(77,10)
c   10 format(/,1x,' ----- this space has zero dimension -----') 
      endif
      
c       write(77,*)norder
c       write(77,*)(mxt(i),i=1,norder)
c       write(77,*)no
c       write(77,*)(nx(i),i=1,norder)
 
c       close(77)
  
   66 format(500f10.5) 
c      do i=1,norder
c        write(77,66)(d1(i,j),j=1,norder)
c      end do
      continue


c          test preusporiadania
c          do i=1,no
c          do j=1,ndim
c          ii=mxt(i)
c          jj=mxt(j)
c          if (dabs(d1(i,j)-xtest(ii,jj)).gt.1d-8) 
c     *write(*,*)'warning:bad reordering of D',i,j   
c          if (d1(i,j).ne.d1(j,i).and.(j.le.no)) 
c     *write(*,*)'warning: matrix D1 is not symmetric!!!'   
c          end do      
c          end do


c                 END Cholesky      
 
c 6666 continue
      deallocate(rcp,r)

!      open(6,file='d1matr.dat',status='unknown',form='unformatted')
!      do i=1,no
!       do j=1,ndim
!        write(6)i,j,d1(i,j)
!       enddo
!      enddo

!      close(6)


      open(6,file='mxt.dat',status='unknown',form='unformatted')

      write(6)ndim,no
      do i=1,ndim
        write(6)i,mxt(i)
      enddo
      close(6)

      write(*,*)' Number of linearly independent states ',no


      deallocate(cdu)
 
      return 
      end subroutine cholesk

      end module choleski
     
