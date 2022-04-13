!c     last modification 30.5.2013

      module hami

      contains

      subroutine ham(ndim,ndimr,no,nor,irow,wr,xr,vr,ipar,jcal,mxtr)

      use choleski

      implicit double precision (a-h,o-z)

      include 'formats_eqm.inc'

      double precision, dimension(:,:), allocatable :: d1,amatr,hami
     *,hamid,cq,vr,d1r,hamir,hamidr,xr
     
      double precision, dimension(:), allocatable ::  work,wr,wi,wro

      integer, dimension(:), allocatable :: mxt,nxt,mxtr,ipoz,irow



!c      allocate(mxt(ndimr))
!c      mxt=0


!c      open(6,file='mxt.dat',status='old',form='unformatted')

!c      ii=0
!c      do while (.not.eof(6))
!c      read(6)i,mm
!c       ii=ii+1
!c       mxt(i)=mm
!c      enddo

!c      write(*,*)' Number of independent states ',ii

!c      close(6)


!c      allocate(xr(ndim,no))
!c      xr=0.d0
!c      allocate(wr(no))
!c      wr=0.d0

      allocate(d1(ndimr,ndim))
      d1=0.d0

      allocate(amatr(ndimr,ndim))
      amatr=0.d0

      ndimtotal=ndim


      open(6,file='d_mat.dat',status='old',form='unformatted')
      read(6)ndimrt,ndimt
      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
        write(*,*)' Dimensions does not match in D_m file'
        stop
      endif
      do iii=1,ndimrt
       read(6)(d1(iii,jjj),jjj=1,ndimt)
      enddo
      
!      do while (.not.eof(6))
!      read(6)i,j,dd
!       d1(i,j)=dd
!      enddo



      close(6)

      open(6,file='a_mat.dat',status='old',form='unformatted')

!      do while (.not.eof(6))
!      read(6)i,j,dd
!       amatr(i,j)=dd
!      enddo

      read(6)ndimrt,ndimt
      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
        write(*,*)' Dimensions does not match in A_m file'
        stop
      endif
      do iii=1,ndimrt
       read(6)(amatr(iii,jjj),jjj=1,ndimt)
      enddo

      close(6)

      write(*,*)ndim,ndimr


      iout=0
      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** matrix  A *************'       
      write(998,*)
      do i=1,no
        write(998,102)(amatr(i,j),j=1,ndim)
      enddo
      endif

      allocate(hami(ndimr,ndimr))
      hami=0.d0

c      do i=1,no
c         ii=i
c        do j=1,no
c           jj=mxt(j)
c           hh=0.d0
c          do k=1,ndim
c            kk=k
c            hh=hh+amatr(ii,kk)*d1(j,k)
c          enddo
c            hami(i,j)=hh
c        enddo 
c      enddo
      
       call dgemm('N','T',ndimr,ndimr,ndim,1.d0,amatr,ndimr,d1
     *,ndimr,0.d0,hami,ndimr)

      iout=0

      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** D *************'       
      write(998,*)
      do i=1,ndim
        write(998,102)(d1(j,i),j=1,no)
      enddo
      endif

      deallocate(d1)


      iout=0
      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** matrix  AD *************'       
      write(998,*)
      do i=1,ndimr
        write(998,102)(hami(i,j),j=1,ndimr)
      enddo
      endif

      do i=1,ndimr
c       write(777,*)i,hami(i,i)
       do j=1,ndimr
          if (dabs(hami(i,j)-hami(j,i)).gt.0.0001d-0) then
           write(*,*)'Non-symmetric AD'
           write(*,'(2i5,3f10.5)')i,j,hami(i,j),hami(j,i),hami(i,j)-hami(j,i)
          endif   
     
       enddo
      enddo


      allocate(mxt(ndimr))
      mxt=0

      open(6,file='mxt.dat',status='old',form='unformatted')

      read(6)ndimm,no

      if (ndimm.ne.ndimr) then
              write(*,*)' Dimensions does not match '

      endif

      do while (.not.eof(6))
      read(6)i,mm
       mxt(i)=mm
      enddo

      close(6)

      allocate(nxt(ndimr))
      nxt=0

      open(6,file='nxt.dat',status='unknown',form='unformatted')
      write(6)ndimm,no

      do i=1,no
        ii=mxt(i)
        nxt(ii)=i
      enddo

      write(6)(nxt(i),i=1,ndimr)
      close(6)

      deallocate(nxt)


      write(*,*)' Loaded number of states ',ndimm,no

c      write(992,*)(mxt(jj),jj=1,ndimm)

      allocate(d1(ndimr,ndimr))
      d1=0.d0

      call read_dmat(d1)

      iout=0
      if (iout.eq.1) then

      write(998,*)
      write(998,*)'******** D1 *************'
      write(998,*)
      do i=1,ndimr
      write(998,102)(d1(j,i),j=1,ndimr)
      enddo
      endif

      call permutuj(ndimr,d1,mxt)
      call permutuj(ndimr,hami,mxt)



      if (iout.eq.1) then

      write(998,*)
      write(998,*)'******** perm D1 *************'
      write(998,*)
      do j=1,ndimr
      write(998,*)j,mxt(j)
      enddo
      write(998,*)
      do i=1,ndimr
      write(998,102)(d1(j,i),j=1,ndimr)
      enddo
      endif

      allocate(d1r(no,no))
      d1r=0.d0

      do i=1,no
       do j=i,no
        d1r(i,j)=d1(i,j)
       enddo
      enddo

      deallocate(d1)


c      allocate(hamir(nor,nor))
c      hamir=0.d0
 
c      call redmat(irow,hamir,hami,no,nor)

c      do i=1,nor
c       do j=1,nor
c         hamir(i,j)=hami(i,j)
c       enddo
c      enddo

c      deallocate(hami)

     

c      if (iout.eq.1) then 
c      write(998,*)    
c      write(998,*)'******** reduced matrix  AD *************'       
c      write(998,*)
c      do i=1,nor
c        write(998,102)(hamir(i,j),j=1,nor)
c      enddo
c      endif

c      allocate(d1r(nor,nor))
c      d1r=0.d0


 
c      call redmat(irow,d1r,d1,no,nor)

c      do i=1,nor
c       do j=1,nor
c         d1r(i,j)=d1(i,j)
c       enddo
c      enddo

c      deallocate(d1)

c      if (iout.eq.1) then 
c      write(998,*)    
c      write(998,*)'******** reduced D *************'       
c      write(998,*)
c      do i=1,nor
c      write(998,102)(d1r(j,i),j=1,nor)
c      enddo
c      endif


      nor=no

      call dpotrf('U',nor,d1r,nor,info)
      write(*,*)' Factorization info ',info

      

      call dpotri('U',nor,d1r,nor,info)

      write(*,*)' Inverse info ',info


      do i=1,nor
       do j=i+1,nor
         d1r(j,i)=d1r(i,j)
       enddo
      enddo


      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** D^-1 *************'       
      write(998,*)
      do i=1,nor
      write(998,102)(d1r(j,i),j=1,nor)
      enddo
      endif


      

c     calculation of D^-1 (AD)

c      allocate(cq(no,no))
c      cq=0.d0

c      open(6,file='cq.dat',status='old',form='unformatted')

c      do while (.not.eof(6))
c      read(6)i,j,dd
c       cq(i,j)=dd
c      enddo

c      close(6)

c      allocate(hamid(no,no))
c      hamid=0.d0



      allocate(hamir(nor,nor))
      hamir=0.d0

      do i=1,nor
       do j=1,nor
         hamir(i,j)=hami(i,j)
       enddo
      enddo

      deallocate(hami)

      allocate(hamidr(nor,nor))
      hamidr=0.d0
      
      write(*,*)'calculating D^-1 AD'

c      call dgemm('N','N',no,no,no,1.d0,cq,no,hami,no,0.d0,hamid,no)

      call dgemm('N','N',nor,nor,nor,1.d0,d1r,nor,hamir,nor,0.d0,
     *hamidr,nor)

      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** D^-1 AD *************'       
      write(998,*)
      do i=1,nor
      write(998,102)(hamidr(j,i),j=1,nor)
      enddo
      endif

      deallocate(d1r,hamir)


      

c     diagonalisation

      no=nor
       
      write(*,*)' Diagonalisation '
      lwork=20*no
      lda=no
      ldvl=no
      ldvr=no
      allocate(work(lwork),wi(no),wr(no),vr(no,no),wro(no))
      wr=0.d0
      wro=0.d0
      wi=0.d0
      work=0.d0
      vr=0.d0

      if (nor.gt.0) then

      CALL DGEEV('N','V',no,hamidr,lda,wr,wi,vl,ldvl,vr,
     *ldvr,work,lwork,info)


        write(*,*)' info=  ',info

        do i=1,no
       if (dabs(wi(i)).gt.1.d-10) write(*,*)' Imaginary part',i,wi(i)
        enddo
        
        allocate(ipoz(no))


      do i=1,no
        ipoz(i)=i
        wro(i)=wr(i)
      enddo

      do i=1,no
        do j=1,no-i
          if (wro(j).gt.wro(j+1)) then
            xe=wro(j)
            wro(j)=wro(j+1)
            wro(j+1)=xe
            ipozz=ipoz(j)
            ipoz(j)=ipoz(j+1)
            ipoz(j+1)=ipozz
      endif
      enddo
      enddo

        write(99,*)' Parity = ',ipar, ' J = ',jcal
        write(99,*)
        write(99,*)(wro(i),i=1,no) 
        write(99,*)


      endif

!      deallocate (work,wi,wro)


      ndim=ndimr

      allocate(d1(no,ndimtotal))
      d1=0.d0
   
      call read_dmatch(d1)

c      open(6,file='dmat.dat',status='old',form='unformatted')
c      do while (.not.eof(6))
c      read(6)i,j,dd
c       d1(i,j)=dd
c      enddo
c      close(6)

c      call dgemm('N','T',no,no,ndim,1.d0,amatr,no,d1,no,0.d0,hami,no)

     
      write(*,*)' X  ndim, no, ndimtotal',ndim,no,ndimtotal

      allocate(xr(ndimtotal,no))
      xr=0.d0


      call dgemm('T','N',ndimtotal,no,no,1.d0,d1,no,vr,no,0.d0,
     *xr,ndimtotal)

      iout=0
      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** X nenorm. *************'       
      write(998,*)
      do i=1,ndimtotal
      write(998,102)(xr(i,j),j=1,no)
      enddo
      endif




      call normal(no,ndimtotal,vr,xr,jcal,mxtr,mxt)


      iout=0
      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** X *************'       
      write(998,*)
      do i=1,ndimtotal
      write(998,102)(xr(i,j),j=1,no)
      enddo
      endif


      write(*,*)'**** OK '

c      deallocate(d1)

!      deallocate (work,wi,wro)
            
      end subroutine ham

      subroutine redmat(irow,amatrp,amatr,ndim,ndimr)

      implicit double precision (a-h,o-z)

      double precision, dimension(:,:), allocatable ::  amatr,amatrp
      integer, dimension(:), allocatable :: irow

      ii=0
      do i=1,ndim
      if (irow(i).ne.0) ii=ii+1
      enddo

      ndimr=ii

c      allocate(amatrp(ndimr,ndimr))
c      amatrp=0.d0

      ii=0
      do i=1,ndim
       if (irow(i).ne.0) then
        ii=ii+1
        jj=0       
        do j=1,ndim
        if (irow(j).ne.0) then 
          jj=jj+1
          amatrp(ii,jj)=amatr(i,j)
        endif
        enddo

      endif
      enddo


      endsubroutine redmat



      subroutine kicksp(ispu,irow,phonbase,phon1,phon2,mxt,ndim,ndimr)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      type(phonbase_typ), dimension (:), allocatable :: phonbase
      type(phon_typ), dimension (:), allocatable :: phon1,phon2

      integer, dimension (:), allocatable :: irow,mxt

      allocate(irow(ndim))
      irow=0

      iii=0

      do i=1,ndim
       ii=mxt(i)
       il=phonbase(ii)%ila
       ilp=phonbase(ii)%ilap
       e1=phon2(il)%enf
       e2=phon1(ilp)%enf

       if (i.eq.i) then
c       if (e1.lt.10.d0.or.e2.lt.10.d0) then 

c       if (i.eq.2*int(i/2)) then 

c       if (il.ne.ispu.and.ilp.ne.ispu) then 
                                 irow(i)=1
                                 iii=iii+1
                         endif

      enddo

      ndimr=iii
     
      endsubroutine kicksp

******************************************************************************      
      subroutine normal(no,ndmx,vr,xc,jcal,mxtr,mxt)

      implicit double precision (a-h,o-z)

c      include 'chole.inc'

      double precision, dimension(:,:), allocatable :: vr,xc,d1
      integer, dimension (:), allocatable :: mxtr,mxt

      xfact=(dfloat(2*jcal+1))**0.5d0
c      xc=0.d0

c      do i=1,no
c        do j=1,ndmx
c          xpom=0.d0
c           do k=1,no
c             xpom=xpom+d1(k,j)*vr(k,i)
c           enddo
c             xc(j,i)=xpom
c         enddo
c      enddo

      do i=1,no
        xpom=0.d0
        do j=1,no
          jj=mxtr(mxt(j))
          xpom=xpom+vr(j,i)*xc(jj,i)
         enddo
        do j=1,no
          vr(j,i)=vr(j,i)/dsqrt(xpom)
        enddo
        do j=1,ndmx
          xc(j,i)=xfact*xc(j,i)/dsqrt(xpom)
        enddo
        write(771,*)' i =',i,xpom
      enddo

      return
      end subroutine normal
**************************************************************
      subroutine permutuj(ndim,d1,mxt)

      implicit double precision (a-h,o-z)

      double precision, dimension(:,:), allocatable :: d1
      double precision, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: mxt,nxt

      allocate (work(ndim))
      work=0.d0
      allocate(nxt(ndim))

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=d1(i,:)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then 
        d1(i,:)=d1(j,:)
        d1(j,:)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
       endif
       enddo
      enddo

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=d1(:,i)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then 
        d1(:,i)=d1(:,j)
        d1(:,j)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
        endif
       enddo
      enddo


 
      deallocate(work)

      return
      end subroutine permutuj


!***************************************************************************

      subroutine read_dmatch(dmatr)

      implicit double precision (a-h,o-z)


      double precision, dimension(:,:), allocatable :: dmatr
      double precision, dimension(:), allocatable :: dmm
      integer, dimension(:), allocatable :: nxtr

      
      dmatr=0.d0

      open(2,file='nxt.dat',status='unknown',form='unformatted')
      read(2)ndim,no
      allocate(nxtr(ndim))
      read(2)(nxtr(i),i=1,ndim)
      close(2)

      write(911,*)' idphon =', ndim
      write(911,*)(nxtr(i),i=1,ndim)


      open(66,file='d_mat.dat',status='old',form='unformatted')

      read(66)ndimrt,ndimt
      allocate(dmm(ndimt))
      dmm=0.0d0


!      do while (.not.eof(66))
!      read(66)i,j,dd
!      if (nxtr(i).ne.0) then 
!      dmatr(nxtr(i),j)=dd
!      endif

!      enddo

!      close(66)

      do iii=1,ndimrt
       read(66)(dmm(jjj),jjj=1,ndimt)
       do jjj=1,ndimt
         if (nxtr(iii).ne.0) then
          dmatr(nxtr(iii),jjj)=dmm(jjj)
         endif
       enddo
      enddo

      close(66)
      deallocate(dmm)


      deallocate(nxtr)

      end subroutine read_dmatch

      

      end module hami
