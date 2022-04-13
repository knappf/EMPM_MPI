!c     last modification 11.5.2018

module hami

use types_eqm
use read_admat

contains

subroutine ham(ndim,ndimr,no,nor,ns,irow,wr,xr,vr,ipar,jcal,mxtr,phonbs,nx,h_corr)

      use choleski
      use cm_ort_svd

      implicit double precision (a-h,o-z)
       
!      include 'types_eqm.inc'
      include 'formats_eqm.inc'

      type(phonbase_typ), dimension (:), allocatable :: phonbs
      double precision, dimension(:,:), allocatable :: d1,amatr,hami,hamid,dmatr,cq,vr,d1r,hamir,hamidr,hamd,xr
      double precision, dimension(:,:), allocatable :: d_orig,h_orig,h_corr
     
      double precision, dimension(:), allocatable ::  work,wr,wi,wro

      integer, dimension(:), allocatable :: mxt,nxt,mxtr,ipoz,irow,nx

      double precision, dimension (:), allocatable :: s
      double precision, dimension (:,:), allocatable :: vt
      integer, dimension (:), allocatable :: ind_red
      logical :: decoup_cm


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

      write(998,*)'-------------- parity =',ipar,' -- J=',jcal,'------------'  

      allocate(d1(ndimr,ndim))
      d1=0.d0

      allocate(amatr(ndimr,ndim))
      amatr=0.d0

      ndimtotal=ndim


!      open(6,file='d_mat.dat',status='old',form='unformatted')
!      read(6)ndimrt,ndimt
!      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
!        write(*,*)' Dimensions does not match in D_m file'
!        stop
!      endif
!      do iii=1,ndimrt
!       read(6)(d1(iii,jjj),jjj=1,ndimt)
!      enddo

      call read_admatr('./scratch/d_mat_',d1,ndimr,ndim)

      iout=0
      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** matrix  D *************'       
      write(998,*)
      do i=1,ndimr
        write(998,102)(d1(i,j),j=1,ndim)
      enddo
      endif

      
!      do while (.not.eof(6))
!      read(6)i,j,dd
!       d1(i,j)=dd
!      enddo



!      close(6)

     
!      open(6,file='a_mat.dat',status='old',form='unformatted')


!      read(6)ndimrt,ndimt
!      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
!        write(*,*)' Dimensions does not match in A_m file'
!        stop
!      endif
!      do iii=1,ndimrt
!       read(6)(amatr(iii,jjj),jjj=1,ndimt)
!      enddo

!      close(6)

      call read_admatr('./scratch/a_mat_',amatr,ndimr,ndim)

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

!c      do i=1,no
!c         ii=i
!c        do j=1,no
!c           jj=mxt(j)
!c           hh=0.d0
!c          do k=1,ndim
!c            kk=k
!c            hh=hh+amatr(ii,kk)*d1(j,k)
!c          enddo
!c            hami(i,j)=hh
!c        enddo 
!c      enddo
      
      call dgemm('N','T',ndimr,ndimr,ndim,1.d0,amatr,ndimr,d1,ndimr,0.d0,hami,ndimr)

      iout=0

      if (iout.eq.1) then 
      write(998,*)    
      write(998,*)'******** D  *************'       
      write(998,*)
      do i=1,ndim
        write(998,102)(d1(j,i),j=1,no)
      enddo
      endif

!      deallocate(d1)


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
       do j=1,ndimr
          if (dabs(hami(i,j)-hami(j,i)).gt.0.001d0) then
           write(*,*)'Non-symmetric AD'
           write(*,'(2i5,3f10.5)')i,j,hami(i,j),hami(j,i),hami(i,j)-hami(j,i)
          endif   
     
       enddo
      enddo

!     SVD part 
      call dmat_spur_set(ndim,ndimr,no,phonbs,d1,mxtr,nx,s,vt,ns)

    allocate(dmatr(ndimr,ndimr))
    dmatr=0.0d0 

    do i=1,ndimr
       do j=1,ndimr
         dmatr(i,j)=d1(i,mxtr(j))
       enddo
    enddo
      
     call reduce_mat(nx,hami,ndimr,no)
     call reduce_mat(nx,dmatr,ndimr,no)

     iout=0
     if (iout.eq.1) then 
     write(998,*)    
     write(998,*)'******** matrix  D after Choleski *************'       
     write(998,*)
     do i=1,no
       write(998,102)(dmatr(i,j),j=1,no)
     enddo

     write(998,*)    
     write(998,*)'******** matrix  AD after Choleski *************'       
     write(998,*)
     do i=1,no
       write(998,102)(hami(i,j),j=1,no)
       
     enddo
     endif

!  copies of original AD and D matrices      
    allocate(h_orig(no,no),d_orig(no,no))
    h_orig=hami
    d_orig=dmatr


      
! transformation of reduced AD
     deallocate(amatr)
     allocate(amatr(no,no))
     amat=0.0d0
   
     call dgemm('N','T',no,no,no,1.d0,hami,no,vt,no,0.d0,amatr,no)
     hami=0.d0
     call dgemm('N','N',no,no,no,1.d0,vt,no,amatr,no,0.d0,hami,no)

   
   do i=1,no
    do j=i,no
!      if (dabs(hami(i,j)-hami(j,i)).gt.0.01d0) write(*,'(2i5,3f10.5)'),i,j,dabs(hami(i,j)-hami(j,i))
    enddo
   enddo

! transformation of reduced D
!   allocate(dmatc_orig(no,no))
!   dmatc_orig=d1
   
   ! deallocate(amat)
   !  allocate(amat(dim_ind,dim_ind))
     amatr=0.0d0
   
     call dgemm('N','T',no,no,no,1.d0,dmatr,no,vt,no,0.d0,amatr,no)
     dmatr=0.d0
     call dgemm('N','N',no,no,no,1.d0,vt,no,amatr,no,0.d0,dmatr,no)



! decoupling spurious subspace     

     decoup_cm=.true.
     if (decoup_cm.eq..TRUE.) then 
            do i=1,ns
          hami(i,i)=hami(i,i)+100000000.0d0
!          dmatr(i,i)=1.0d0
!             do j=ns+1,no
             do j=ns+1,no
!              hami(i,j)=0.0d0
!              hami(j,i)=0.0d0 
!              dmatr(i,j)=0.0d0
!              dmatr(j,i)=0.0d0
             enddo

!             hami(i,i)=1000000.d0
!             dmatr(i,i)=1.0d0
           enddo         
      endif
   
   iout=0
        if (iout.eq.1) then
         write(998,*)
         write(998,*)'******** transfromed matrix  D *************'
         write(998,*)
         do i=1,no
           write(998,'(1000f11.6)')(dmatr(i,j),j=1,no)
         enddo
        endif

        do i=1,no
          do j=i,no
        
            if (dabs(dmatr(i,j)-dmatr(j,i)).gt.0.01d0) then 
            write(*,*)' Transformed D asymmetric!'
            write(*,'(2i5,3f10.5)')i,j,dabs(dmatr(i,j)-dmatr(j,i))
            endif
          enddo
         enddo


       iout=0
        if (iout.eq.1) then
         write(998,*)
         write(998,*)'******** transfromed matrix  AD *************'
         write(998,*)
         do i=1,no
           write(998,'(1000f11.6)')(hami(i,j),j=1,no)
         enddo
        endif      
   

   !  Generalized EGV problem
   
   ! inverse
         call dpotrf('U',no,dmatr,no,info)
         write(*,*)' Factorization info ',info
   
         call dpotri('U',no,dmatr,no,info)
   
         write(*,*)' Inverse info ',info
         do i=1,no
           do j=i+1,no
             dmatr(j,i)=dmatr(i,j)
           enddo
          enddo
   
   ! D^-1 AD
   !      deallocate(dmat)
         allocate(hamd(no,no))
         call dgemm('N','N',no,no,no,1.d0,dmatr,no,hami,no,0.d0,hamd,no)
   
   !  diag
         write(*,*)' Diagonalisation '
         lwork=20*no
         allocate(work(lwork),wi(no),wr(no),vr(no,no),wro(no))
         wr=0.d0
         wro=0.d0
         wi=0.d0
         work=0.d0
         vr=0.d0
   
!         if (no.gt.0) then
   
         CALL DGEEV('N','V',no,hamd,no,wr,wi,vl,no,vr,no,work,lwork,info)
         write(*,*)' info=  ',info
   
   
!   write(99,*)' '
!   write(99,*)'Eigenvalues '
!   write(99,*)' '
!   write(99,'(1000f15.10)')(wr(j),j=1,no)
!   write(99,*)' C (new basis)'
!   do i=1,no
!    write(99,'(1000f15.10)')(vr(i,j),j=1,no)
!   enddo
   
  
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
!   write(99,*)(wro(i),i=1,no-ns) 
   write(99,*)(wro(i),i=1,no) 
   write(99,*)

   
   
!   transformation of C to original basis

amatr=0.0d0
do i=1,ns
!vt(i,i)=1.0 !?????? preco
enddo

call dgemm('T','N',no,no,no,1.d0,vt,no,vr,no,0.d0,amatr,no)   !? spravne
vr=amatr
!!!!!!!!!!!!!
! reduce rows of D
deallocate(hamd)
allocate(hamd(no,ndim))

ii=0
do i=1,ndimr
 if (nx(i).ne.0) then
   ii=ii+1
!  do j=1,dim_base
      hamd(ii,:)=d1(i,:)
!  enddo
 endif
enddo

!write(99,*)' '
!write(99,*)' D matrix reduced rows'
!do i=1,dim_ind
! write(99,'(1000f15.10)')(hamd(i,j),j=1,dim_base)
!enddo

write(*,*)' check dim_ind =',ii
write(*,*)'X calculation'

allocate(xr(ndim,no))
!  X=DC

call dgemm('T','N',ndim,no,no,1.0d0,hamd,no,vr,no,0.d0,xr,ndim)

!write(998,*)' '
!write(998,*)' X matrix'
!do i=1,ndim
! write(998,'(1000f15.10)')(xr(i,j),j=1,no)
!enddo


call normalize_c(no,ndim,ndimr,vr,xr,mxtr,nx,jcal)


! test  C^T (AD) C
deallocate(amatr)
allocate(amatr(no,no))
amatr=0.0d0

call dgemm('N','N',no,no,no,1.d0,h_orig,no,vr,no,0.d0,amatr,no)
h_orig=0.d0
call dgemm('T','N',no,no,no,1.d0,vr,no,amatr,no,0.d0,h_orig,no)

iout=0
if (iout.eq.1) then
 write(998,*)
 write(998,*)'******** check  C^T (AD) C *************'
 write(998,*)
 do i=1,no
   write(998,'(1000f11.6)')(h_orig(i,j),j=1,no)
 enddo
endif  

!  copy of <spur| H | phys > correction 

allocate(h_corr(no,ns))

do i=1,no
  do j=1,ns
    h_corr(i,j)=h_orig(i,j)
  enddo
enddo





amatr=0.d0

call dgemm('N','N',no,no,no,1.d0,d_orig,no,vr,no,0.d0,amatr,no)
h_orig=0.d0
call dgemm('T','N',no,no,no,1.d0,vr,no,amatr,no,0.d0,d_orig,no)

iout=0
if (iout.eq.1) then
 write(998,*)
 write(998,*)'******** check  C^T (D) C *************'
 write(998,*)
 do i=1,no
   write(998,'(1000f11.6)')(d_orig(i,j),j=1,no)
 enddo
endif      


deallocate(h_orig,d_orig,amatr) 


!write(998,*)' '
!write(998,*)' C matrix normalized'
!do i=1,no
! write(998,'(1000f15.10)')(vr(i,j),j=1,no)
!enddo

!ns=0

!write(998,*) 'energies '
!write(998,'(1000f15.10)')(wr(i),i=1,no-ns)
!write(998,*)' '
!write(998,*)' X matrix normalized'
!do i=1,ndim
! write(998,'(1000f15.10)')(xr(i,j),j=1,no-ns)
!enddo

!allocate(ind_red_ind(dim_ind))

!ii=0
!do i=1,dim_baser
!if (nx(i).eq.1) then
!ii=ii+1
!ind_red_ind(ii)=ind_red(i)
!endif
!enddo


return            
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

!c      allocate(amatrp(ndimr,ndimr))
!c      amatrp=0.d0

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

!      include 'types_eqm.inc'

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
!c       if (e1.lt.10.d0.or.e2.lt.10.d0) then 

!c       if (i.eq.2*int(i/2)) then 

!c       if (il.ne.ispu.and.ilp.ne.ispu) then 
                                 irow(i)=1
                                 iii=iii+1
                         endif

      enddo

      ndimr=iii
     
      endsubroutine kicksp

!******************************************************************************      
      subroutine normal(no,ndmx,vr,xc,jcal,mxtr,mxt)

      implicit double precision (a-h,o-z)

!c      include 'chole.inc'

      double precision, dimension(:,:), allocatable :: vr,xc,d1
      integer, dimension (:), allocatable :: mxtr,mxt

      xfact=(dfloat(2*jcal+1))**0.5d0
!c      xc=0.d0

!c      do i=1,no
!c        do j=1,ndmx
!c          xpom=0.d0
!c           do k=1,no
!c             xpom=xpom+d1(k,j)*vr(k,i)
!c           enddo
!c             xc(j,i)=xpom
!c         enddo
!c      enddo

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
!        write(771,*)' i =',i,xpom
      enddo

      return
      end subroutine normal
!**************************************************************
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

!**********************************************************************
subroutine reduce_mat(nx,mat,ndim,ndim_red)
  implicit none
  double precision, dimension (:,:),allocatable :: mat,matc
  integer, dimension(:), allocatable :: nx
  integer i,j,ii,jj,ndim,ndim_red
  
  allocate(matc(ndim,ndim))
  matc=mat
  deallocate(mat)
  allocate(mat(ndim_red,ndim_red))
  
  
  ii=0
  do i=1,ndim
   if (nx(i).ne.0) then
     ii=ii+1
     jj=0
     do j=1,ndim
      if (nx(j).ne.0) then
         jj=jj+1
         mat(ii,jj)=matc(i,j)
       endif
      enddo
  
    endif
  enddo
  deallocate(matc)
  
  end subroutine reduce_mat
  
!***************************************************************************
subroutine normalize_c(idim_ind,idim_base,idim_baser,vr,xamp,ind_red,nx,ijj)

  implicit double precision (a-h,o-z)

!c      include 'chole.inc'

  double precision, dimension(:,:), allocatable :: vr,xamp,xampr
  integer, dimension (:), allocatable :: ind_red,nx

  allocate(xampr(idim_ind,idim_ind))

  do j=1,idim_ind
  ii=0
  do i=1,idim_baser
    if (nx(i).eq.1) then
      ii=ii+1
      xampr(ii,j)=xamp(ind_red(i),j)
    endif
  enddo
 enddo



  xfact=1.d0*(dfloat(2*ijj+1))**0.5d0

  do i=1,idim_ind
    xpom=0.d0
    do j=1,idim_ind
      xpom=xpom+vr(j,i)*xampr(j,i)
     enddo
    do j=1,idim_ind
      vr(j,i)=vr(j,i)/dsqrt(xpom)
    enddo
    do j=1,idim_base
      xamp(j,i)=xamp(j,i)/dsqrt(xpom)
!          xampr(j,i)=xampr(j,i)/dsqrt(xpom)
    enddo
!    write(771,*)' norm i =',i,xpom
  enddo

 xamp=xamp*xfact
 do j=1,idim_ind
  ii=0
  do i=1,idim_baser
    if (nx(i).eq.1) then
      ii=ii+1
      xampr(ii,j)=xamp(ind_red(i),j)
    endif
  enddo
 enddo



!     xampr=xampr*xfact
!    test of normalization

do k=1,idim_ind
do i=1,idim_ind
 xpom=0.0d0
 do j=1,idim_ind
   xpom=xpom+vr(j,k)*xampr(j,i)
 enddo
! if (dabs(xpom).gt.1.d-10) write(771,*) 'ovrl ',k,i,xpom,xpom/xfact
enddo
enddo

  return
end subroutine normalize_c
!**************************************************************


      

      end module hami
