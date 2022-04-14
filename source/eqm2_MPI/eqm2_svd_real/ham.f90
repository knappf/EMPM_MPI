!c     last modification 11.5.2018

module hami

use types_eqm
use read_admat

contains

subroutine ham_geev(ndim,ndimr,no,nor,ns,irow,wr,xr,vr,ipar,jcal,mxtr,phonbs,nx,h_corr)

  use choleski
  use cm_ort_svd

  implicit double precision (a-h,o-z)
   
!      include 'types_eqm.inc'
  include 'formats_eqm.inc'

  type(phonbase_typ), dimension (:), allocatable :: phonbs
  real(kind=4), dimension(:,:), allocatable :: d1,amatr,hami,hamid,dmatr,cq,vr,d1r,hamir,hamidr,hamd,xr
  real(kind=4), dimension(:,:), allocatable :: d_orig,h_orig,h_corr
 
  real(kind=8), dimension(:), allocatable ::  work,wr,wi,wro
  real(kind=8), dimension(:,:), allocatable :: hami_dbl, dmatr_dbl

  integer, dimension(:), allocatable :: mxt,nxt,mxtr,ipoz,irow,nx,ifail

  real(kind=4), dimension (:), allocatable :: s
  real(kind=4), dimension (:,:), allocatable :: vt
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
  d1=0.0

  allocate(amatr(ndimr,ndim))
  amatr=0.0

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

!  call read_admatr('./scratch/d_mat_',d1,ndimr,ndim)
  time_est = SECNDS(0.0)
  write(*,*)' Reading D'
  call read_admatr_OMP('./scratch/d_mat_',d1,ndimr,ndim)
  time_run = SECNDS(time_est)
  print '("Loading time  = ",f15.5," s.")',time_run

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

!  call read_admatr('./scratch/a_mat_',amatr,ndimr,ndim)
  time_est = SECNDS(0.0)
  write(*,*)' Reading A' 
  call read_admatr_OMP('./scratch/a_mat_',amatr,ndimr,ndim)
  time_run = SECNDS(time_est)
  print '("Loading time  = ",f15.5," s.")',time_run

  write(*,*)'A,D dimensions :', ndim,ndimr

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
  hami=0.0

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
  
  write(*,*)' AD product'
  call sgemm('N','T',ndimr,ndimr,ndim,1.0,amatr,ndimr,d1,ndimr,0.0,hami,ndimr)

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
      if (abs(hami(i,j)-hami(j,i)).gt.0.01) then
       write(*,*)'Non-symmetric AD'
       write(*,'(2i5,3f10.5)')i,j,hami(i,j),hami(j,i),hami(i,j)-hami(j,i)
      endif   
 
   enddo
  enddo

!     SVD part 
  call dmat_spur_set(ndim,ndimr,no,phonbs,d1,mxtr,nx,s,vt,ns)

allocate(dmatr(ndimr,ndimr))
dmatr=0.0 

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
 amatr=0.0

 call sgemm('N','T',no,no,no,1.0,hami,no,vt,no,0.0,amatr,no)
 hami=0.0
 call sgemm('N','N',no,no,no,1.0,vt,no,amatr,no,0.0,hami,no)


!do i=1,no
!do j=i,no
!      if (dabs(hami(i,j)-hami(j,i)).gt.0.01d0) write(*,'(2i5,3f10.5)'),i,j,dabs(hami(i,j)-hami(j,i))
!enddo
!enddo

! transformation of reduced D
!   allocate(dmatc_orig(no,no))
!   dmatc_orig=d1

! deallocate(amat)
!  allocate(amat(dim_ind,dim_ind))

 call sgemm('N','T',no,no,no,1.0,dmatr,no,vt,no,0.0,amatr,no)
 dmatr=0.0
 call sgemm('N','N',no,no,no,1.0,vt,no,amatr,no,0.0,dmatr,no)



! decoupling spurious subspace     

 decoup_cm=.true.
 if (decoup_cm.eq..TRUE.) then 
        do i=1,ns
      hami(i,i)=hami(i,i)+1000000000.0
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
       write(998,102)(dmatr(i,j),j=1,no)
     enddo
    endif

    do i=1,no
      do j=i,no
    
        if (abs(dmatr(i,j)-dmatr(j,i)).gt.0.01) then 
        write(*,*)' Transformed D asymmetric!'
        write(*,'(2i5,3f10.5)')i,j,abs(dmatr(i,j)-dmatr(j,i))
        endif
      enddo
     enddo


   iout=0
    if (iout.eq.1) then
     write(998,*)
     write(998,*)'******** transfromed matrix  AD *************'
     write(998,*)
     do i=1,no
       write(998,102)(hami(i,j),j=1,no)
     enddo
    endif      

    do i=1,no
      do j=i,no
    
        if (abs(hami(i,j)-hami(j,i)).gt.0.01) then 
        write(*,*)' Transformed H asymmetric!'
        write(*,'(2i5,3f10.5)')i,j,abs(hami(i,j)-hami(j,i))
        endif
      enddo
     enddo    


 !  Generalized EGV problem
   
  
    write(*,*)' Diagonalisation '
    lwork=20*no
    allocate(work(lwork),ifail(no),wr(no),wro(no),hami_dbl(no,no),dmatr_dbl(no,no))
    wr=0.0d0
    wro=0.0d0
    ifail=0
    work=0.0d0
  
    hami_dbl=0.0d0
    dmatr_dbl=0.d0

!    do i=1,no
!      do j=1,no
        hami_dbl=dble(hami)
        dmatr_dbl=dble(dmatr)
!      enddo
!    enddo
     deallocate(dmatr)

!    call dsygvx(1,'V', 'V', 'U', no, hami, no, dmatr, no, vl, vu, il, iu, abstol, m, wr, vr, no, work, lwork, lwork, ifail, info)

!    call ssygv(1, 'V','U', no , hami, no, dmatr, no, wr, work, lwork, info)

    call dsygv(1, 'V','U', no , hami_dbl, no, dmatr_dbl, no, wr, work, lwork, info)

!    write(*,*)' Number of eigevalues is interval (',vl,vu,')  :',m
     write(*,*)'Diagonalization completed!'


write(99,*)' Parity = ',ipar, ' J = ',jcal
write(99,*)
!   write(99,*)(wro(i),i=1,no-ns) 
write(99,*)(wr(i),i=1,no) 
write(99,*)

hami=real(hami_dbl)
deallocate(hami_dbl,dmatr_dbl)

!   transformation of C to original basis

amatr=0.0
!do i=1,ns
!vt(i,i)=1.0 !?????? preco
!enddo

allocate(vr(no,no))
vr=0.0

call sgemm('T','N',no,no,no,1.0,vt,no,hami,no,0.0,amatr,no)   !? spravne
vr=amatr
!!!!!!!!!!!!!
! reduce rows of D
!deallocate(hamd)
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

deallocate(amatr,d1,hami)

write(*,*)' check dim_ind =',ii
write(*,*)'X calculation'

allocate(xr(ndim,no))
xr=0.0
!  X=DC
call sgemm('T','N',ndim,no,no,1.0,hamd,no,vr,no,0.0,xr,ndim)
write(*,*)'C normalization'

iout=0
if (iout.eq.1) then
write(998,*)' '
write(998,*)' X matrix'
do i=1,ndim
 write(998,102)(xr(i,j),j=1,no)
enddo
endif 


call normalize_c(no,ndim,ndimr,vr,xr,mxtr,nx,jcal)


! test  C^T (AD) C
!deallocate(amatr)
allocate(amatr(no,no))
amatr=0.0

call sgemm('N','N',no,no,no,1.0,h_orig,no,vr,no,0.0,amatr,no)
h_orig=0.0
call sgemm('T','N',no,no,no,1.0,vr,no,amatr,no,0.0,h_orig,no)

iout=0
if (iout.eq.1) then
write(998,*)
write(998,*)'******** check  C^T (AD) C *************'
write(998,*)
do i=1,no
write(998,102)(h_orig(i,j),j=1,no)
enddo
endif  

!  copy of <spur| H | phys > correction 

!allocate(h_corr(no,ns))

!do i=1,no
!do j=1,ns
!h_corr(i,j)=h_orig(i,j)
!enddo
!enddo

!with the new diag routine the spurious states are the last of h and not the
!first as before

jc=no-ns+1
allocate(h_corr(no,jc:no))

do i=1,no
  do j=jc,no!1,ns
    h_corr(i,j)=h_orig(i,j)
  enddo
enddo





amatr=0.0

call sgemm('N','N',no,no,no,1.0,d_orig,no,vr,no,0.0,amatr,no)
h_orig=0.0
call sgemm('T','N',no,no,no,1.0,vr,no,amatr,no,0.0,d_orig,no)

iout=0
if (iout.eq.1) then
write(998,*)
write(998,*)'******** check  C^T (D) C *************'
write(998,*)
do i=1,no
write(998,102)(d_orig(i,j),j=1,no)
enddo
endif      


deallocate(h_orig,d_orig,amatr) 


!write(998,*)' '
!write(998,*)' C matrix normalized'
!do i=1,no
! write(998,'(1000f15.10)')(vr(i,j),j=1,no)
!enddo

!ns=0

! control output of X
iout=0
if (iout.eq.1) then
 write(998,*) 'energies '
 write(998,102)(wr(i),i=1,no-ns)
 write(998,*)' '
 write(998,*)' X matrix normalized'
do i=1,ndim
 write(998,102)(xr(i,j),j=1,no-ns)
enddo
endif

!allocate(ind_red_ind(dim_ind))

!ii=0
!do i=1,dim_baser
!if (nx(i).eq.1) then
!ii=ii+1
!ind_red_ind(ii)=ind_red(i)
!endif
!enddo


return            
end subroutine ham_geev



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
  real(kind=4), dimension (:,:),allocatable :: mat,matc
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

    implicit real(kind=4) (a-h,o-z)
  
  !c      include 'chole.inc'
  
    real(kind=4), dimension(:,:), allocatable :: vr,xamp,xampr
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
  
  
  
    xfact=1.0*(float(2*ijj+1))**0.5
  
    do i=1,idim_ind
      xpom=0.0
      do j=1,idim_ind
        xpom=xpom+vr(j,i)*xampr(j,i)
       enddo
      do j=1,idim_ind
        vr(j,i)=vr(j,i)/sqrt(xpom)
      enddo
      do j=1,idim_base
        xamp(j,i)=xamp(j,i)/sqrt(xpom)
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
   xpom=0.0
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
