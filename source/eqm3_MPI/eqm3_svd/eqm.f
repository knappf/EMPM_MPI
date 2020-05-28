c     last modification 21.6.2010      
      
      program eqm 

      use phonon_base
!      use metricmat
      use choleski
      use admatr
      use hami

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      type(phonbase_typ), dimension (:), allocatable :: phonbs
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      integer, dimension (:), allocatable :: phonus, phonmus

c     choleski arrays
      double precision, dimension(:,:), allocatable ::
     *dd,cq,d1,cdu,xr,vr
      double precision, dimension(:), allocatable :: wr
      
      integer, dimension(:), allocatable :: nx,mxt,irow,mxtr
      
      character*30 namex,names,namec

c


!      CALL OMP_SET_NUM_THREADS(8)

      nf=2
      ip=1
      j=1
      idim1=10000
      idim2=10000
      idimbs=2000000
      
      
      
      nlam=0
      namex='2phonon/2f_x.dat'
      namec='2phonon/2f_c.dat'
      names='2phonon/2f_states.dat'
      
      open(12,file=namex,status='unknown',form='unformatted')
      open(22,file=namec,status='unknown',form='unformatted')
      open(13,file=names,status='unknown',form='unformatted')
      open(99,file='2phon.log',status='unknown',form='formatted')

      jmax=12
      jmin=0

      do ip=-1,1,2

      if (ip.eq.1) jcale=0
      if (ip.eq.-1) jcale=1

      


      do jcal=jmin,jmax
!       do jcal=jcale,jcale

  

      write(*,*)
      write(*,*)'----------------------------------------------'
      write(*,*)' Parity = ',ip,'   J = ',jcal

      call phonbase(nf,ip,jcal,phonus,phonmus,idim1,idim2,
     *idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr)

      write(*,*)' Dimension = ',idphon
      write(*,*)' Truncated dimension = ',idphontr
      idphontot=idphon

      if (idphon.gt.0) then 

!      if (jcal.ne.1) then 
!       call dmat(nf,ip,jcal,phonbs,idphon,phon1,phon2)

!       call dmat_test(idphontr)

!       stop

       call admat(nf,ip,jcal,phonbs,idphon,no,phon1,phon2,mxt)

!      endif

!        stop


!c       idphon=no

       call cholesk(idphontr,no,idphontr,dd,nx,mxt,D1,cdu)
       deallocate(dd,cdu,d1)

!c      write(*,*)' Number of linearly independent states ',no

!c      deallocate(dd,cdu,d1,nx)

!c      stop


!c     ispu=14
      
!c     call kicksp(ispu,irow,phonbs,phon1,phon2,mxt,no,nor)

      write(*,*)'OK'
     
      write(*,*)
  
      call ham(idphon,idphontr,no,nor,irow,wr,xr,vr,ip,jcal,mxtr)
 
      
      
      do i=1,no
       write(13)nlam+i,ip,jcal,wr(i)
      enddo

 
      write(*,*)' X ',idphontr,no
      write(12)ip,jcal,no,idphontot
      write(22)ip,jcal,no,no

      write(122,*)(mxt(i),i=1,no)
      write(122,*)(mxtr(i),i=1,idphontr)
      write(122,*)(mxtr(mxt(i)),i=1,idphontr)


      write(*,*)' idphontot =',idphontot      
      do j=1,no
c       write(12)(phonbs(mxtr(mxt(i)))%ila,phonbs(mxtr(mxt(i)))%ilap,
c     *xr(i,j),i=1,idphontr)


       write(12)(phonbs(i)%ila,phonbs(i)%ilap,
     *xr(i,j),i=1,idphontot)

       write(22)(phonbs(mxtr(mxt(i)))%ila,phonbs(mxtr(mxt(i)))%ilap,
     *vr(i,j),i=1,no)

       
c      write(122,*)(phonbs(mxtr(mxt(i)))%ila,phonbs(mxtr(mxt(i)))%ilap,
c     *i=1,idphontr)

      enddo

      nlam=nlam+no

      write(*,*)'*************'


c      deallocate(xr,wr)
      
c      deallocate(irow)

      deallocate(nx,mxt,wr,xr,vr)


      endif

      deallocate(mxtr)

c      deallocate(phon1,phon2,phonus,phonmus,phonbs)

      write(*,*)'----------------------------------------------'
      enddo  ! cycle over J
      enddo  ! cycle over parity

      write(*,*) ' Number of 2 phonon states ', nlam
      
c      close(66)
c      close(33)
      close(99)
      close(12)
      close(22)
    

      end
