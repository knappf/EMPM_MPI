!     last modification 11.5.2018      
      
      program eqm 

      use types_eqm  
      use phonon_base
!      use metricmat
      use choleski
!      use admatr
      use read_admat
      use hami
      use cm_ort_svd

      implicit double precision (a-h,o-z)

!      include 'types_eqm.inc'

      type(phonbase_typ), dimension (:), allocatable :: phonbs
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      integer, dimension (:), allocatable :: phonus, phonmus
      type (amp2_typ), dimension(:,:), allocatable :: xcc

!     choleski arrays
      double precision, dimension(:,:), allocatable :: dd,cq,d1,cdu,xr,vr,h_corr
      double precision, dimension(:), allocatable :: wr
      
      integer, dimension(:), allocatable :: nx,mxt,irow,mxtr
      integer, dimension (:,:), allocatable ::ndbx
      integer, dimension (:,:,:), allocatable ::  ipozx
      integer, dimension (:), allocatable :: ndcc
      
      character*30 namex,names,namec



!      CALL OMP_SET_NUM_THREADS(8)

      nf=2  
      nlam=0

      if (nf.eq.2) then 
       namex='2phonon/2f_x.dat'
       namec='2phonon/2f_c.dat'
       names='2phonon/2f_states.dat'
      endif

!      call readx(namex2,xcc,ndcc,nd1f,nd2f)
!      call redrsumx(nd2f,nd1f,xcc,ndcc,ipozx,ndbx)

      idim1=10000
      idim2=1000000
      idimbs=50000000
      
      open(12,file=namex,status='unknown',form='unformatted')   ! file with X amplitudes 
      open(22,file=namec,status='unknown',form='unformatted')   ! file with C amplitudes 
      open(13,file=names,status='unknown',form='unformatted')
      open(99,file='2phon.log',status='unknown',form='formatted')
      open(993,file='h_corr.dat',status='unknown',form='formatted')

      
      call execute_command_line('rm 1_phon_dens_list.dat')



      open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')
      jmax=0
      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       if (ijj.gt.jmax) jmax=ijj
      enddo
      close(3)

!     J and Parity intervals       
      ipmin=-1
      ipmax=1
      jmin=0
      jmax=2*jmax

      write(*,*)'2-phonon calculation parity interval: ','<',ipmin,',',ipmax,'     >'
      write(*,*)'                         J  interval: ','<',jmin,',',jmax,'     >'  

      do ipar=ipmin,ipmax,2 
        do jcal=jmin,jmax,1
            write(*,*)
            write(*,*)'----------------------------------------------'
            write(*,*)' Parity = ',ipar,'   J = ',jcal
            call phonbase(nf,ipar,jcal,phonus,phonmus,idim1,idim2,idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr,n_spur)
        enddo
      enddo  

!      CALL execute_command_line('./run_dens2.sh' )           
!      write(*,*)'2-phon dens. calculated'
!      stop


      do ipar=ipmin,ipmax,2 !,-2        
      do jcal=jmin,jmax,1
     

  
      write(*,*)
      write(*,*)'----------------------------------------------'
      write(*,*)' Parity = ',ipar,'   J = ',jcal

 
      call phonbase(nf,ipar,jcal,phonus,phonmus,idim1,idim2,idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr,n_spur)
!     stop
      write(*,*)' Dimension = ',idphon
      write(*,*)' Truncated dimension = ',idphontr
      idphontot=idphon

      if (idphon.gt.0) then 
!      CALL execute_command_line('./run_dens2.sh' )      

! calculation of matrices A and D
!       call admat(nf,ipar,jcal,phonbs,idphon,no,phon1,phon2,mxt,xcc,ndcc,ipozx,ndbx)
       open(23,file='AD_J_Pi.dat',status='unknown',form='formatted')     
       write(23,*)ipar,jcal    
       close(23)
       CALL execute_command_line('./run_admat2.sh' )
 
!       stop

!  Choleski anal. of spurious subspace
       call read_sub_dmat_str(dd,n_spur,idphon)
!       call read_sub_dmat(dd,n_spur)
       write(*,*)
       write(*,*)' Cholesky analysis of spurious subspace'
!       call cholesk_dec(n_spur,dd,nx,mxt)

       call cholesk(1,n_spur,no,noo,dd,nx,mxt)

       open(3,file='chol_spur_subspace.dat',form='formatted',status='unknown')
       write(3,*)n_spur
       do k=1,n_spur
        write(3,*)nx(k)
       enddo
       close(3)

       deallocate(dd,nx,mxt)

!  Choleski anal. of D matrix
       write(*,*)
       write(*,*)' Cholesky analysis of D matrix '
       call cholesk(0,idphontr,no,idphontr,dd,nx,mxt)

      write(*,*)

!      stop

!      call ham(idphon,idphontr,no,nor,ns,irow,wr,xr,vr,ipar,jcal,mxtr,phonbs,nx,h_corr)
      call ham_geev(idphon,idphontr,no,nor,ns,irow,wr,xr,vr,ipar,jcal,mxtr,phonbs,nx,h_corr)

      write(*,*)
      write(*,*)' Number of CM spurious states ',ns
           
!      ns=0
      jc=no-ns+1
      
      do i=1,no
        write(13)nlam+i,ipar,jcal,wr(i)
!         do j=1,ns
!             write(993,*)nlam+i,nlam+j,h_corr(i,j)
!         enddo
         do j=jc,no
             write(993,*)nlam+i,nlam+j,h_corr(i,j)
         enddo


      enddo
     


!      do i=1,no
!       write(13)nlam+i,ipar,jcal,wr(i)
!      enddo

 
      write(*,*)' X ',idphontr,no,no-ns
      write(12)ipar,jcal,no,idphontot
      write(12)(phonbs(i)%ila,phonbs(i)%ilap,i=1,idphontot)  

!      write(22)ipar,jcal,no,no
!      write(22)(phonbs(mxt(i))%ila,phonbs(mxt(i))%ilap,i=1,no)

      mxt=0
      ii=0
      do i=1,idphontr
        if (nx(i).eq.1) then
            ii=ii+1
            mxt(ii)=mxtr(i)
        endif    
      enddo


     write(22)ipar,jcal,no,no
     write(22)(phonbs(mxt(i))%ila,phonbs(mxt(i))%ilap,i=1,no)

    !  write(122,*)(mxt(i),i=1,no)
    !  write(122,*)(mxtr(i),i=1,idphontr)
    !  write(122,*)(mxtr(mxt(i)),i=1,idphontr)


      write(*,*)' idphontot =',idphontot      
      do j=1,no
       write(12)(xr(i,j),i=1,idphontot)
       write(22)(vr(i,j),i=1,no)
      enddo

      nlam=nlam+no

!c      deallocate(xr,wr)
      
!c      deallocate(irow)

      deallocate(dd,nx,mxt,wr,xr,vr,h_corr)


      endif

      deallocate(mxtr)

!c      deallocate(phon1,phon2,phonus,phonmus,phonbs)

      write(*,*)'----------------------------------------------'

      CALL execute_command_line('rm ./scratch/a_mat*')
      CALL execute_command_line('rm ./scratch/d_mat*')
      enddo  ! cycle over J
      enddo  ! cycle over parity

      write(*,*) 'Number of 2-phonon states ', nlam
      
!c      close(66)
!c      close(33)
      close(99)
      close(12)
      close(22)
      close(892)
    

      end
