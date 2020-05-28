module read_admat

    use types_eqm
    contains

    subroutine read_sub_dmat_str(dmatr,ndim,ndimt)

        double precision, dimension(:,:), allocatable :: dmatr
        double precision, dimension(:), allocatable :: dmm
        character (len=15) :: row_num  

        if (.not.allocated(dmatr))  allocate(dmatr(ndim,ndim))
        dmatr=0.d0
  
        allocate(dmm(ndimt))
        dmm=0.0d0

        do iii=1,ndim
          write(row_num,'(i15.15)')iii
          open(66,file='./scratch/d_mat_'//row_num,status='old',form='unformatted')
           read(66)(dmm(jjj),jjj=1,ndimt)
           dmatr(iii,1:ndim)=dmm(1:ndim)
          close(66) 
        enddo
  
        deallocate(dmm)

        return

    end subroutine read_sub_dmat_str

    subroutine read_dmat(dmatr)

        implicit double precision (a-h,o-z)
  
  !c      include 'formats_eqm.inc'
  
        double precision, dimension(:,:), allocatable :: dmatr
         double precision, dimension(:), allocatable :: dmm
  
        integer, dimension(:), allocatable :: nxtr
        character (len=15) :: row_number 
  !c      double precision, dimension(:), allocatable :: work,e
  
    
  
        open(2,file='nxtr.dat',status='unknown',form='unformatted')
        read(2)ndim_tot
        allocate(nxtr(ndim_tot))
        nxtr=0
        read(2)(nxtr(i),i=1,ndim_tot)
        close(2)
  
        write(991,*)'ndim_tot =', ndim_tot
        write(991,*)(nxtr(i),i=1,ndim_tot)

        ndim_tr=0
        do i=1,ndim_tot
            if (nxtr(i).ne.0) ndim_tr=ndim_tr+1
        enddo

!        allocate(dmatr(ndim_tr,ndim_tot))
        dmatr=0.d0
  
  !c      iold=-100
  
   
!        read(66)ndimrt,ndimt

        allocate(dmm(ndim_tot))
        dmm=0.0d0
  
  !      if (ndimrt.ne.ndimr.or.ndimt.ne.ndim) then
  !        write(*,*)' Dimensions does not match in d_m file'
  !        stop
  !      endif

        
        do i=1,ndim_tr
         write(row_number,'(i15.15)')i
         open(66,file='./scratch/d_mat_'//row_number,status='old',form='unformatted') 
         read(66)(dmm(j),j=1,ndim_tot)
         do j=1,ndim_tot
           if (nxtr(j).ne.0) then
            dmatr(i,nxtr(j))=dmm(j)
           endif
         enddo
         close(66)
        enddo
  
 
        deallocate(dmm)
  
  
  
  
  !      do while (.not.eof(66))
  !      read(66)i,j,dd
  !      if (nxtr(j).ne.0) then 
  !      dmatr(i,nxtr(j))=dd
  !      endif
  
  !      enddo
  
  !      close(66)
  
        deallocate(nxtr)
        
        end subroutine read_dmat

        subroutine read_admatr(path,matr,ndim_tr,ndim_tot)

            double precision, dimension(:,:), allocatable :: matr
            double precision, dimension(:), allocatable :: dmm
            character (len=15) :: row_num
            character (len=16) :: path   
    
            if (.not.allocated(matr))  allocate(matr(ndim_tr,ndim_tot))
            matr=0.d0
      
            allocate(dmm(ndim_tot))
            dmm=0.0d0
    
            do i=1,ndim_tr
              write(932,*)i,ndim_tr
              write(row_num,'(i15.15)')i
!              open(66,file='./scratch/d_mat_'//row_num,status='old',form='unformatted')
              open(66,file=path//row_num,status='old',form='unformatted')
               read(66)(dmm(j),j=1,ndim_tot)
               matr(i,1:ndim_tot)=dmm(1:ndim_tot)
              close(66) 
            enddo
      
            deallocate(dmm)
    
            return

        end subroutine read_admatr

        subroutine read_admatr_OMP(path,matr,ndim_tr,ndim_tot)

            double precision, dimension(:,:), allocatable :: matr
            double precision, dimension(:), allocatable :: dmm
            character (len=15) :: row_num
            character (len=16) :: path   
    
            if (.not.allocated(matr))  allocate(matr(ndim_tr,ndim_tot))
            matr=0.d0
      
!            allocate(dmm(ndim_tot))
!            dmm=0.0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,row_num,ifile)
!$OMP DO 
            do i=1,ndim_tr
!              write(932,*)i,ndim_tr
              write(row_num,'(i15.15)')i
              ifile=8000+i
              open(ifile,file=path//row_num,status='old',form='unformatted')
!               read(ifile)(dmm(j),j=1,ndim_tot)
               read(ifile)(matr(i,j),j=1,ndim_tot)
!               matr(i,1:ndim_tot)=dmm(1:ndim_tot)
              close(ifile) 
            enddo
!$OMP END DO
!$OMP END PARALLEL      
      
!            deallocate(dmm)
    
            return

        end subroutine read_admatr_OMP

    subroutine readx(fname,xcc,ndcc,nd1f,nd2f)

        implicit double precision (a-h,o-z)
  
  !      include 'types_eqm.inc'
  
        type (amp2_typ), dimension(:,:), allocatable :: xcc
        integer, dimension (:), allocatable :: ndcc
  
        character(len=30)fname
  
  !c      allocate(xcc(ndimi,ndimj))
  !c      allocate(ndcc(ndimi))
  
        open (3,file='1phonon/1f_states.dat',status='old',form='unformatted')
  
        do while (.not.eof(3))
         read(3)i,ipart,ijj,en
        enddo
        close(3)
  
        nd1f=i
        write(*,*)'Number of 1-phonon states', nd1f
  
  
        open (3,file='2phonon/2f_states.dat',status='old',form='unformatted')
  
        do while (.not.eof(3))
         read(3)i,ipart,ijj,en
        enddo
        close(3)
  
        nd2f=i
        write(*,*)'Number of 2-phonon states', nd2f
  
  
        open(2,file=fname,status='old',form='unformatted')
  
        ilamp=0
  
        ndimj=40000
  
        allocate(xcc(nd2f,ndimj))
        allocate(ndcc(nd2f))
  
  
  
        do while (.not.eof(2))
  
        read(2)ipar,ijj,no,idphon
  
  !c      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)
  
  
        do ilam=1,no
  
           if (idphon.gt.ndimj) then
        write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
                 stop
             endif
  
         read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is,xcc(ilam+ilamp,i)%am,i=1,idphon)
         ndcc(ilam+ilamp)=idphon
        enddo
  
        ilamps=ilamp
        ilamp=ilamp+no
  
  !      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11
  
  !      deallocate(xcc,ndcc)
        no=0
        idphon=0
  
        enddo
  
    11  continue
  
        close(2)
  
        return
        end subroutine readx
  !***********************************************************
        subroutine redrsumx(iamax,ifonmx,xx,ndx,ipozbr,ndbr)
  
        implicit double precision (a-h,o-z)
  
        include 'formats_eqm.inc'
  !      include 'types_eqm.inc'
  
        type (amp2_typ), dimension(:,:), allocatable :: xx
  
        integer, dimension (:,:,:), allocatable :: ipozbr
        integer, dimension (:,:), allocatable :: ndbr
        integer, dimension (:), allocatable :: ndx
  
  
        character(len=30)fname
  
        ndipo=1000
  
        if (.not.allocated(ipozbr)) allocate(ipozbr(iamax,ifonmx,ndipo))
        ipozbr=0
  
        if (.not.allocated(ndbr)) allocate(ndbr(iamax,ifonmx))
        ndbr=0
  
        do ia=1,iamax
  
  !c      write(*,*)'ia =',ia,ndx(ia)
  
        do i=1,ndx(ia)
         ibt=xx(ia,i)%is
  !c       write(*,*)i,ibt,ndbr(ia,ibt)
         ndbr(ia,ibt)=ndbr(ia,ibt)+1
  
          if (ndbr(ia,ibt).gt.ndipo) then
                 write(*,*)' Increase dimension in redrsumx'
                 stop
              endif
         ipozbr(ia,ibt,ndbr(ia,ibt))=i
        enddo
  
        enddo
  
  
        return
        end subroutine redrsumx
  
  !***********************************************************************    

end module read_admat
