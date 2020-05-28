!     last modification 7.1.2020  

      module phonon_base 

      use types_eqm  

      contains

      subroutine phonbase(nf,ipare,j,phonus,phonmus,idim1,idim2,idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr,n_spur,myid)

      implicit double precision (a-h,o-z)

    !  include 'types_eqm.inc'

      integer, dimension (:), allocatable :: phonus, phonmus
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      type(phonbase_typ), dimension (:), allocatable :: phonbs,phonbs_reor
      integer, dimension(:), allocatable :: mxtr,nxtr
      integer :: myid
      character*30 name1f,name2f


      if (nf.eq.2) name1f='1phonon/1f_states.dat'
      if (nf.eq.2) name2f='1phonon/1f_states.dat'

      if (nf.eq.3) name1f='1phonon/1f_states.dat'
      if (nf.eq.3) name2f='2phonon/2f_states.dat'
      
      if (allocated(phon1).eq..TRUE.) then 
       deallocate(phon1,phon2,phonus,phonmus)
       deallocate(phonbs)
      endif


      allocate(phon1(0:idim1),phon2(0:idim2),phonus(0:idim1),phonmus(0:idim2))
      allocate(phonbs(idimbs))


!   ALL PHONONS     
      do i=1,idim1
        phon1(i)%us=1
      enddo

      do i=1,idim2
        phon2(i)%us=1
      enddo

      open(881,file='tda_r_overl.dat',status='unknown',form='formatted')
      do while (.not.eof(881))
       read(881,*)iii,eee,rr
      enddo
     close(881)
 
 !     if (iii.eq.0) ilamcm=5000
       ilamcm=iii
       if (myid.eq.0) write(*,*)' CM phonon is number ',ilamcm

      do i=1,idim1
      phonus(i)=i
      enddo

      do i=1,idim2
      phonmus(i)=i
      enddo

      open (3,file=name1f,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       phon1(i)%par=ipar
       phon1(i)%j=ijj
       phon1(i)%enf=en
!       phon1(i)%us=1
      enddo

      nnphontr=i

      close(3)

      open (3,file=name2f,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       phon2(i)%par=ipar
       phon2(i)%j=ijj
       phon2(i)%enf=en
!       phon2(i)%us=1
      enddo

      nnphonmtr=i

      close(3)

      i=0

 

!      if (myid.eq.0) write(992,*)'  i  spur   1-ph: lamb    en_lam)    (n-1)-ph: gamma    en_gamma'

      do ilam=1,nnphonmtr

       do ilamp=1,nnphontr

        ilamu=phonmus(ilam)    ! used phonons
        iparl=phon2(ilamu)%par
        jl=phon2(ilamu)%j
                
          ilamup=phonus(ilamp)
          iparlp=phon1(ilamup)%par
          jlp=phon1(ilamup)%j

          eila=phon2(ilamu)%enf
          eilap=phon1(ilamup)%enf

          ipart=iparl*iparlp

          jmin=abs(jl-jlp)
          jmax=jl+jlp

      if (ipart.eq.ipare) then 
      if (j.le.jmax.and.j.ge.jmin) then         
 
          i=i+1

          if (i.gt.idimbs) then 
                  write(*,*)'Small dimension of array phonbs'
                  stop
                endif
          
          phonbs(i)%ila=ilamu
          phonbs(i)%ilap=ilamup  
          phonbs(i)%spur=0
!          if (ilamup.eq.ilamcm.or.eila.gt.1000.0d0) phonbs(i)%spur=1  ! how to define spurius states?
!          if (ilamup.eq.ilamcm) phonbs(i)%spur=1
          if (ilamu.eq.ilamcm.or.ilamup.eq.ilamcm) phonbs(i)%spur=1
      
      endif 
      endif 
 
        end do
        
      end do

      idphon=i

!      if (myid.eq.0) write(*,*)' Total dimension ',idphon

      !   reordering of basis
    if (1.eq.1) then
      allocate(phonbs_reor(idphon))
      ii=0
      do i=1,idphon
        if (phonbs(i)%spur.eq.1) then
          ii=ii+1
          phonbs_reor(ii)%ila=phonbs(i)%ila
          phonbs_reor(ii)%ilap=phonbs(i)%ilap
          phonbs_reor(ii)%spur=1
        endif
      enddo

      n_spur=ii  ! number of selected spurious states


      ii=0
      do i=1,idphon
        if (phonbs(i)%spur.eq.0) then
          ii=ii+1
          phonbs_reor(ii+n_spur)%ila=phonbs(i)%ila
          phonbs_reor(ii+n_spur)%ilap=phonbs(i)%ilap
          phonbs_reor(ii+n_spur)%spur=0
        endif
      enddo


      do i=1,idphon
        phonbs(i)%spur=phonbs_reor(i)%spur
        phonbs(i)%ila=phonbs_reor(i)%ila
        phonbs(i)%ilap=phonbs_reor(i)%ilap

!        if (myid.eq.0) write(992,'(i5,i5,5x,i5,3x,f15.5,5x,i5,3x,f15.5)')i,phonbs(i)%spur,phonbs(i)%ilap,phon1(phonbs(i)%ilap)%enf,phonbs(i)%ila,phon2(phonbs(i)%ila)%enf
      enddo

      deallocate(phonbs_reor)

    endif

     
      allocate(mxtr(idphon),nxtr(idphon))
      mxtr=0
      nxtr=0

!    write(*,*)' Energetic truncation? '
!     read(*,*)xetrun
!     write(*,*)xetrun
     open(27,file='en_trun.dat',status='old',form='formatted')
     read(27,*)xetrun
!     xetrun=100.0d0
!     write(*,*)' Energetic truncation of 2-ph subspace',xetrun
     close(27)



      jj=0

 
!      if (myid.eq.0) write(992,*)'reduced basis  '


      do ii=1,idphon
         ilamu=phonbs(ii)%ila
         ilamup=phonbs(ii)%ilap
         eila=phon2(ilamu)%enf
         eilap=phon1(ilamup)%enf
         ius1=phon1(ilamup)%us
         ius2=phon2(ilamu)%us

!         if ((phonbs(ii)%spur.eq.1).or.(phonbs(ii)%spur.eq.0.and.(eila+eilap).le.xetrun)) then   
!         if ((phonbs(ii)%spur.eq.1).or.(phonbs(ii)%spur.eq.0.and.eila.le.xetrun)) then 
         if ((phonbs(ii)%spur.eq.1.and.ilamu.ge.ilampup).or.(phonbs(ii)%spur.eq.0.and.eila.le.xetrun.and.ilamu.ge.ilamup)) then
!       if ((eila+eilap).le.xetrun.and.ilamu.le.ilamup.and.ius1.eq.1.and.ius2.eq.1.and.phonbs(ii)%spur.eq.0)  then 
!         if ((eila+eilap).le.xetrun.and.ilamu.le.ilamup.and.ius1.eq.1.and.ius2.eq.1)  then 

               jj=jj+1
               mxtr(jj)=ii
               nxtr(ii)=jj
!               if (myid.eq.0) write(992,'(i5,i5,5x,i5,3x,f15.5,5x,i5,3x,f15.5)')ii,phonbs(ii)%spur,ilamup,eilap,ilamu,eila
        endif


      enddo

      idphontr=jj

!      if (myid.eq.0) close(992)
!      if (myid.eq.0) write(*,*) ' Truncated dimension ',idphontr

!      open(2,file='mxtr.dat',status='unknown',form='unformatted')
!      write(2)idphontr
!      write(2)(mxtr(i),i=1,idphontr)
!      close(2)

!      open(2,file='nxtr.dat',status='unknown',form='unformatted')
!      write(2)idphon
!      write(2)(nxtr(i),i=1,idphon)
!      close(2)

      deallocate(nxtr)
 
      return

      end subroutine phonbase

      end module phonon_base
