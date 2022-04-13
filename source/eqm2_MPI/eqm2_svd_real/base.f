c     last modification 14.4.2010  

      module phonon_base 

      contains

      subroutine phonbase(nf,ipare,j,phonus,phonmus,idim1,idim2,
     *idimbs,idphon,idphontr,phonbs,phon1,phon2,mxtr)

      implicit double precision (a-h,o-z)

      include 'types_eqm.inc'

      integer, dimension (:), allocatable :: phonus, phonmus
      type(phon_typ), dimension (:), allocatable :: phon1,phon2
      type(phonbase_typ), dimension (:), allocatable :: phonbs
      integer, dimension(:), allocatable :: mxtr,nxtr

      character*30 name1f,name2f


      if (nf.eq.2) name1f='1phonon/1f_states.dat'
      if (nf.eq.2) name2f='1phonon/1f_states.dat'
      
      if (allocated(phon1).eq..TRUE.) then 
       deallocate(phon1,phon2,phonus,phonmus)
       deallocate(phonbs)
      endif


      allocate(phon1(0:idim1),phon2(0:idim2)
     *,phonus(0:idim1),phonmus(0:idim2))
      allocate(phonbs(idimbs))


      open(2,file='0hom_phon.dat',status='old',form='formatted')

!      phonus=0
!      phonmus=0

      do i=1,idim1
      phon1(i)%us=0
      enddo

      do i=1,idim2
      phon2(i)%us=0
      enddo

      do while (.not.eof(2))
      read(2,*)i,ee
!       write(*,*)i
       phon1(i)%us=1
       phon2(i)%us=1
      enddo

      close(2)


      open(2,file='1hom_phon.dat',status='old',form='formatted')


      do while (.not.eof(2))
      read(2,*)i,ee
!       write(*,*)i
       phon1(i)%us=1
       phon2(i)%us=1
      enddo

      close(2)

      open(2,file='2hom_phon.dat',status='old',form='formatted')

!      phonus=0
!      phonmus=0


      do while (.not.eof(2))
      read(2,*)i,ee
!       write(*,*)i
       if (ee.lt.15.0d0) then 
        phon1(i)%us=1
        phon2(i)%us=1
       endif
      enddo

      close(2)

      open(2,file='3hom_phon.dat',status='old',form='formatted')

!      phonus=0
!      phonmus=0


      do while (.not.eof(2))
      read(2,*)i,ee
!       write(*,*)i
       if (ee.lt.15.0d0) then
        phon1(i)%us=1
        phon2(i)%us=1
       endif
      enddo

      close(2)



!   ALL PHONONS     
      do i=1,idim1
      phon1(i)%us=0
      enddo

      do i=1,idim2
      phon2(i)%us=1
      enddo

      open (3,file=name1f,status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       if (en.le.100.0d0) phon1(i)%us=1
       phon1(i)%us=1
      enddo


      close(3)


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

      do ilamp=1,nnphontr

       do ilam=1,nnphonmtr

        ilamu=phonmus(ilam)    ! used phonons
        iparl=phon2(ilamu)%par
        jl=phon2(ilamu)%j
        
!        do ilamp=1,nnphontr

         
        
          ilamup=phonus(ilamp)
          iparlp=phon1(ilamup)%par
          jlp=phon1(ilamup)%j

          ipart=iparl*iparlp

          jmin=abs(jl-jlp)
          jmax=jl+jlp

          if (ipart.ne.ipare) goto 6
          if (j.gt.jmax.or.j.lt.jmin) goto 6
!          if (phon2(ilamu)%us.ne.1) goto 6
!          if (phon1(ilamup)%us.ne.1) goto 6          
 
          i=i+1

          if (i.gt.idimbs) then 
                  write(*,*)'Small dimension of array phonbs'
                  stop
                endif
          
          phonbs(i)%ila=ilamu
          phonbs(i)%ilap=ilamup  
          phonbs(i)%spur=0
      



          
  6       continue         
        end do
        
      end do

      idphon=i

      write(*,*)' Total dimension ',idphon

     
      allocate(mxtr(idphon),nxtr(idphon))
      mxtr=0
      nxtr=0

!      xetrun=25.0

!      write(*,*)' Energetic truncation? '
!      read(*,*)xetrun
!      write(*,*)xetrun
      xetrun=10000000000.0d0
!      xetrun=105.0

      jj=0
      do ii=1,idphon
         ilamu=phonbs(ii)%ila
         ilamup=phonbs(ii)%ilap
         eila=phon2(ilamu)%enf
         eilap=phon1(ilamup)%enf
         ius1=phon1(ilamup)%us
         ius2=phon2(ilamu)%us
         if (ius1.ne.0.and.ius2.ne.0) then
         write(991,*)ilamu,eila,ilamup,eilap,ius1,ius2
         endif

         if (eila.lt.0.d0) goto 661  ! surious states
         if (eilap.lt.0.d0) goto 661


       if ((eila+eilap).le.xetrun.and.ilamu.le.ilamup
     *.and.ius1.eq.1.and.ius2.ne.0.and.(ius1+ius2).le.3) then 
       
               jj=jj+1
               mxtr(jj)=ii
               nxtr(ii)=jj
         write(992,*)jj,ilamu,eila,ilamup,eilap
       endif

  661 continue       

      enddo

      idphontr=jj

      write(*,*) ' Truncated dimension ',idphontr

      open(2,file='mxtr.dat',status='unknown',form='unformatted')
      write(2)idphontr
      write(2)(mxtr(i),i=1,idphontr)
      close(2)

      open(2,file='nxtr.dat',status='unknown',form='unformatted')
      write(2)idphon
      write(2)(nxtr(i),i=1,idphon)
      close(2)

      deallocate(nxtr)
 
c      deallocate(phon1,phon2,phonus,phonmus,phonbs)

      return

      end subroutine phonbase

      end module phonon_base
