*     
*
*     last modification 5.10.2010   
*     
*      
      program fullham 

      use fullmat   

      implicit double precision (a-h,o-z)

      include 'types_fullham.inc'

      integer,  dimension (:), allocatable :: ndmpho
      integer,  dimension (:,:), allocatable :: nbs
      type (phon_typ), dimension (:,:), allocatable :: phonbs



      write(*,*)' Parity ? '
      read(*,*)iparc
c      iparc=1
      write(*,*)' J ? '
      read (*,*)jcalc
c      jcalc=0
      write(*,*)' Nf of full diagonalization?'
      read(*,*)nfon



      write(*,*)' Construction of total H matrix'
      write(*,*)' Parity = ',iparc
      write(*,*)' Angular momentum = ',jcalc

      do iff=0,nfon
        call bstot(iff,iparc,jcalc,nbs,ndmpho,phonbs)
      enddo

      write(*,*)'Dimension of 0 phonon space =  ',ndmpho(0)
      write(*,*)'Dimension of 1 phonon space =  ',ndmpho(1)
      write(*,*)'Dimension of 2 phonon space =  ',ndmpho(2)
      write(*,*)'Dimension of 3 phonon space =  ',ndmpho(3)

!      write(*,*)' Nf of full diagonalization?'
!      read(*,*)nfon

      if (nfon.eq.0) then 
                ndmpho(1)=0
                ndmpho(2)=0        
                ndmpho(3)=0
      endif

      if (nfon.eq.1) then 
              ndmpho(2)=0
              ndmpho(3)=0
      endif

      if (nfon.eq.2) then
          ndmpho(3)=0
      endif


!      if (ndmpho(1).gt.0) then 
!        write(*,*)'Calculating nondiagonal part 01'
!        call nondiag0(1,nbs,ndmpho,phonbs,iparc,jcalc)
!       endif 

      if (ndmpho(3).gt.0.and.ndmpho(1).gt.0) then 
       write(*,*)'Calculating nondiagonal part 13'
!      call nondiag2(3,nbs,ndmpho,phonbs,iparc,jcalc)
      endif
!      stop

      if (ndmpho(3).gt.0.and.ndmpho(2).gt.0) then
       write(*,*)'Calculating nondiagonal part 23'
!      call nondiag1r(3,nbs,ndmpho,phonbs,iparc,jcalc)
       endif

        
      if (ndmpho(2).gt.0.and.ndmpho(1).gt.0) then 
       write(*,*)'Calculating nondiagonal part 12'
!       call nondiag1r(2,nbs,ndmpho,phonbs,iparc,jcalc)
      endif
      
     
      if (ndmpho(2).gt.0.and.ndmpho(0).gt.0) then
       write(*,*)'Calculating nondiagonal part 02'
       call nondiag2(2,nbs,ndmpho,phonbs,iparc,jcalc)
      endif

      call fulhamt(ndmpho,nbs)

      end


