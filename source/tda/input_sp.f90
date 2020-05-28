
 module input_sp


  contains 

  subroutine inp_sp(levn,levp,jmax)

      implicit none

      include 'types_tda_cp.inc'

      integer :: i,ii,jmax
      integer, dimension(:), allocatable:: il_n,il_p   
      type(level_typ),dimension(:), allocatable :: levn,levp     

      allocate(il_n(0:1000))
      il_n=0

      allocate(il_p(0:1000))
      il_p=0

      jmax=0
   
      open(1,file='singpart_coup.dat',status='old',form='formatted')

      read(1,*)
!      do i=1,ipnmx
!       read(1,16)nt,lt,jt,erspnt,ersppt
       i=0
       do while (.not.eof(1))
        i=i+1
        read(1,'(3i5,7x,f15.5,3i5,7x,f15.5)')levn(i)%n,levn(i)%l,levn(i)%j,levn(i)%en,levp(i)%n,levp(i)%l,levp(i)%j,levp(i)%en

        if (levn(i)%j.gt.jmax) jmax=levn(i)%j
        if (levp(i)%j.gt.jmax) jmax=levp(i)%j
 
        levp(i)%nr=il_p(levp(i)%l)
        levn(i)%nr=il_n(levn(i)%l)

        levp(i)%n=2*levp(i)%nr+levp(i)%l
        levn(i)%n=2*levn(i)%nr+levn(i)%l


        if (levn(i)%j.eq.2*levn(i)%l-1.or.levn(i)%l.eq.0) il_n(levn(i)%l)=il_n(levn(i)%l)+1
        if (levp(i)%j.eq.2*levp(i)%l-1.or.levp(i)%l.eq.0) il_p(levp(i)%l)=il_p(levp(i)%l)+1

 
        levp(i)%n=2*levp(i)%nr+levp(i)%l
        levn(i)%n=2*levn(i)%nr+levn(i)%l

!       levn(i)%n=nt
!       levn(i)%l=lt
!       levn(i)%j=jt
!       levn(i)%nr=(nt-lt)/2
!       levn(i)%en=alfa*erspnt
!       if (ikeyhf.eq.1) levn(i)%en=levn(i)%en+(dfloat(nt)+3.d0/2.d0)*homcm
!      write(96,*)i,orbit((nt-lt)/2),lt,jt
!       write(96,'(i5,i5,8x,i3,a3,i3)')i,nt,(nt-lt)/2,orbit(lt),jt
      enddo
      write(*,*)' Number of levels loaded ',i
      write(*,*)' Maximal value 2J',jmax

      close(1)

      open(1,file='singpart_coup.dat',status='unknown',form='formatted')

      write(1,*)'neutrons n l 2j                            protons  n l 2j'

      do ii=1,i
       write(1,'(3i5,7x,f15.5,3i5,7x,f15.5)')levn(ii)%n,levn(ii)%l,levn(ii)%j,levn(ii)%en,levp(ii)%n,levp(ii)%l,levp(ii)%j,levp(ii)%en
      enddo


      close(1)

   end subroutine inp_sp

    subroutine input_kin(kin_p,kin_n)
        implicit none

        include 'types_tda_cp.inc'

        integer :: id,i,j

        double precision, allocatable, dimension(:,:) :: kin_p, kin_n


        open(165,file='kin_nat_orb.dat',status='unknown',form='unformatted')
        read(165)id
        allocate(kin_p(id,id),kin_n(id,id))
        kin_n=0.d0
        kin_p=0.d0

        read(165)((kin_p(i,j),i=1,id),j=1,id)
        read(165)((kin_n(i,j),i=1,id),j=1,id)
        close(165)
   

    return
    end subroutine input_kin

end module input_sp



