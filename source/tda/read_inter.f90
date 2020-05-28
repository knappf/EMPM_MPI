module read_inter
    implicit none
    contains


    subroutine read_v(file_int,nlevmaxp,nlevmaxn,jmax,vp,vn,vpn)
    implicit none

    include 'types_tda_cp.inc'

    double precision :: ee,qei
    double precision, dimension(:,:,:,:,:),allocatable :: vp,vn,vpn
    integer (kind=1) :: ip,ij,it
    integer (kind=2) :: ia,ib,ic,id
    integer :: iaa,ibb,icc,idd,itzt,jmin,jmax,nlevmax,jtotmax,ifazab,ifazcd,nlevmaxp,nlevmaxn,nlevdim
    integer :: i,ii,iii,ll,jj,k
    double precision :: vv
    double precision :: factab,factcd,xnorm

    type(level_typ),dimension(:),allocatable :: lev


    character(len=30)file_int

    nlevdim=10000
    allocate(lev(nlevdim))
    lev%n=0
 !   levn%n=0
 !   levp%n=0
    jmax=0

    iii=0
    open(1,file='NAT_p.out',status='unknown',form='formatted')
    read(1,*)
      do k=1,nlevmaxp
!         read(1,*) lhfp(ii)%index,lhfp(ii)%l,lhfp(ii)%j2,
!     &    lhfp(ii)%ei,lhfp(ii)%qei,lhfp(ii)%vi

       read(1,*)ii,ll,jj,ee,qei
       iii=iii+1

!        lev(i)%n=nn
      i=2*ii-1
      lev(i)%l=ll
      lev(i)%j=jj
      lev(i)%en=ee
      lev(i)%it=-1

!      levp(ii)%l=ll
!      levp(ii)%j=jj
!      levp(ii)%en=ee
       if (jj.gt.jmax) jmax=jj

      enddo
     close(1)

!     nlevmaxp=iii
     write(*,*)' Number of proton  levels',nlevmaxp

     iii=0
     open(1,file='NAT_n.out',status='unknown',form='formatted')
     read(1,*)
      do k=1,nlevmaxn
!         read(1,*) lhfp(ii)%index,lhfp(ii)%l,lhfp(ii)%j2,
!     &    lhfp(ii)%ei,lhfp(ii)%qei,lhfp(ii)%vi

       read(1,*)ii,ll,jj,ee,qei
       iii=iii+1

!        lev(i)%n=nn
      i=2*ii
      lev(i)%l=ll
      lev(i)%j=jj
      lev(i)%en=ee
      lev(i)%it=1

!      levn(ii)%l=ll
!      levn(ii)%j=jj
      !levn(ii)%en=ee

      if (jj.gt.jmax) jmax=jj

      end do
     close(1)

!    nlevmaxn=iii

    write(*,*)' Number of neutron levels',nlevmaxn





    open(2,status='unknown',file=file_int,form='unformatted')
 


    write(*,*)' Reading proton-proton interaction m.e. '
    write(*,*)
    jmin=0

    allocate(vp(0:jmax,1:nlevmaxp,1:nlevmaxp,1:nlevmaxp,1:nlevmaxp))
    vp=0.d0

    itzt=-1   ! Proton-proton interaction
    jtotmax=0
    do while (.not.eof(2))
      read(2)it,ij,ia,ib,ic,id,vv


      if (it.eq.itzt) then
         if (ij/2.le.jmax) then

         if (ij.gt.jtotmax) jtotmax=ij
         if (jtotmax/2.gt.jmax) then
                 write(*,*)' Jmax larger then in input'
                 stop
         endif

         ifazab=(-1)**((lev(ia)%j+lev(ib)%j-ij)/2)
         ifazcd=(-1)**((lev(ic)%j+lev(id)%j-ij)/2)

         iaa=int(ia/2)+1
         ibb=int(ib/2)+1
         icc=int(ic/2)+1
         idd=int(id/2)+1

         if (iaa.le.nlevmaxp.and.ibb.le.nlevmaxp.and.icc.le.nlevmaxp.and.idd.le.nlevmaxp) then

         factab=1.d0
         factcd=1.d0
         if (iaa.eq.ibb) factab=dsqrt(2.d0)
         if (icc.eq.idd) factcd=dsqrt(2.d0)
         xnorm=factab*factcd
         vv=vv*xnorm


         vp(ij/2,iaa,ibb,icc,idd)=vv
         vp(ij/2,icc,idd,iaa,ibb)=vv
         vp(ij/2,ibb,iaa,icc,idd)=dfloat(-1*ifazab)*vv
         vp(ij/2,icc,idd,ibb,iaa)=dfloat(-1*ifazab)*vv
         vp(ij/2,iaa,ibb,idd,icc)=dfloat(-1*ifazcd)*vv
         vp(ij/2,idd,icc,iaa,ibb)=dfloat(-1*ifazcd)*vv
         vp(ij/2,ibb,iaa,idd,icc)=dfloat(ifazcd*ifazab)*vv
         vp(ij/2,idd,icc,ibb,iaa)=dfloat(ifazcd*ifazab)*vv

         endif
        endif
       endif
    enddo

    jtotmax=jtotmax/2

    write(*,*)'J_max of loaded PP interaction m.e ',jtotmax
    write(*,*)
     
    rewind(2)

    jmin=0

    write(*,*)' Reading neutron-neutron interaction m.e. '
    allocate(vn(0:jmax,1:nlevmaxn,1:nlevmaxn,1:nlevmaxn,1:nlevmaxn))
    vn=0.d0

    itzt=1   ! Neutron-neutron interaction
    jtotmax=0
    do while (.not.eof(2))

!        read(2,100)it,ip,ij,ia,ib,ic,id,vpl(1)!,vpl(2),vpl(3)
!     *,vpl(4),vpl(5),vpl(6),vpl(7),vpl(8)

      read(2)it,ij,ia,ib,ic,id,vv

!        vv=vpl(ius)


      if (it.eq.itzt) then
        if (ij/2.le.jmax) then

         if (ij.gt.jtotmax) jtotmax=ij
         if (jtotmax/2.gt.jmax) then
                 write(*,*)' Jmax larger then in input'
                 stop
         endif

         ifazab=(-1)**((lev(ia)%j+lev(ib)%j-ij)/2)
         ifazcd=(-1)**((lev(ic)%j+lev(id)%j-ij)/2)

         iaa=int(ia/2)
         ibb=int(ib/2)
         icc=int(ic/2)
         idd=int(id/2)

         if (iaa.le.nlevmaxn.and.ibb.le.nlevmaxn.and.icc.le.nlevmaxn.and.idd.le.nlevmaxn) then

         factab=1.d0
         factcd=1.d0
         if (iaa.eq.ibb) factab=dsqrt(2.d0)
         if (icc.eq.idd) factcd=dsqrt(2.d0)
         xnorm=factab*factcd
         vv=vv*xnorm


         vn(ij/2,iaa,ibb,icc,idd)=vv
         vn(ij/2,icc,idd,iaa,ibb)=vv
         vn(ij/2,ibb,iaa,icc,idd)=dfloat(-1*ifazab)*vv
         vn(ij/2,icc,idd,ibb,iaa)=dfloat(-1*ifazab)*vv
         vn(ij/2,iaa,ibb,idd,icc)=dfloat(-1*ifazcd)*vv
         vn(ij/2,idd,icc,iaa,ibb)=dfloat(-1*ifazcd)*vv
         vn(ij/2,ibb,iaa,idd,icc)=dfloat(ifazcd*ifazab)*vv
         vn(ij/2,idd,icc,ibb,iaa)=dfloat(ifazcd*ifazab)*vv

          endif
        endif
       endif
    enddo

    jtotmax=jtotmax/2

    write(*,*)'J_max of loaded NN interaction m.e ',jtotmax
    write(*,*)


    write(*,*)' Reading m.e. of proton-neutron interaction '
    write(*,*)
    rewind(2)
    jmin=0


    allocate(vpn(0:jmax,1:nlevmaxp,1:nlevmaxn,1:nlevmaxp,1:nlevmaxn))
    vpn=0.d0

    itzt=0   ! Proton-neutron interaction
    jtotmax=0
    do while (.not.eof(2))

!        read(2,100)it,ip,ij,ia,ib,ic,id,vpl(1)!,vpl(2),vpl(3)
!     *,vpl(4),vpl(5),vpl(6),vpl(7),vpl(8)
!        write(222,100)it,ip,ij,ia,ib,ic,id,vpl(1)

      read(2)it,ij,ia,ib,ic,id,vv

!        vv=vpl(ius)


      if (it.eq.itzt) then
       if (ij/2.le.jmax) then
         if (ij.gt.jtotmax) jtotmax=ij
         if (jtotmax/2.gt.jmax) then
                 write(*,*)' Jmax larger then in input'
                 stop
         endif


         iaa=int(ia/2)+1
         ibb=int(ib/2)
         icc=int(ic/2)+1
         idd=int(id/2)

         if (iaa.le.nlevmaxp.and.ibb.le.nlevmaxn.and.icc.le.nlevmaxp.and.idd.le.nlevmaxn) then

         vpn(ij/2,iaa,ibb,icc,idd)=vv
         vpn(ij/2,icc,idd,iaa,ibb)=vv

         endif
        endif
       endif
    enddo

    jtotmax=jtotmax/2

    write(*,*)' J_max of loaded interaction PN  m.e ',jtotmax
    write(*,*)

    close(2)

    
    
    end subroutine read_v

    subroutine read_f(file_fp,file_fn,file_fpn,nlevmaxp,nlevmaxn,jmax,fp,fn,fpn)
    implicit none
    
    double precision, dimension(:,:,:,:,:), allocatable :: fp,fn,fpn

    integer :: nlevmaxp,nlevmaxn,jmax
    character(len=30)file_fp,file_fn,file_fpn
    
    double precision :: vint

    integer(kind=1) :: j_f
    integer(kind=2) :: i_i,i_j,i_k,i_l
    integer :: i,j,k,l,ijt


    allocate(fp(0:jmax,nlevmaxp,nlevmaxp,nlevmaxp,nlevmaxp))
    allocate(fn(0:jmax,nlevmaxn,nlevmaxn,nlevmaxn,nlevmaxn))
    allocate(fpn(0:jmax,nlevmaxp,nlevmaxp,nlevmaxn,nlevmaxn))

    fp=0.d0
    fn=0.d0
    fpn=0.d0

    open(1,file=file_fp,status='old',form='unformatted')
   
    do while (.not.eof(1))

      read(1)j_f,i_i,i_j,i_k,i_l,vint
       ijt=j_f
       i=i_i
       j=i_j
       k=i_k
       l=i_l

    
    if (i.le.nlevmaxp.or.k.le.nlevmaxp.or.j.le.nlevmaxp.or.l.le.nlevmaxp) fp(ijt,i,j,k,l)=vint

    enddo
    close(1)

    open(1,file=file_fn,status='old',form='unformatted')

    do while (.not.eof(1))
      read(1)j_f,i_i,i_j,i_k,i_l,vint
       ijt=j_f
       i=i_i
       j=i_j
       k=i_k
       l=i_l

    if (i.le.nlevmaxn.or.k.le.nlevmaxn.or.j.le.nlevmaxn.or.l.le.nlevmaxn) fn(ijt,i,j,k,l)=vint
   
    enddo
    close(1)


      
   open(1,file=file_fpn,status='old',form='unformatted')
      
       
   do while (.not.eof(1))
      read(1)j_f,i_i,i_j,i_k,i_l,vint
         ijt=j_f
         i=i_i
         j=i_j
         k=i_k
         l=i_l
      
         if (i.le.nlevmaxp.or.j.le.nlevmaxp.or.k.le.nlevmaxn.or.l.le.nlevmaxn) fpn(ijt,i,j,k,l)=vint
      
   enddo
      
   close(1)
            

    end subroutine read_f



end module read_inter