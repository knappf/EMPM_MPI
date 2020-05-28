module inter_stor
!
use int_arr
use technical
USE geom

contains 
!
!------------------------------------------
subroutine read_inter(iflag,hbarom,iz,in,imax,jmax,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn,v_pn,irowc_pn,irowe_pn)
implicit none

integer :: imax,n_dim,i_poz,i_pozr,i_j,k_l,ndimv,ndimv2,ndimr,i_pp,in_pp,in_nn,in_pn,iz,in,iflag
integer (kind=1) :: ip,ij,it,ipart,jcp,jmax
integer (kind=2) :: i1,i2,i3,i4
integer, dimension (:,:), allocatable :: i_pair 
integer (kind=2), dimension (:), allocatable :: i_bs,j_bs
double precision :: hbarom
real(kind=8):: v_elem,v_elemcm 

integer, dimension (:), allocatable :: n_pp,n_nn,n_pn
double precision, dimension (:,:), allocatable :: v_pp,v_nn,v_pn
integer, dimension (:,:), allocatable :: icol_nn,irow_nn,icol_pp,irow_pp,icol_pn,irow_pn
integer, dimension (:,:), allocatable :: irowc_nn,irowc_pp,irowc_pn,irowe_nn,irowe_pp,irowe_pn

!imax=28 ! maximum number of s.p levels
!jmax=24
!call create_pair(imax,n_dim,i_pair,i_bs,j_bs)

if (iflag.eq.0) open (7,file='vlk.dat',status='old', form='unformatted')
if (iflag.eq.1) open (7,file='vlk_hfb.dat',status='old', form='unformatted')
if (iflag.eq.2) open (7,file='vlk_dd.dat',status='old', form='unformatted')
if (iflag.eq.3) open (7,file='vlk_dd_hfb.dat',status='old', form='unformatted')
if (iflag.eq.4) open (7,file='vlk_nat_orb.dat',status='old', form='unformatted')

i_poz=0
i_pozr=0
ndimv=int(4.0d0*(247110.26983*exp(0.03520*imax)))
!write(*,*)'ndimv',ndimv
!ndimv=30000000
!ndimv2=4000000
ndimv2=int(dfloat(ndimv)/2.0d0)
!ndimr=100000
ndimr=imax*imax
!ndimv=30000000
!ndimv2=4000000
!ndimv2=int(ndimv/3)
!ndimr=100000
!allocate(val_ar(ndimv),i_col(ndimv),i_row(ndimr))
if (.not.allocated(v_pp)) then
allocate(v_pp(0:jmax,ndimv2),v_nn(0:jmax,ndimv2),v_pn(0:jmax,ndimv))
allocate(n_pp(0:jmax),n_nn(0:jmax),n_pn(0:jmax))
allocate(icol_nn(0:jmax,ndimv2),icol_pp(0:jmax,ndimv2),icol_pn(0:jmax,ndimv)) 
allocate(irowc_nn(0:jmax,ndimr))
allocate(irowe_nn(0:jmax,ndimr))
allocate(irowc_pp(0:jmax,ndimr))
allocate(irowe_pp(0:jmax,ndimr))
allocate(irowc_pn(0:jmax,ndimr))
allocate(irowe_pn(0:jmax,ndimr))
endif

if (.not.allocated(irow_pp)) then 
allocate(irow_nn(0:jmax,ndimv2),irow_pp(0:jmax,ndimv2),irow_pn(0:jmax,ndimv))
endif


v_nn=0.0d0
v_pp=0.0d0
v_pn=0.0d0
n_pp=0
n_nn=0
n_pn=0
icol_nn=0
irow_nn=0
irowc_nn=0
irowe_nn=0

icol_pp=0
irow_pp=0
irowc_pp=0
irowe_pp=0

icol_pn=0
irow_pn=0
irowc_pn=0
irowe_pn=0


in_nn=0
in_pp=0
in_pn=0

do while (.not.eof(7))
if (iflag.eq.0) then 
!  read(7)it,ip,ij,i1,i2,i3,i4,v_elem,v_elemcm
  read(7)it,ij,i1,i2,i3,i4,v_elem,v_elemcm
  v_elem=v_elem-v_elemcm*hbarom/dble(iz+in)
endif
if (iflag.ne.0) read(7)it,ij,i1,i2,i3,i4,v_elem


if((i1.le.2*imax.and.i2.le.2*imax).and.(i3.le.2*imax.and.i4.le.2*imax)) then



!write(77,'(7i5,2f15.10)')it,ip,ij,i1,i2,i3,i4,v_elem,v_elemcm


 
if (it==1) then 
in_nn=in_nn+1
n_nn(ij/2)=n_nn(ij/2)+1

if (n_nn(ij/2) > ndimv2) then
write(*,*)'Increase ndimv2!'
stop
endif


 v_nn(ij/2,n_nn(ij/2))=v_elem
 i_j=imax*(i1/2-1)+i2/2
 k_l=imax*(i3/2-1)+i4/2
 icol_nn(ij/2,n_nn(ij/2))=k_l
 irow_nn(ij/2,n_nn(ij/2))=i_j
!write(77,'(7i4,f15.10)')it,ip,ij,i1,i2,i3,i4,v_elem
endif

if (it==-1) then 
in_pp=in_pp+1

 n_pp(ij/2)=n_pp(ij/2)+1

if (n_pp(ij/2) > ndimv2) then
write(*,*)'Increase ndimv2!'
stop
endif


 v_pp(ij/2,n_pp(ij/2))=v_elem
 i_j=imax*((i1+1)/2-1)+(i2+1)/2
 k_l=imax*((i3+1)/2-1)+(i4+1)/2
 icol_pp(ij/2,n_pp(ij/2))=k_l
 irow_pp(ij/2,n_pp(ij/2))=i_j
!write(88,'(7i4,f15.10)')it,ip,ij,i1,i2,i3,i4,v_elem
endif

if (it==0) then
in_pn=in_pn+1

 n_pn(ij/2)=n_pn(ij/2)+1

if (n_pn(ij/2) > ndimv) then
write(*,*)'Increase ndimv!'
stop
endif

 v_pn(ij/2,n_pn(ij/2))=v_elem
 i_j=imax*((i1+1)/2-1)+(i2)/2
 k_l=imax*((i3+1)/2-1)+(i4)/2
 icol_pn(ij/2,n_pn(ij/2))=k_l
 irow_pn(ij/2,n_pn(ij/2))=i_j
!write(99,'(7i4,f15.10)')it,ip,ij,i1,i2,i3,i4,v_elem
endif

!if (it == ipart .and. ij==jcp) then
!write(77,'(7i4,f15.10)')it,ip,ij,i1,i2,i3,i4,v_elem
!call crs_store(i_j,k_l,v_elem,i_poz,i_pozr,val_ar,i_col,i_row)
!endif


endif
enddo

if (iflag.le.4) write(*,*)'Number of V     nn,pp,np matrix elements:',in_nn,in_pp,in_pn
!if (iflag.eq.2) write(*,*)'Number of V_DD  nn,pp,np matrix elements:',in_nn,in_pp,in_pn
close(7)

call crs_store(jmax,n_nn,irow_nn,icol_nn,v_nn,irowc_nn,irowe_nn)
call crs_store(jmax,n_pp,irow_pp,icol_pp,v_pp,irowc_pp,irowe_pp)
call crs_store(jmax,n_pn,irow_pn,icol_pn,v_pn,irowc_pn,irowe_pn)

return
end subroutine read_inter

!--------------------------------

subroutine crs_store(jmax,nelem,i_row,i_col,val_ar,i_rowc,i_rowe)

implicit none

integer :: i,i_prev,i_poz,i_pozr,jcal
integer(kind=1) :: jmax
!double precision :: val 

integer, dimension (:), allocatable :: nelem
integer, dimension (:,:), allocatable :: i_col,i_row,i_rowc,i_rowe
double precision, dimension (:,:), allocatable :: val_ar


do jcal=0,jmax
   i_prev=1
   i_rowc(jcal,1)=i_prev
  do i_poz=1,nelem(jcal)

    i=i_row(jcal,i_poz)
    if (i.ne.i_prev) then 
       i_rowc(jcal,i)=i_poz    
       i_rowe(jcal,i_prev)=i_poz-1
       i_prev=i 
    endif
    if (i_poz==nelem(jcal)) i_rowe(jcal,i_row(jcal,i_poz))=i_poz
       
  enddo

enddo

return
end subroutine crs_store
!-------------------------------------------------
double precision function velem(it,ij,i1,i2,i3,i4,imax,n_nn,icol_nn,vnn,irowc_nn,irowe_nn,n_pp,icol_pp,vpp,irowc_pp,irowe_pp,n_pn,icol_pn,vpn,irowc_pn,irowe_pn)
implicit none

integer :: it,ij,i1,i2,i3,i4,imax,i_j,k_l
integer :: i_min,i_max,i

integer, dimension (:), allocatable :: n_pp,n_nn,n_pn
double precision, dimension (:,:), allocatable :: vpp,vnn,vpn
integer, dimension (:,:), allocatable :: icol_nn,icol_pp,icol_pn
integer, dimension (:,:), allocatable :: irowc_nn,irowc_pp,irowc_pn,irowe_nn,irowe_pp,irowe_pn


velem=0.0d0

if (it==1) then 
i_j=imax*(i1/2-1)+i2/2
k_l=imax*(i3/2-1)+i4/2

i_min=irowc_nn(ij,i_j)
i_max=irowe_nn(ij,i_j)

if (i_min > 0) then

do i=i_min,i_max
if (icol_nn(ij,i) == k_l) then 
 velem=vnn(ij,i)
 return
endif

enddo

endif
endif

if (it==-1) then 
i_j=imax*((i1+1)/2-1)+(i2+1)/2
k_l=imax*((i3+1)/2-1)+(i4+1)/2

i_min=irowc_pp(ij,i_j)
i_max=irowe_pp(ij,i_j)

if (i_min > 0) then 
 
do i=i_min,i_max
if (icol_pp(ij,i) == k_l) then 
 velem=vpp(ij,i)
 return
endif

enddo

endif

endif

if (it==0) then 
i_j=imax*((i1+1)/2-1)+(i2)/2
k_l=imax*((i3+1)/2-1)+(i4)/2

i_min=irowc_pn(ij,i_j)
i_max=irowe_pn(ij,i_j)

if (i_min > 0) then

do i=i_min,i_max
if (icol_pn(ij,i) == k_l) then 
 velem=vpn(ij,i)
 return
endif

enddo

endif

endif
return 
end function velem
!
double precision function vnn(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp,ipom,ifaz
integer :: ii,jj,kk,ll
double precision :: fact 

ii=i
jj=j
kk=k
ll=l
ifaz=1
fact=1.0d0
if (ii==jj) fact=dsqrt(2.d0)
if (kk==ll) fact=fact*dsqrt(2.d0)

if (i > j) then 
! ipom=ii
! ii=jj
! jj=ipom
 
 ii=j
 jj=i
 if (ihf.eq.0) ifaz=-1*(-1)**((levn(i)%j2+levn(j)%j2)/2-Jp)
 if (ihf.eq.1) ifaz=-1*(-1)**((lhfn(i)%j2+lhfn(j)%j2)/2-Jp)
 if (ihf.eq.2) ifaz=-1*(-1)**((lnon(i)%j2+lnon(j)%j2)/2-Jp)
endif

if (k > l) then 
! ipom=kk
! kk=ll
! ll=ipom
kk=l
ll=k
if (ihf.eq.0) ifaz=-1*ifaz*(-1)**((levn(k)%j2+levn(l)%j2)/2-Jp)
if (ihf.eq.1) ifaz=-1*ifaz*(-1)**((lhfn(k)%j2+lhfn(l)%j2)/2-Jp)
if (ihf.eq.2) ifaz=-1*ifaz*(-1)**((lnon(k)%j2+lnon(l)%j2)/2-Jp)

endif

if (ii > kk) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j
endif 

if (ii == kk .AND. jj > ll ) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j

endif 


vnn=dfloat(ifaz)*fact*velem(1,Jp,2*ii,2*jj,2*kk,2*ll,imax,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn,v_pn,irowc_pn,irowe_pn)
return 
end function vnn
!
double precision function vpp(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp,ipom
integer :: ii,jj,kk,ll,ifaz
double precision :: fact

ii=i
jj=j
kk=k
ll=l
ifaz=1
fact=1.0d0

if (ii==jj) fact=dsqrt(2.d0)
if (kk==ll) fact=fact*dsqrt(2.d0)

if (i > j) then 
 !ipom=ii
 !ii=jj
 !jj=ipom
 ii=j
 jj=i
 
 if (ihf.eq.0) ifaz=-1*(-1)**((levp(i)%j2+levp(j)%j2)/2-Jp)
 if (ihf.eq.1) ifaz=-1*(-1)**((lhfp(i)%j2+lhfp(j)%j2)/2-Jp)
 if (ihf.eq.2) ifaz=-1*(-1)**((lnop(i)%j2+lnop(j)%j2)/2-Jp)
endif

if (k > l) then 
! ipom=kk
! kk=ll
! ll=ipom
kk=l
ll=k
if (ihf.eq.0) ifaz=-1*ifaz*(-1)**((levp(k)%j2+levp(l)%j2)/2-Jp)
if (ihf.eq.1) ifaz=-1*ifaz*(-1)**((lhfp(k)%j2+lhfp(l)%j2)/2-Jp)
if (ihf.eq.2) ifaz=-1*ifaz*(-1)**((lnop(k)%j2+lnop(l)%j2)/2-Jp)
endif

if (ii > kk) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j
endif 

if (ii == kk .AND. jj > ll ) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom


!ii=k
!kk=i
!jj=l
!ll=j


endif 


vpp=dfloat(ifaz)*fact*velem(-1,Jp,2*ii-1,2*jj-1,2*kk-1,2*ll-1,imax,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn,v_pn,irowc_pn,irowe_pn)
return 
end function vpp
!
double precision function vpn(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp,ipom
integer :: ii,jj,kk,ll

ii=i
jj=j
kk=k
ll=l


if (ii > kk) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j
endif 

if (ii == kk .AND. jj > ll ) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom


!ii=k
!kk=i
!jj=l
!ll=j


endif 



vpn=velem(0,Jp,2*ii-1,2*jj,2*kk-1,2*ll,imax,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn,v_pn,irowc_pn,irowe_pn)
return 
end function vpn


double precision function vnn_dd(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp,ipom,ifaz
integer :: ii,jj,kk,ll
double precision :: fact 

ii=i
jj=j
kk=k
ll=l
ifaz=1
fact=1.0d0
!if (ii==jj) fact=dsqrt(2.d0)
!if (kk==ll) fact=fact*dsqrt(2.d0)

if (i > j) then 
! ipom=ii
! ii=jj
! jj=ipom
 
 ii=j
 jj=i
 if (ihf.eq.0) ifaz=-1*(-1)**((levn(i)%j2+levn(j)%j2)/2-Jp)
 if (ihf.eq.1) ifaz=-1*(-1)**((lhfn(i)%j2+lhfn(j)%j2)/2-Jp)
endif

if (k > l) then 
! ipom=kk
! kk=ll
! ll=ipom
kk=l
ll=k
if (ihf.eq.0) ifaz=-1*ifaz*(-1)**((levn(k)%j2+levn(l)%j2)/2-Jp)
if (ihf.eq.1) ifaz=-1*ifaz*(-1)**((lhfn(k)%j2+lhfn(l)%j2)/2-Jp)
endif

if (ii > kk) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j
endif 

if (ii == kk .AND. jj > ll ) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j

endif 


vnn_dd=dfloat(ifaz)*fact*velem(1,Jp,2*ii,2*jj,2*kk,2*ll,imax,n_nn_dd,icol_nn_dd,v_nn_dd,irowc_nn_dd,irowe_nn_dd,n_pp_dd,icol_pp_dd,v_pp_dd,irowc_pp_dd,irowe_pp_dd,n_pn_dd,icol_pn_dd,v_pn_dd,irowc_pn_dd,irowe_pn_dd)
return 
end function vnn_dd
!
double precision function vpp_dd(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp,ipom
integer :: ii,jj,kk,ll,ifaz
double precision :: fact

ii=i
jj=j
kk=k
ll=l
ifaz=1
fact=1.0d0

!if (ii==jj) fact=dsqrt(2.d0)
!if (kk==ll) fact=fact*dsqrt(2.d0)

if (i > j) then 
 !ipom=ii
 !ii=jj
 !jj=ipom
 ii=j
 jj=i
 
 if (ihf.eq.0) ifaz=-1*(-1)**((levp(i)%j2+levp(j)%j2)/2-Jp)
 if (ihf.eq.1) ifaz=-1*(-1)**((lhfp(i)%j2+lhfp(j)%j2)/2-Jp)
endif

if (k > l) then 
! ipom=kk
! kk=ll
! ll=ipom
kk=l
ll=k
if (ihf.eq.0) ifaz=-1*ifaz*(-1)**((levp(k)%j2+levp(l)%j2)/2-Jp)
if (ihf.eq.1) ifaz=-1*ifaz*(-1)**((lhfp(k)%j2+lhfp(l)%j2)/2-Jp)
endif

if (ii > kk) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j
endif 

if (ii == kk .AND. jj > ll ) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom


!ii=k
!kk=i
!jj=l
!ll=j


endif 


vpp_dd=dfloat(ifaz)*fact*velem(-1,Jp,2*ii-1,2*jj-1,2*kk-1,2*ll-1,imax,n_nn_dd,icol_nn_dd,v_nn_dd,irowc_nn_dd,irowe_nn_dd,n_pp_dd,icol_pp_dd,v_pp_dd,irowc_pp_dd,irowe_pp_dd,n_pn_dd,icol_pn_dd,v_pn_dd,irowc_pn_dd,irowe_pn_dd)
return 
end function vpp_dd
!
double precision function vpn_dd(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp,ipom
integer :: ii,jj,kk,ll

ii=i
jj=j
kk=k
ll=l


if (ii > kk) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom

!ii=k
!kk=i
!jj=l
!ll=j
endif 

if (ii == kk .AND. jj > ll ) then
ipom=ii
ii=kk
kk=ipom
ipom=jj
jj=ll
ll=ipom


!ii=k
!kk=i
!jj=l
!ll=j


endif 



vpn_dd=velem(0,Jp,2*ii-1,2*jj,2*kk-1,2*ll,imax,n_nn_dd,icol_nn_dd,v_nn_dd,irowc_nn_dd,irowe_nn_dd,n_pp_dd,icol_pp_dd,v_pp_dd,irowc_pp_dd,irowe_pp_dd,n_pn_dd,icol_pn_dd,v_pn_dd,irowc_pn_dd,irowe_pn_dd)
return 
end function vpn_dd


double precision function vpn_tda(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp

vpn_tda=vpn(i,j,k,l,Jp)+3.0d0*vpn_dd(i,j,k,l,Jp)

return 
end function vpn_tda

double precision function vnn_tda(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp

vnn_tda=vnn(i,j,k,l,Jp)+3.0d0*vnn_dd(i,j,k,l,Jp)

return 
end function vnn_tda


double precision function vpp_tda(i,j,k,l,Jp)

implicit none

integer :: i,j,k,l,Jp

vpp_tda=vpp(i,j,k,l,Jp)+3.0d0*vpp_dd(i,j,k,l,Jp)

return 
end function vpp_tda

subroutine precal_C1_C2(i1,cgg1_int,cgg2_int)

integer :: i1,i,j,k
real (kind=8):: clebsch
real (kind=8), dimension (:,:,:), allocatable:: cgg1_int,cgg2_int

allocate (cgg1_int(1:i1,1:i1,0:i1))
cgg1_int=0.0d0
allocate (cgg2_int(1:i1,1:i1,0:i1))
cgg2_int=0.0d0


do i=1,i1,2
 do j=1,i1,2
  do k=0,i1
    clebsch=cg_int(i,1,j,1,2*k)
    cgg1_int(i,j,k)=clebsch
    clebsch=cg_int(i,1,j,-1,2*k)
    cgg2_int(i,j,k)=clebsch
  enddo
 enddo
enddo

end subroutine precal_C1_C2



end module inter_stor
