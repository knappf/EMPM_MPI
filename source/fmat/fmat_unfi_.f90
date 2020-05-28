!*     Program FMAT generates F-matrix elements of two-body interaction
!*     input rom HFB code       

!*     last modification 21.8.2019
      program fmatr

      use anglib   ! angular momentum staff

      implicit double precision (a-h,o-z)

      include 'formats_fmat.inc'
      include 'types_fmat.inc'

      integer (kind=1) :: ip,ij,it
      integer (kind=2) :: ia,ib,ic,id
      real (kind=8) :: vv,fmat

      character(len=30)file_int,file_sp,file_fp,file_fn,file_fpn,file_fnp

      type(level_typ),dimension(:),allocatable :: lev,levn,levp

      integer, dimension(:),allocatable :: isp_prot_ord,isp_neut_ord
      double precision, dimension(:), allocatable :: vpl

      double precision, dimension(:,:,:,:,:), allocatable :: vp,vn,vpn

      LOGICAL :: swapped

      

      neerno=0 ! error info 

      allocate(vpl(10))
      vpl=0.d0


      open(3,status='unknown',file='Input_files',form='formatted')

      write(*,*)' Name of interaction file'
      read(3,*)
      read(3,'(A )')file_int
      write(*,*)file_int
      write(*,*)
      write(*,*)' Name of output proton interaction file'
      read(3,*)
      read(3,'(A )')file_fp
      write(*,*)file_fp
      write(*,*)
      write(*,*)' Name of output neutron interaction file'
      read(3,*)
      read(3,'(A )')file_fn
      write(*,*)file_fn
      write(*,*)
      write(*,*)' Name of output proton-neutron interaction file'
      read(3,*)
      read(3,'(A )')file_fpn
      write(*,*)file_fpn
      write(*,*)
      write(*,*)' Name of output neutron-proton interaction file'
      read(3,*)
      read(3,'(A )')file_fnp
      write(*,*)file_fnp
      write(*,*)

      close(3)

      nlevdim=10000
      allocate(lev(nlevdim),levn(nlevdim),levp(nlevdim))
      lev%n=0
      levn%n=0
      levp%n=0
      jmax=0

      iii=0
      open(1,file='NAT_p.out',status='unknown',form='formatted')
      read(1,*)
        do while (.not.eof(1))
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

        levp(ii)%l=ll
        levp(ii)%j=jj
        levp(ii)%en=ee
         if (jj.gt.jmax) jmax=jj

        enddo
       close(1)

       nlevmaxp=iii
       write(*,*)' Number of proton  levels',nlevmaxp

       iii=0
       open(1,file='NAT_n.out',status='unknown',form='formatted')
       read(1,*)
        do while (.not.eof(1))
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

        levn(ii)%l=ll
        levn(ii)%j=jj
        levn(ii)%en=ee

        if (jj.gt.jmax) jmax=jj

        end do
       close(1)

      nlevmaxn=iii

      write(*,*)' Number of neutron levels',nlevmaxn
      write(*,*)'Maximal J',jmax
      

      nlev=nlevmaxn
      if (nlevmap.gt.nlev) nlev=nlevmaxp


      allocate(isp_neut_ord(nlevmaxp))
      allocate(isp_prot_ord(nlevmaxn))

      do i = 1, nlevmaxp
        isp_prot_ord(i)=i        
      end do
      
      do i = 1, nlevmaxn
        isp_neut_ord(i)=i        
      end do
      
!   sortring of single-particle states "holes" first
      do j=nlevmaxp-1,1,-1
        swapped=.false.
        do i=1,j
         if (levp(j)%en.lt.levp(j+1)%en) then
                itemp=isp_prot_ord(i)              
                isp_prot_ord(i)=isp_prot_ord(i+1)
                isp_prot_ord(i+1)=itemp
                swapped=.true.
         endif
         IF (.NOT. swapped) EXIT        
        enddo
      enddo

      do i = 1, nlevmaxp
              
      end do

      do j=nlevmaxn-1,1,-1
        swapped=.false.
        do i=1,j
         if (levn(j)%en.lt.levn(j+1)%en) then
                itemp=isp_neut_ord(i)              
                isp_neut_ord(i)=isp_neut_ord(i+1)
                isp_neut_ord(i+1)=itemp
                swapped=.true.
         endif
         IF (.NOT. swapped) EXIT        
        enddo
      enddo


      open(1,file='singpart_coup.dat',status='unknown',form='formatted')

      write(1,*)'neutrons n l 2j    hole/part               protons  n l 2j    hole/par'
      do i=1,nlev
!       write(1,12)lev(2*i)%n,lev(2*i)%l,lev(2*i)%j,lev(2*i)%en,lev(2*i-1)%en
       write(1,'(3i5,7x,f15.5,3i5,7x,f15.5)')levn(isp_neut_ord(i))%n,levn(isp_neut_ord(i))%l,levn(isp_neut_ord(i))%j,levn(isp_neut_ord(i))%en,levp(isp_prot_ord(i))%n,levp(isp_prot_ord(i))%l,levp(isp_prot_ord(i))%j,levp(isp_prot_ord(i))%en

      enddo


!     reordering of single-particle basis?


      write(*,*)'Proton states included'
      read(*,*)nlevminpc,nlevmaxpc
      write(*,*)nlevminpc,nlevmaxpc
      write(*,*)'Neutron states included'
      read(*,*)nlevminnc,nlevmaxnc
      write(*,*)nlevminnc,nlevmaxnc

      jmax=0
      do i=1,nlevmaxpc
      if (levp(isp_prot_ord(i))%j.gt.jmax) jmax=levp(isp_prot_ord(i))%j
      enddo

      do i=1,nlevmaxnc
       if (levn(isp_neut_ord(i))%j.gt.jmax) jmax=levn(isp_neut_ord(i))%j
      enddo
   
      write(*,*)'F matrix J_max',jmax
 

      open(2,status='unknown',file=file_int,form='unformatted')
      open(3,status='unknown',file=file_fp,form='unformatted')


      write(*,*)' Calculation of proton-proton interaction '
      write(*,*)
      jmin=0

      allocate(vp(0:jmax,1:nlevmaxpc,1:nlevmaxpc,1:nlevmaxpc,1:nlevmaxpc))
      vp=0.d0
 
      itzt=-1   ! Proton-proton interaction
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

           iaa=int(ia/2)+1
           ibb=int(ib/2)+1
           icc=int(ic/2)+1
           idd=int(id/2)+1

           if (iaa.le.nlevmaxpc.and.ibb.le.nlevmaxpc.and.icc.le.nlevmaxpc.and.idd.le.nlevmaxpc) then 

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

      
      nlevmin=nlevminpc
      nlevmax=nlevmaxpc


      do isi=jmin,jmax

       do ii=nlevmin,nlevmax
         write(1234,*)isi,jmax,ii,nlevmax
          jii=lev(2*ii-1)%j
        do jj=nlevmin,nlevmax
            jjj=lev(2*jj-1)%j
         do kk=nlevmin,nlevmax
             jkk=lev(2*kk-1)%j
          do ll=nlevmin,nlevmax
               jll=lev(2*ll-1)%j

               ipari=(-1)**(lev(ii)%l+lev(kk)%l)

             fmat=0.d0  
             do j=jmin,jmax

             ifaz=(-1)**((jjj+jkk)/2-isi-j) 
            sixjc=sixj(jii,jjj,2*isi,jll,jkk,2*j) 
            fmat=fmat+dfloat((2*j+1)*ifaz)*sixjc*vp(j,ii,kk,jj,ll)

             enddo
          
!        if (dabs(fmat).gt.1.d-10) write(3)-1,ipari,isi,ii,jj,kk,ll,fmat 
         if (dabs(fmat).gt.1.d-10) write(3)int(isi,1),int(ii,2),int(jj,2),int(kk,2),int(ll,2),fmat
          enddo
         enddo
        enddo
       enddo



      enddo



!      close(2)
      close(3)
      deallocate(vp)

      
      rewind(2)
 !     open(2,status='unknown',file=file_int,form='unformatted')
!c      open(3,status='unknown',file=file_fn,form='formatted')
      open(3,status='unknown',file=file_fn,form='unformatted')
      

      write(*,*)' Calculation of neutron-neutron interaction '
      write(*,*)
      jmin=0

     
      allocate(vp(0:jmax,1:nlevmaxnc,1:nlevmaxnc,1:nlevmaxnc,1:nlevmaxnc))
      vp=0.d0
 
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

           if (iaa.le.nlevmaxnc.and.ibb.le.nlevmaxnc.and.icc.le.nlevmaxnc.and.idd.le.nlevmaxnc) then

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

      write(*,*)'J_max of loaded NN interaction m.e ',jtotmax
      write(*,*)

      nlevmin=nlevminnc
      nlevmax=nlevmaxnc

      do isi=jmin,jmax

       do ii=nlevmin,nlevmax
          jii=lev(2*ii)%j
        do jj=nlevmin,nlevmax
            jjj=lev(2*jj)%j
         do kk=nlevmin,nlevmax
             jkk=lev(2*kk)%j
          do ll=nlevmin,nlevmax
               jll=lev(2*ll)%j

               ipari=(-1)**(lev(ii)%l+lev(kk)%l)

             fmat=0.d0  
             do j=jmin,jmax

             ifaz=(-1)**((jjj+jkk)/2-isi-j) 
            sixjc=sixj(jii,jjj,2*isi,jll,jkk,2*j) 
            fmat=fmat+dfloat((2*j+1)*ifaz)*sixjc*vp(j,ii,kk,jj,ll)

             enddo
          
!c          if (dabs(fmat).gt.1.d-10) 
!c     *write(3,100)1,ipari,isi,ii,jj,kk,ll,fmat 

!          if (dabs(fmat).gt.1.d-10) write(3)1,ipari,isi,ii,jj,kk,ll,fmat 
          if (dabs(fmat).gt.1.d-10) write(3)int(isi,1),int(ii,2),int(jj,2),int(kk,2),int(ll,2),fmat

 
          enddo
         enddo
        enddo
       enddo



      enddo


      close(3)
      deallocate(vp)

!****************************     
      rewind(2)
!      open(2,status='unknown',file=file_int,form='unformatted')
!c      open(3,status='unknown',file=file_fpn,form='formatted')
!c      open(4,status='unknown',file=file_fnp,form='formatted')

      open(3,status='unknown',file=file_fpn,form='unformatted')
      open(4,status='unknown',file=file_fnp,form='unformatted')


      write(*,*)' Calculation of proton-neutron interaction '
      write(*,*)
      jmin=0


      allocate(vp(0:jmax,1:nlevmaxpc,1:nlevmaxnc,1:nlevmaxpc,1:nlevmaxnc))
      vp=0.d0
 
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

           if (iaa.le.nlevmaxpc.and.ibb.le.nlevmaxnc.and.icc.le.nlevmaxpc.and.idd.le.nlevmaxnc) then
           
           vp(ij/2,iaa,ibb,icc,idd)=vv
           vp(ij/2,icc,idd,iaa,ibb)=vv
    
           endif
          endif
         endif
      enddo

      jtotmax=jtotmax/2 

      write(*,*)' J_max of loaded interaction PN  m.e ',jtotmax
      write(*,*)

      nlevmin=1

      do isi=jmin,jmax

       do ii=nlevminpc,nlevmaxpc
          jii=lev(2*ii-1)%j
        do jj=nlevminpc,nlevmaxpc
            jjj=lev(2*jj-1)%j
         do kk=nlevminnc,nlevmaxnc
             jkk=lev(2*kk)%j
          do ll=nlevminnc,nlevmaxnc
               jll=lev(2*ll)%j

               ipari=(-1)**(lev(ii)%l+lev(kk)%l)

             fmat=0.d0  
             do j=jmin,jmax

             ifaz=(-1)**((jjj+jkk)/2-isi-j) 
            sixjc=sixj(jii,jjj,2*isi,jll,jkk,2*j) 
            fmat=fmat+dfloat((2*j+1)*ifaz)*sixjc*vp(j,ii,kk,jj,ll)

             enddo
          
          if (dabs(fmat).gt.1.d-10) then
!c        write(3,100)0,ipari,isi,ii,jj,kk,ll,fmat 
!c        write(4,100)0,ipari,isi,kk,ll,ii,jj,fmat 

!        write(3)0,ipari,isi,ii,jj,kk,ll,fmat 
!        write(4)0,ipari,isi,kk,ll,ii,jj,fmat 
        write(3)int(isi,1),int(ii,2),int(jj,2),int(kk,2),int(ll,2),fmat
        write(4)int(isi,1),int(kk,2),int(ll,2),int(ii,2),int(jj,2),fmat

          endif
 
          enddo
         enddo
        enddo
       enddo


      enddo


      close(2)
      close(3)
      close(4)
      deallocate(vp)


      end
