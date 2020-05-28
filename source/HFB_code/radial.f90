      subroutine radial_wf(irad_hf_nat)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision, allocatable, save :: argu2(:)
       double precision, allocatable, save :: argup(:)
       double precision, allocatable, save :: argun(:)
       double precision, allocatable, save :: rad_denp(:)
       double precision, allocatable, save :: rad_denn(:)
       double precision, allocatable, save :: rad_den2(:)

       double precision, allocatable, save :: OrbNOp(:,:)
       double precision, allocatable, save :: OrbNOn(:,:)

 
       character*2 :: hfnat

      if (.not.allocated(lnop)) allocate(lnop(id),lnon(id))
      
       if (irad_hf_nat.eq.1) then 
            tran_NOp=tran_p
            tran_NOn=tran_n    
            hfnat='HF'
            lnop=lhfp
            lnon=lhfn
       endif

       if (irad_hf_nat.eq.2) then
      tran_NOn=0.d0
      tran_NOp=0.d0
       do i=1,id
        do j=1,id
         do k=1,id
          tran_NOp(i,j)=tran_NOp(i,j)+tran_p(k,j)*tran_HFp(i,k)
          tran_NOn(i,j)=tran_NOn(i,j)+tran_n(k,j)*tran_HFn(i,k)
          hfnat='NO'
         enddo
        enddo
       enddo

      endif

       bos1=dsqrt(0.5d0*(mp+mn)*hbarom/(hbarc**2.d0))

       allocate(OrbNOp(0:igrid2,id),OrbNOn(0:igrid2,id))
       OrbNOp=0.d0
       OrbNOn=0.d0
        do i=0,igrid2
         radi=dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_NOp(k,j)*R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_NOn(k,j)*R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo
          OrbNOp(i,j) = val_p
          OrbNOn(i,j) = val_n
         enddo
        enddo
 
      
      open (1,file='u_'//hfnat//'_p.dat',status='unknown',form='formatted')
      do i=0,igrid2
        write(1,'(1000f10.5)')dble(i)*sizebox/dble(igrid2),((OrbNOp(i,k)*dble(i)*sizebox/dble(igrid2))**2.0d0,k=1,id)
      enddo
      close(1)

      open (1,file='R_'//hfnat//'_p.dat',status='unknown',form='formatted')
      do i=0,igrid2
        write(1,'(1000f10.5)')dble(i)*sizebox/dble(igrid2),((OrbNOp(i,k))**2.0d0,k=1,id)
      enddo
      close(1)

      open (1,file='u_'//hfnat//'_n.dat',status='unknown',form='formatted')
      do i=0,igrid2
        write(1,'(1000f10.5)')dble(i)*sizebox/dble(igrid2),((OrbNOn(i,k)*dble(i)*sizebox/dble(igrid2))**2.0d0,k=1,id)
      enddo
      close(1)

      open (1,file='R_'//hfnat//'_n.dat',status='unknown',form='formatted')
      do i=0,igrid2
        write(1,'(1000f10.5)')dble(i)*sizebox/dble(igrid2),((OrbNOn(i,k))**2.0d0,k=1,id)
      enddo
      close(1)

!      open(1,file='OrbNOp.out',status='unknown',form='formatted')
!         write(1,*) id,igrid2
!         write(1,*) (dble(i)*sizebox/dble(igrid2),i=0,igrid2)
!         do k=1,id
!           write(1,*) (OrbNOp(i,k),i=0,igrid2)
!         enddo
!      close(1)
!      open(1,file='OrbNOn.out',status='unknown',form='formatted')
!         write(1,*) id,igrid2
!         write(1,*) (dble(i)*sizebox/dble(igrid2),i=0,igrid2)
!         do k=1,id
!          write(1,*) (OrbNOn(i,k),i=0,igrid2)
!         enddo
!      close(1)

       deallocate(OrbNOp,OrbNOn)

       r4tot=0.d0
       r2tot=0.d0
       r4Z=0.d0
       r2Z=0.d0
       r4N=0.d0
       r2N=0.d0
       allocate(argu2(0:igrid2),argup(0:igrid2),argun(0:igrid2))
       argu2=0.d0
       argup=0.d0
       argun=0.d0
       allocate(rad_denp(0:igrid2),rad_denn(0:igrid2),rad_den2(0:igrid2))
!       allocate(rad_wf_p(id),rad_wf_n(id))

       rad_denp=0.d0
       rad_denn=0.d0
       rad_den2=0.d0
        do i=0,igrid2
         val1=0.d0
         val2=0.d0
         rad_wf_p=0.d0
         rad_wf_n=0.d0
         radi=dble(i)*sizebox/dble(igrid2)
         do j=1,id
          val_p=0.d0
          val_n=0.d0
          do k=1,id
           val_p=val_p+tran_NOp(k,j)*R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
           val_n=val_n+tran_NOn(k,j)*R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
          enddo



          val1=val1+lnop(j)%vi**2.d0*val_p**2.d0*dble(lnop(j)%j2+1)
          val2=val2+lnon(j)%vi**2.d0*val_n**2.d0*dble(lnon(j)%j2+1)
         enddo
         rad_denp(i)=val1/(4.d0*pi)
         rad_denn(i)=val2/(4.d0*pi)
         rad_den2(i)=(val1+val2)/(4.d0*pi)

!         write(99,'(1000f10.5)')radi,(rad_wf_p(j),j=1,id)
!         write(98,'(1000f10.5)')radi,(rad_wf_n(j),j=1,id)
        enddo

       open(47,file='rad_density_'//hfnat//'.out',status='unknown',form='formatted')
        rad_sum=0.d0
        dx=sizebox/dble(igrid2)
        do i=0,igrid2
         radi=dble(i)*sizebox/dble(igrid2)
         rad_sum=rad_sum+rad_den2(i)*dx*radi**2.d0
         write(47,*) radi,rad_den2(i)
        enddo
        write(47,*) 'Int rho = ',rad_sum*4.d0*pi
        close(47)

        open(47,file='rad_density_n_'//hfnat//'.out',status='unknown',form='formatted')
          rad_sum=0.d0
          dx=sizebox/dble(igrid2)
          do i=0,igrid2
           radi=dble(i)*sizebox/dble(igrid2)
           rad_sum=rad_sum+rad_denn(i)*dx*radi**2.d0
           write(47,*) radi,rad_denn(i)
          enddo
          write(47,*) 'Int rho N = ',rad_sum*4.d0*pi
          close(47)
          
          open(47,file='rad_density_p_'//hfnat//'.out',status='unknown',form='formatted')
          rad_sum=0.d0
          dx=sizebox/dble(igrid2)
          do i=0,igrid2
           radi=dble(i)*sizebox/dble(igrid2)
           rad_sum=rad_sum+rad_denp(i)*dx*radi**2.d0
           write(47,*) radi,rad_denp(i)
          enddo
          write(47,*) 'Int rho Z = ',rad_sum*4.d0*pi
          close(47) 


       do i=0,igrid2
        radi=dble(i)*sizebox/dble(igrid2)
        argu2(i)=radi**2.d0*rad_den2(i)
        argup(i)=radi**2.d0*rad_denp(i)
        argun(i)=radi**2.d0*rad_denn(i)
       enddo
       dx=sizebox/dble(igrid2)
       do i=0,igrid2-1
        radi=(dble(i)+0.5d0)*sizebox/dble(igrid2)
        r4tot=r4tot+dx*0.5d0*(argu2(i)+argu2(i+1))*radi**2.d0
        r2tot=r2tot+dx*0.5d0*(argu2(i)+argu2(i+1))
        r4Z=r4Z+dx*0.5d0*(argup(i)+argup(i+1))*radi**2.d0
        r2Z=r2Z+dx*0.5d0*(argup(i)+argup(i+1))
        r4N=r4N+dx*0.5d0*(argun(i)+argun(i+1))*radi**2.d0
        r2N=r2N+dx*0.5d0*(argun(i)+argun(i+1))
       enddo

       open(1,file='radial_'//hfnat//'.out',status='unknown',form='formatted')
        write(1,*) 'sqrt(<r^2>_tot)=',dsqrt(r4tot/r2tot),'fm'
        write(1,*) 'sqrt(<r^2>_p)=',dsqrt(r4Z/r2Z),'fm'
        write(1,*) 'sqrt(<r^2>_n)=',dsqrt(r4N/r2N),'fm'
       close(1)
       deallocate(argu2,argup,argun)
       deallocate(rad_denp,rad_denn,rad_den2)

       return
      end
