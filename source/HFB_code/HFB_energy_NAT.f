      subroutine HFB_energy_NAT

       USE technical
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'


       call read_inter(4,hbarom,az,an,imax,int(jmax,1)
     &,n_nn,icol_nn,v_nn,irowc_nn,irowe_nn,n_pp
     &,icol_pp,v_pp,irowc_pp,irowe_pp,n_pn,icol_pn
     &,v_pn,irowc_pn,irowe_pn)

       E_kin=0.d0
       E_prot=0.d0
       E_neut=0.d0
       E_pn=0.d0
       E_HFB=0.d0

       do i=1,id
!        do j=1,id
         E_kin=E_kin+kin_p(i,i)*dble(lnop(i)%j2+1)*lnop(i)%vi**2.d0
     &              +kin_n(i,i)*dble(lnon(i)%j2+1)*lnon(i)%vi**2.d0
!        enddo
       enddo

       do i=1,id
!        do j=1,id
!         if(levp(i)%j2.eq.levp(j)%j2) then
         do k=1,id
!          do l=1,id
!           if(levp(k)%j2.eq.levp(l)%j2) then
            do Jp=0,jmax
            E_prot=E_prot+0.5d0*Vpp(i,k,i,k,Jp)
     &           *lnop(i)%vi**2.d0*lnop(k)%vi**2.d0*dble(2*Jp+1)
            E_neut=E_neut+0.5d0*Vnn(i,k,i,k,Jp)
     &           *lnon(i)%vi**2.d0*lnon(k)%vi**2.d0*dble(2*Jp+1)
            E_pn=E_pn+1.d0*Vpn(i,k,i,k,Jp)
     &           *lnop(i)%vi**2.d0*lnon(k)%vi**2.d0*dble(2*Jp+1)
            enddo
!           endif
!          enddo
         enddo
!         endif
!        enddo
       enddo

       E_HFB = E_kin + E_prot + E_neut + E_pn
       
!       write(*,*) 'HFB energy                 =',E_HFB,'  MeV'
!       write(*,*) 'Kinetic energy             =',E_kin,'  MeV'
!       write(*,*) 'Proton interaction energy  =',E_prot,'  MeV'
!       write(*,*) 'Neutron interaction energy =',E_neut,'  MeV'
!       write(*,*) 'P-N interaction energy     =',E_pn,'  MeV'

       open(1,file='Energy_NAT.out',status='unknown',form='formatted')
        write(1,*) E_HFB
       close(1)

       return
      end
