      subroutine HFB_energy

       USE technical
       use int_arr
       USE inter_stor

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       E_kin=0.d0
       E_prot=0.d0
       E_neut=0.d0
       E_pn=0.d0
       E_pair=0.d0
       E_HFB=0.d0

       do i=1,id
        do j=1,id
         E_kin=E_kin+kin_p(i,j)*rhop_HFB(j,i)
     &                *dsqrt(dble((levp(i)%j2+1)*(levp(j)%j2+1)))
     &              +kin_n(i,j)*rhon_HFB(j,i)
     &                *dsqrt(dble((levn(i)%j2+1)*(levn(j)%j2+1)))
        enddo
       enddo

       do i=1,id
        do j=1,id
         if(levp(i)%j2.eq.levp(j)%j2) then
         do k=1,id
          do l=1,id
           if(levp(k)%j2.eq.levp(l)%j2) then
            E_prot=E_prot+0.5d0*(Fpp(i,j,k,l)+Fpp_DD(i,j,k,l))
     &                                    *rhop_HFB(l,k)*rhop_HFB(j,i)
     &        *dsqrt(dble((levp(i)%j2+1)*(levp(k)%j2+1)))
            E_neut=E_neut+0.5d0*(Fnn(i,j,k,l)+Fnn_DD(i,j,k,l))
     &                                    *rhon_HFB(l,k)*rhon_HFB(j,i)
     &        *dsqrt(dble((levn(i)%j2+1)*(levn(k)%j2+1)))
            E_pn=E_pn+1.d0*(Fpn(i,j,k,l)+Fpn_DD(i,j,k,l))
     &        *rhop_HFB(j,i)*rhon_HFB(l,k)
     &        *dsqrt(dble((levp(i)%j2+1)*(levn(k)%j2+1)))
           endif
          enddo
         enddo
         endif
        enddo
       enddo

       do i=1,id
        do j=1,id
         do k=1,id
          do l=1,id
          
           if (if_dd.ne.0) then
           if(ifp_hfb) E_pair=E_pair+0.25d0
     &                       *(Vpp(i,j,k,l,0)+Vpp_DD(i,j,k,l,0)/3.d0)
     &                                 *kapp_HFB(i,j)*kapp_HFB(l,k)
     &        *dsqrt(dble((levp(i)%j2+1)*(levp(k)%j2+1)))
           if(ifn_hfb) E_pair=E_pair+0.25d0
     &                       *(Vnn(i,j,k,l,0)+Vnn_DD(i,j,k,l,0)/3.d0)
     &                                 *kapn_HFB(i,j)*kapn_HFB(l,k)
     &        *dsqrt(dble((levn(i)%j2+1)*(levn(k)%j2+1)))
           else 
                      if(ifp_hfb) E_pair=E_pair+0.25d0
     &                       *(Vpp(i,j,k,l,0))
     &                                 *kapp_HFB(i,j)*kapp_HFB(l,k)
     &        *dsqrt(dble((levp(i)%j2+1)*(levp(k)%j2+1)))
           if(ifn_hfb) E_pair=E_pair+0.25d0
     &                       *(Vnn(i,j,k,l,0))
     &                                 *kapn_HFB(i,j)*kapn_HFB(l,k)
     &        *dsqrt(dble((levn(i)%j2+1)*(levn(k)%j2+1)))
           
           
           endif
          enddo
         enddo
        enddo
       enddo

       E_HFB = E_kin + E_prot + E_neut + E_pn + E_pair

!       write(*,*) 'HFB energy                 =',E_HFB,'  MeV'
!       write(*,*) 'Kinetic energy             =',E_kin,'  MeV'
!       write(*,*) 'Proton interaction energy  =',E_prot,'  MeV'
!       write(*,*) 'Neutron interaction energy =',E_neut,'  MeV'
!       write(*,*) 'P-N interaction energy     =',E_pn,'  MeV'
!       write(*,*) 'Pairing energy             =',E_pair,'  MeV'

       return
      end
