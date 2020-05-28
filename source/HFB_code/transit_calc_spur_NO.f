      subroutine transit_calc_spur_NO(iipar,Jc,phh,TD,ww,idi)

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       TYPE(twoquas_type) :: phh(idi)
       double precision :: TD(idi,idi),ww(idi)

       if(iipar.eq.-1.and.Jc.eq.1) then
        open(11,file='E1_is_NO.out',form='formatted',status='unknown')
        open(12,file='E1_iv_NO.out',form='formatted',status='unknown')
        open(13,file='Squeeze_NO.out',form='formatted',status='unknown')
        open(14,file='E1_ph_NO.out',form='formatted',status='unknown')
        write(11,*) 'isE1: E_i      B_i'
        write(12,*) 'ivE1: E_i      B_i'
        write(13,*) 'Squeeze: E_i      B_i'
        write(14,*) 'physE1: E_i      B_i'
        do i=2,idi
         E1_is=0.d0
         E1_iv=0.d0
         S1=0.d0
         E1_ph=0.d0
         do j=1,idi
          if(if_QTDA.eq.0) ip1=phh(j)%q2
          if(if_QTDA.eq.0) ih1=phh(j)%q1
          if(if_QTDA.eq.1) ip1=phh(j)%q1
          if(if_QTDA.eq.1) ih1=phh(j)%q2
          itz=phh(j)%tz
          if(itz.eq.-1) ff1=(lnop(ip1)%ui*lnop(ih1)%vi
     &                           -lnop(ip1)%vi*lnop(ih1)%ui)
          if(itz.eq.1) ff1=(lnon(ip1)%ui*lnon(ih1)%vi
     &                           -lnon(ip1)%vi*lnon(ih1)%ui)
          if(if_QTDA.eq.0.and.itz.eq.-1) ff1=ff1*
     &                  (-1)**((lnop(ip1)%j2-lnop(ih1)%j2)/2)
          if(if_QTDA.eq.0.and.itz.eq.+1) ff1=ff1*
     &                  (-1)**((lnon(ip1)%j2-lnon(ih1)%j2)/2)
          if(ip1.eq.ih1) ff1=ff1/dsqrt(2.d0)
          if(itz.eq.-1) E1_is=0.5d0*trE1_p(ih1,ip1)*ff1*TD(j,i)+E1_is
          if(itz.eq.+1) E1_is=0.5d0*trE1_n(ih1,ip1)*ff1*TD(j,i)+E1_is
          if(itz.eq.-1) E1_iv=-0.5d0*trE1_p(ih1,ip1)*ff1*TD(j,i)+E1_iv
          if(itz.eq.+1) E1_iv=0.5d0*trE1_n(ih1,ip1)*ff1*TD(j,i)+E1_iv
          if(itz.eq.-1) S1=trS1_p(ih1,ip1)*ff1*TD(j,i)+S1
          if(itz.eq.+1) S1=trS1_n(ih1,ip1)*ff1*TD(j,i)+S1
          if(itz.eq.-1) E1_ph=1.0d0*trE1_p(ih1,ip1)*ff1*TD(j,i)+E1_ph
         enddo
         write(11,*) ww(i),E1_is**2.d0
         write(12,*) ww(i),E1_iv**2.d0
         write(13,*) ww(i),S1**2.d0
         write(14,*) ww(i),E1_ph**2.d0
        enddo
        close(11)
        close(12)
        close(13)
        close(14)
       endif

       return
      end
