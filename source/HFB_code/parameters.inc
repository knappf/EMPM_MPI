
       parameter (pi=3.141592653589793d0)
       parameter (hbarc=197.3269680000d0)        ! hbar*c [MeV*fm]
       double precision :: mp,mn
       parameter (mp=938.2723100000000d0)        ! proton mass
       parameter (mn=939.5656300000000d0)        ! neutron mass
       parameter (charge2=1.4398840000d0)        ! electric charge to the power of 2

   

       INTEGER :: AZ, AN        ! proton and neutron number
       INTEGER :: noscmax       ! the maximal oscillator shell (defines the configuration space: #shells = noscmax + 1) for HFB
       INTEGER :: nosc_TDA       ! the maximal oscillator shell (defines the configuration space: #shells = noscmax + 1) for TDA
       INTEGER :: if_QTDA       ! = 0 for particle-hole TDA, = 1 for QTDA
       INTEGER :: if_ort        ! = 0 without TDA orthogonalization, = 1 with TDA orthogonalization
       INTEGER :: igrid,igrid2
       INTEGER :: min_p,max_p,min_n,max_n
       INTEGER :: if_NAT


       DOUBLE PRECISION :: precis  ! parameter to stop the iteration step
       DOUBLE PRECISION :: hbarom  ! hbar*omega [MeV]
       DOUBLE PRECISION :: quenp,quenn,quenpn,F_0
       DOUBLE PRECISION :: VDD_S0,VDD_S1   ! the strength of the Density Dependent interaction
       DOUBLE PRECISION :: vso             ! the strength of the spin-orbit term
       DOUBLE PRECISION :: sizebox

       logical :: ifp_hfb, ifn_hfb

       common /input/ AZ, AN, noscmax, nosc_TDA, if_QTDA, if_ort, if_NAT
       common /spacetda/ min_p,max_p,min_n,max_n
       common /hb/ hbarom
       common /pre/ precis
       common /logi/ ifp_hfb, ifn_hfb
       common /qu/ quenp,quenn,quenpn,F_0
       common /DD/ VDD_S0,VDD_S1,sizebox,igrid,igrid2
       common /spinorb/ vso

!       b = sqrt(m*omega/hbar) = sqrt(m*c^2*hbar*omega)/(hbar*c) = sqrt(938.5*41*A^(-1/3) MeV^2)/(197.327 MeV*fm)
