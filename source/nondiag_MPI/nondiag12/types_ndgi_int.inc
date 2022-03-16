      TYPE level_typ

        sequence

        integer                 :: n
        integer                 :: l
        integer                 :: j         ! 2*j        
        real                    :: en        ! s.p.-energy
      
      END TYPE level_typ

      TYPE amp_typ

        sequence 

        integer                 :: par        
        integer                 :: hol      
        double precision        :: am

      END TYPE amp_typ

      TYPE rho_typ

        sequence 

        integer*2                :: ilap
        integer*2                 :: j
        integer*2                 :: i1        
        integer*2                 :: i2      
!        double precision        :: ro
        real(kind=4)              :: ro 

      END TYPE rho_typ

      TYPE rho2_typ

        sequence

        integer                 :: ilap
        integer(kind=1)         :: j
        integer(kind=2)         :: i1
        integer(kind=2)         :: i2
!        double precision        :: ro
        real                    :: ro
      END TYPE rho2_typ




     

      
