      TYPE level_typ

        sequence

       	integer                 :: npr
        integer                 :: n
        integer                 :: l
        integer                 :: j         ! 2*j
        integer                 :: nr
        integer                 :: bl
        integer                 :: it       ! 2*tz

        double precision        :: en        ! s.p.-energy

      END TYPE level_typ 

      TYPE ph_typ

        sequence

        integer                 :: par
        integer                 :: hol

      END TYPE ph_typ

      TYPE amp_typ

        sequence

        integer                 :: par
        integer                 :: hol
        double precision   :: am

      END TYPE amp_typ
      
