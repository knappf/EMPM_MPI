      TYPE level_typ

        sequence

	integer                 :: npr 
        integer                 :: n
        integer                 :: l
        integer                 :: j         ! 2*j        
        double precision        :: en        ! s.p.-energy
      
      END TYPE level_typ

      TYPE amp_typ

        sequence 

        integer                 :: par        
        integer                 :: hol      
        double precision        :: am

      END TYPE amp_typ

      TYPE ro_typ

       sequence

       integer*2                 :: ib
       integer*2                 :: isi
       integer*2                 :: i1
       integer*2                 :: i2
       real (kind=4)             :: rho

      END TYPE ro_typ



     

      
