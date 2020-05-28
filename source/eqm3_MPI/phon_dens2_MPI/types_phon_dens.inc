        TYPE level_typ

!        sequence

        integer                 :: npr
        integer                 :: n
        integer                 :: l
        integer                 :: j         ! 2*j
        double precision        :: en        ! s.p.-energy

      END TYPE level_typ

      TYPE amp2_typ

!        sequence
      
        integer                 :: is
        integer                 :: ig     
        double precision   :: am

      END TYPE amp2_typ


      TYPE roc_typ       
!      sequence
       integer                    :: ib
       integer (kind=1)           :: is
       integer (kind=2)           :: i1
       integer (kind=2)           :: i2
       double precision           :: rho
  
      END TYPE roc_typ
      
      TYPE phon_typ
!      sequence 
        integer                 :: par        
        integer                 :: jj      
        double precision        :: enf
      END TYPE phon_typ

      TYPE rho_typ

!      sequence

      integer*2                 :: ilap
      integer*2                 :: j
      integer*2                 :: i1
      integer*2                 :: i2
!      double precision        :: ro
      real (kind=4)             :: ro                      
      
      END TYPE rho_typ

     TYPE ro_real_typ

!       sequence

       integer                   :: ib
       integer (kind=1)           :: is
       integer (kind=2)           :: i1
       integer (kind=2)           :: i2

       real(kind=4)        :: rho

      END TYPE ro_real_typ





     

      
