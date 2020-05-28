module types_eqm

    TYPE level_typ

!    sequence

    integer                 :: npr
    integer                 :: n
    integer                 :: l
    integer                 :: j         ! 2*j        
    double precision        :: en        ! s.p.-energy
  
  END TYPE level_typ

  
  TYPE phon_typ

!        sequence 
    integer                 :: par        
    integer                 :: j      
    double precision        :: enf
    integer                 :: us
    integer                 :: ihom

  END TYPE phon_typ

  TYPE phonbase_typ

!        sequence

    integer                 :: ila
    integer                 :: ilap
    integer                 :: spur
  END TYPE phonbase_typ



  TYPE amp_typ

!        sequence 

    integer                 :: par        
    integer                 :: hol      
!        double precision        :: am
    real                    :: am 

  END TYPE amp_typ

  TYPE rho_typ

!        sequence 

    integer*2                 :: ilap
    integer*2                 :: j
    integer*2                 :: i1        
    integer*2                 :: i2      
!        double precision        :: ro
    real                      :: ro

  END TYPE rho_typ

  TYPE rho2_typ

!        sequence

    integer                 :: ilap
    integer (kind=1)        :: j
    integer (kind=2)        :: i1
    integer (kind=2)        :: i2
    !real(kind=4)            :: ro
  END TYPE rho2_typ

  TYPE amp2_typ

!      sequence

  integer                 :: is
  integer                 :: ig
  double precision        :: am
!      real                    :: am 

  END TYPE amp2_typ

end module types_eqm
