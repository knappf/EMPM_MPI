module eofmod    !<-- please save this in a separate file, compile, and link with other programs
implicit none
contains

logical function EOF( unit )
    integer, optional :: unit
    integer :: ios
    character*3 :: ftm,unf

    if ( .not. present( unit ) ) stop "no unit given"
    inquire( unit, formatted=ftm,unformatted=unf) 
    read( unit, iostat=ios)
     if (ftm.eq.'YES') read( unit, *, iostat=ios )
!     if (unf.eq.'YES') read( unit, iostat=ios )
    if ( ios == -1 ) then
        EOF = .true.
    else
        EOF = .false.  
        backspace( unit )
    endif
end function

end module

