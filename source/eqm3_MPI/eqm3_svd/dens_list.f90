
module dens_list
    contains

    subroutine create_dens_list

        implicit none

        integer :: i,imax
        logical, dimension(:), allocatable :: dens_calc

        imax=1000000
        allocate(dens_calc(imax))
        dens_calc=.false.

     
        open(62,file='2_phon_dens_list.dat',status='unknown',form='formatted')
        do while (.not.eof(62))
            read(62,*)i
            dens_calc(i)=.true.
        enddo
        close(62)



        open(62,file='2_phon_dens_calc.dat',status='old',form='formatted')
        do while (.not.eof(62))
            read(62,*)i
            dens_calc(i)=.false.
        enddo
        close(62)
        
        call execute_command_line('rm 2_phon_dens_calc_new.dat')

        open(62,file='2_phon_dens_calc_new.dat',status='unknown',form='formatted')
        
        do i=1,imax
            if (dens_calc(i).eq..true.) write(62,*)i
        enddo

        close(62)

    end subroutine create_dens_list

end module dens_list