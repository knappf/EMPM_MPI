program convert_to_real
      implicit double precision (a-h,o-z)

      call readfn_conv
end
!******************************************************************************
     subroutine readfn_conv

      implicit double precision (a-h,o-z)

      TYPE rho_typ

        sequence

        integer*2                 :: ilap
        integer*2                 :: j
        integer*2                 :: i1
        integer*2                 :: i2
        double precision        :: ro

      END TYPE rho_typ


      double precision, dimension(:,:,:,:), allocatable :: ron
      type(rho_typ), dimension(:), allocatable :: ronn

      character*30 name_old,name_new

      name_old='V_phon_p_p.dat'
      name_new='Vr_phon_p_p.dat'

 
      ifile=3

      open(ifile,file=name_old,status='old',form='unformatted')

      open(ifile+1,file=name_new,status='unknown',form='unformatted')


      ndro=5000000
      ndgg=0

      if (.not.allocated(ronn)) allocate (ronn(ndro))


       do while (.not.eof(ifile))

       read(ifile)igg,ndgg


       read(ifile)(ronn(ii)%ilap,ii=1,ndgg)
       read(ifile)(ronn(ii)%j,ii=1,ndgg)
       read(ifile)(ronn(ii)%i1,ii=1,ndgg)
       read(ifile)(ronn(ii)%i2,ii=1,ndgg)
       read(ifile)(ronn(ii)%ro,ii=1,ndgg)

       write(ifile+1)(ronn(ii)%ilap,ii=1,ndgg)
       write(ifile+1)(ronn(ii)%j,ii=1,ndgg)
       write(ifile+1)(ronn(ii)%i1,ii=1,ndgg)
       write(ifile+1)(ronn(ii)%i2,ii=1,ndgg)
       write(ifile+1)(real(ronn(ii)%ro),ii=1,ndgg)



       enddo

       close(ifile)
       close(ifile+1)
       return

       end subroutine readfn_conv

!***********************************************************************


