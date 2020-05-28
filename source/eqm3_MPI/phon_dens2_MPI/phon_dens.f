!Q*     
!*     Program phon_dens computes phonon densities in proton-neutron 
!*     J-coupled formalism.

!*     last update 12.7.2010

      program phon_dens

      use rdens
      use input_sp
!c      use rdens_al

      implicit double precision (a-h,o-z)

!      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'
      include 'types_phon_dens.inc'

      type(level_typ),dimension(:), allocatable :: levn,levp

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous

      character(len=16) fname

      
*     loading of input data 

      xrotrunc=1.d-10
      
      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      close(1)

      allocate(levn(ipnmx+1000),levp(ippmx+1000))

      call inp_sp(levn,levp)

      jmaxn=0
      do i=1,ippmx
       if (levp(i)%j.gt.jmaxn) jmaxn=levp(i)%j
      enddo

      do i=1,ipnmx
       if (levn(i)%j.gt.jmaxn) jmaxn=levn(i)%j
      enddo

      write(*,*)'jmax =',jmaxn


      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
!c       jphon(i)=ijj
      enddo

      close(3)


      ifmx=i ! treba upravit aby sa to nacitalo 
      jmax=jmaxn

      jamax=2*jmax
      jbmax=jamax

      isimax=jmax

      call selphon(nf,phonbs,nphon,iphous,ns1)

      call roo(jamax,jbmax,
     *ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,ipnmx,ippmn,ippmx,isimax,levn,levp)

!      call roo(14,14,9)
!       call rou(10,10,5)
!c       call test_amp1(10,10,5)
!c       call test_amp2(10,10,5)
!c      call ro1(10,10,5)
!c      call ro2(10,10,5)
!        stop

!      fname='2phonon/2f_c.dat'
!      call cutx(fname)

!      stop


      do iph=1,-1,-2
      do ityp=-1,1,2


      if (ityp.eq.-1.and.iph.eq.1) then 
              imin=ippmn
              imax=ippmx
      endif

      if (ityp.eq.-1.and.iph.eq.-1) then 
              imin=ihpmn
              imax=ihpmx
      endif

      if (ityp.eq.1.and.iph.eq.1) then 
              imin=ipnmn
              imax=ipnmx
      endif

      if (ityp.eq.1.and.iph.eq.-1) then 
              imin=ihnmn
              imax=ihnmx
      endif


!      call ropar(imin,imax,ityp,iph,14,14,9)

      enddo
      enddo

      
      end 
!*     
!*     END of the main program 
!* 

      




      
      
