!     last update 14.6.2013

      program nondiag
      
      use phoninteracn
      use input_sp


      implicit double precision (a-h,o-z)
      
      include 'types_ndgi_int.inc'
      include 'input_ndgi_int.inc'
      include 'formats_ndgi_int.inc'
      
      type(level_typ),dimension(:), allocatable :: levn,levp 
      
     
      xrotrunc=1.d-10
      
      write(*,*)'Loading of input '

      open(1,file='input_eqm_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      read(1,26)alfa,beta
      read(1,*)
      read(1,15)iparmn,iparmx
      read(1,15)jmin,jmax
      close(1)


      allocate(levn(ipnmx+1000),levp(ippmx+1000))
      write(*,*)'Parity J ?'
      read(*,*)ipar,jcal

      call inp_sp(levn,levp)


       call vintn1(2)!(2,ipar,jcal,levn,levp)
       call vintn2
!c       call vintn1(3,levn,levp)


      end
