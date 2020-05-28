      subroutine checkf

      implicit double precision (a-h,o-z)

      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'
      include 'types_phon_int.inc'
      

      double precision, dimension(:,:,:,:,:),allocatable ::
     *redefi

      integer, dimension(:), allocatable :: jphon
      type(level_typ),dimension(:), allocatable :: levn,levp
 
      allocate (jphon(100))
      jphon=0
      

      allocate(redefi(0:15,0:50,0:50,1:10,1:10))

      redefi=0.d0


      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      read(1,26)alfa,beta
      read(1,*)
      read(1,15)iparmn,iparmx
      read(1,15)jminn,jmaxn
      close(1)

      allocate(levn(ipnmx),levp(ippmx))


      open(1,file='singpart_coup.dat',status='old',form='formatted')

      read(1,*)
      read(1,*)

      do i=1,ipnmx
       read(1,16)nt,lt,jt,erspnt,ersppt
       levn(i)%n=nt
       levn(i)%l=lt
       levn(i)%j=jt
      enddo
      
      rewind(1)
      
      read(1,*)
      read(1,*)

      do i=1,ippmx
       read(1,16)nt,lt,jt,erspnt,ersppt
       levp(i)%n=nt
       levp(i)%l=lt
       levp(i)%j=jt
      enddo

      close(1)
      

      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
      enddo

      close(3)

      ibt=-10
      

      open(333,file='Vint_phon_p_n.dat',status='old',form='unformatted')

       do while (.not.eof(333))

111    continue 
       read(333)igg
       write(*,*)igg
       if (igg.eq.10000000) goto 112

    
       do while (ibt.ne.0)
       read(333)ibt,isit,i1t,i2t,rot

c       if (igg.eq.ig) then

       if (ibt.ne.0) then         
c       if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx
c     *.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
c         write(*,*)' Small dimensions of ron in read1',ibt,isit,i1t,i2t
c         stop
c       endif

        redefi(isit,igg,ibt,i1t,i2t)=rot       
       endif

       if (ibt.eq.0) then 
        ibt=-10
        goto 111
       endif
 

c       endif

       enddo

       enddo

112    close(333)
    
      ilamax=44
      i1mn=1
      i2mn=1
      i1mx=6
      i2mx=6
      do isi=0,15
      do ila=1,ilamax
       do ilap=1,ilamax
         do i1=i1mn,i1mx
          do i2=i2mn,i2mx
          j1=levn(i1)%j
          j2=levn(i2)%j
c          ifaz=(-1)**(jphon(ila)+jphon(ilap)+(j1+j2)/2)
          ifaz=(-1)**((j1+j2)/2)
          
      if ((redefi(isi,ila,ilap,i1,i2)
     *+dfloat(ifaz)*redefi(isi,ilap,ila,i2,i1)).gt.1.d-10) then
      write(9999,*)isi,ila,ilap,i1,i2
      write(9999,*)jphon(ila),jphon(ilap),j1,j2
      write(9999,*)redefi(isi,ila,ilap,i1,i2),
     *(dfloat((-1)*ifaz))*redefi(isi,ilap,ila,i2,i1)
      endif
          enddo
         enddo
        enddo
       enddo
       enddo


      end subroutine checkf

*********************************************************************
      subroutine checkfmat

      implicit double precision (a-h,o-z)

      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'
      include 'types_phon_int.inc'
      

      double precision, dimension(:,:,:,:,:),allocatable ::
     *fmat

      integer, dimension(:), allocatable :: jphon
      type(level_typ),dimension(:), allocatable :: levn,levp
 
      allocate (jphon(100))
      jphon=0
      

      allocate(fmat(0:15,0:50,0:50,0:50,0:50))

      fmat=0.d0


      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      read(1,26)alfa,beta
      read(1,*)
      read(1,15)iparmn,iparmx
      read(1,15)jminn,jmaxn
      close(1)

      allocate(levn(ipnmx),levp(ippmx))


      open(1,file='singpart_coup.dat',status='old',form='formatted')

      read(1,*)
      read(1,*)

      do i=1,ipnmx
       read(1,16)nt,lt,jt,erspnt,ersppt
       levn(i)%n=nt
       levn(i)%l=lt
       levn(i)%j=jt
      enddo
      
      rewind(1)
      
      read(1,*)
      read(1,*)

      do i=1,ippmx
       read(1,16)nt,lt,jt,erspnt,ersppt
       levp(i)%n=nt
       levp(i)%l=lt
       levp(i)%j=jt
      enddo

      close(1)
      

      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
      enddo

      close(3)
      

      open(333,file='fmat_p.dat',status='old',form='formatted')

      do while (.not.eof(333))

      read(333,10)itt,ipt,ijt,i,j,k,l,vint
       fmat(ijt,i,j,k,l)=vint

      enddo



112    close(333)
    
      ilamax=44
      i1mn=1
      i2mn=1
      i1mx=ipnmx
      i2mx=15
      do isi=0,15
      do i=i1mn,i1mx
       do j=i1mn,i1mx
         do k=i1mn,i1mx
          do l=i1mn,i1mx
          j1=levn(i)%j
          j2=levn(j)%j
          j3=levn(k)%j
          j4=levn(l)%j

          ifaz=(-1)**((j1+j2+j3+j4)/2)
          
      if ((fmat(isi,i,j,k,l)
     *-dfloat(ifaz)*fmat(isi,j,i,l,k)).gt.1.d-10) then
      write(8888,*)isi,i,j,k,l
      write(8888,*)j1,j2,j3,j4
      write(8888,*)fmat(isi,i,j,k,l),
     *(dfloat(ifaz))*fmat(isi,j,i,l,k)
      endif
          enddo
         enddo
        enddo
       enddo
       enddo


      end subroutine checkfmat

****************************************************************************
      subroutine checkro

      implicit double precision (a-h,o-z)

      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'
      include 'types_phon_int.inc'
      

      double precision, dimension(:,:,:,:,:),allocatable ::
     *redefi

      integer, dimension(:), allocatable :: jphon
      type(level_typ),dimension(:), allocatable :: levn,levp
 
      allocate (jphon(100))
      jphon=0
      

      allocate(redefi(0:15,0:50,0:50,1:20,1:20))

      redefi=0.d0


      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
            
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      read(1,26)alfa,beta
      read(1,*)
      read(1,15)iparmn,iparmx
      read(1,15)jminn,jmaxn
      close(1)

      allocate(levn(ipnmx),levp(ippmx))


      open(1,file='singpart_coup.dat',status='old',form='formatted')

      read(1,*)
      read(1,*)

      do i=1,ipnmx
       read(1,16)nt,lt,jt,erspnt,ersppt
       levn(i)%n=nt
       levn(i)%l=lt
       levn(i)%j=jt
      enddo
      
      rewind(1)
      
      read(1,*)
      read(1,*)

      do i=1,ippmx
       read(1,16)nt,lt,jt,erspnt,ersppt
       levp(i)%n=nt
       levp(i)%l=lt
       levp(i)%j=jt
      enddo

      close(1)
      

      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
      enddo

      close(3)

      ibt=-10
      

      open(333,file='1phonon/1f_rpp.dat'
     *,status='old',form='unformatted')

       do while (.not.eof(333))

111    continue 
       read(333)igg
       write(*,*)igg
       if (igg.eq.10000000) goto 112

    
       do while (ibt.ne.0)
       read(333)ibt,isit,i1t,i2t,rot

c       if (igg.eq.ig) then

       if (ibt.ne.0) then         
c       if (ibt.gt.ndla.or.isit.gt.nsi.or.i1t.lt.n1mn.or.i1t.gt.n1mx
c     *.or.i2t.lt.n2mn.or.i2t.gt.n2mx) then
c         write(*,*)' Small dimensions of ron in read1',ibt,isit,i1t,i2t
c         stop
c       endif

        redefi(isit,igg,ibt,i1t,i2t)=rot       
       endif

       if (ibt.eq.0) then 
        ibt=-10
        goto 111
       endif
 

c       endif

       enddo

       enddo

112    close(333)
    
      ilamax=44
      i1mn=1
      i2mn=1
      i1mx=3
      i2mx=3
      do isi=0,15
      do ila=1,ilamax
       do ilap=1,ilamax
         do i1=i1mn,i1mx
          do i2=i2mn,i2mx
          j1=levn(i1)%j
          j2=levn(i2)%j
          ifaz=(-1)**(jphon(ila)+jphon(ilap)+(j1+j2)/2)
          
      if ((redefi(isi,ila,ilap,i1,i2)
     *+dfloat(ifaz)*redefi(isi,ilap,ila,i2,i1)).gt.1.d-10) then
      write(7777,*)isi,ila,ilap,i1,i2
      write(7777,*)jphon(ila),jphon(ilap),j1,j2
      write(7777,*)redefi(isi,ila,ilap,i1,i2),
     *(dfloat((-1)*ifaz))*redefi(isi,ilap,ila,i2,i1)
      endif
          enddo
         enddo
        enddo
       enddo
       enddo


      end subroutine checkro
      
