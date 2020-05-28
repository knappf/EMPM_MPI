!*     Program tda_pn_coup.f computes TDA phonons in  J-coupled proton-neutron 
!*     formalism.
!*
!*     last update 26.2.2020 

      program tda_pn_coup
      
     

 
      use cmconst
      use input_sp
      use read_inter

      implicit double precision (a-h,o-z)

 
      include 'formats_tda_cp.inc' !Definiton of used formats 
      include 'types_tda_cp.inc'      
      include 'input_tda_cp.inc'


      integer, dimension(:), allocatable :: ipozi

      double precision, dimension (:), allocatable :: work,e 
      double precision, dimension(:,:), allocatable :: amtr,xsumph,tbase,hami,amtrold,kin_p,kin_n


      double precision, dimension(:,:,:,:,:), allocatable :: fp,fn,fpn
      double precision, dimension(:,:,:,:,:),allocatable :: vp,vn,vpn

      type(level_typ),dimension(:), allocatable :: levn,levp
      type(ph_typ),dimension(:), allocatable :: iphn,iphp

      character*1, dimension(0:40) :: orbit 
      character(len=30)file_fp,file_fn,file_fpn,file_v
      character*5 tnf
      character*1 tcm


      integer(kind=1) :: j_f
      integer(kind=2) :: i_i,i_j,i_k,i_l

        
      xtrunc=0.000000000001d0
      xrotrunc=0.000000000001d0

      allocate(xsumph(10000,-1000:1000))

      xsumph=0.d0
     
      orbit(0)='s'
      orbit(1)='p'
      orbit(2)='d'
      orbit(3)='f'
      orbit(4)='g'
      orbit(5)='h'
      orbit(6)='i'
      orbit(7)='j'
      orbit(8)='k'
      orbit(9)='l'
      orbit(10)='m'
      orbit(11)='n'
      orbit(13)='o'

      open(881,file='tda_r_overl.dat',status='unknown',form='formatted')
      open(77,file='tda_pn_cp.log',status='unknown',form='formatted') 
    
      open(94,file='phon_struct.dat',status='unknown',form='formatted')
      open(96,file='sp_levord.dat',status='unknown',form='formatted')
    

      open(743,file='1phonon/1ph_cont_cp.dat',status='unknown',form='formatted')
                
!*     loading of input data 
      
      write(*,*)'Loading of input '

      open(1,file='input_tda_coup.dat',status='old',form='formatted')
      
      
      read(1,15)ia,iz
      read(1,15)ihnmn,ihnmx
      read(1,15)ihpmn,ihpmx
      read(1,15)ipnmn,ipnmx
      read(1,15)ippmn,ippmx
      read(1,222)tcm

      alfa=1.0d0   !  reduction constant  
      iparmn=-1    !  parity range
      iparmx=1    
      jminn=0      ! minimal J
      beta=0.0d0
!      read(1,26)alfa,beta
!      read(1,*)
!      read(1,15)iparmn,iparmx
!      read(1,15)jminn,jmaxn
!      read(1,*)
!      read(1,15)ihnmnc,ihnmxc
!      read(1,15)ihpmnc,ihpmxc
!      read(1,15)ipnmnc,ipnmxc
!      read(1,15)ippmnc,ippmxc


      close(1)

      hom=1.d0
      homg=hom*beta
      homcm=beta*hom/dfloat(ia)
      ezerocm=3.d0*homg/2.d0

      allocate(levn(ipnmx+1000),levp(ippmx+1000))
      allocate(iphp(ippmx*ippmx),iphn(ipnmx*ipnmx))

      
      call inp_sp(levn,levp,jmaxn)
!      call input_kin(kin_p,kin_n)   ! natural orbitals only

      jmaxn=0
      do i=1,ippmx
       if (levp(i)%j.gt.jmaxn) jmaxn=levp(i)%j
      enddo

      do i=1,ipnmx
       if (levn(i)%j.gt.jmaxn) jmaxn=levn(i)%j
      enddo

      
      write(*,*)'Range of parity in calculation:'
      write(*,*)iparmn,iparmx
      write(*,*)'Range of J in calculation:'
      write(*,*)jminn,jmaxn

      nhomx=100
!      write(*,*)' maximum nhomx '
!      write(*,*)nhomx

      write(*,*)'Apply CM orthogonalization? y/n'
      write(*,*)tcm
 222  format(a1)
!      read(*,222)tcm

      idimn=0
      idimp=0

!* loading of proton and neutron matrix elements      
      open(34,status='unknown',file='Input_files_tda',form='formatted')

      read(34,*)
      read(34,*)
      read(34,*)
      read(34,*)
      read(34,*)
      read(34,'(A )')file_fp
      read(34,*)
      read(34,'(A )')file_fn
      read(34,*)
      read(34,'(A )')file_fpn
      read(34,*)
      read(34,'(A )')file_v
     
      close(34)
            
      

!      call read_v(file_v,ippmx,ipnmx,jmax,vp,vn,vpn)    ! natural orbitals only
      call read_f(file_fp,file_fn,file_fpn,ippmx,ipnmx,jmaxn,fp,fn,fpn) 

      open(33,file='1phonon/1f_states.dat',status='unknown',form='unformatted')


      open(2,file='1phonon/1f_cp.dat',status='unknown',form='unformatted')

      open(4,file='1phonon/1f_cn.dat',status='unknown',form='unformatted')


      nlam=0



      do ipar=iparmn,iparmx,2   ! cycle over parity
      do ijj=jminn,jmaxn        ! cycle over J
     

      hom=1.d0
!      homg=hom*beta
!      homcm=beta*hom/dfloat(ia)
!      ezerocm=3.d0*homg/2.d0      
      

      call phbase(ipar,ijj,nhomx,levn,levp,iphn,iphp,idphp,idphn)



          write(*,*)'Parity =   ',ipar
          write(*,*)'     J =   ',ijj
          write(*,*)'dimension of proton subspace = ',idphp
          write(*,*)'dimension of neutron subspace= ',idphn      
                
          idimn=idimn+idphp
          idimp=idimp+idphn


          if ((idphp+idphn).ne.0) then  


          allocate(amtr(idphp+idphn,idphp+idphn))
          amtr=0.d0

          do i=1,idphp
            ip=iphp(i)%par
            ih=iphp(i)%hol
            jp=levp(ip)%j
            jh=levp(ih)%j

            do j=1,idphp
              ipp=iphp(j)%par
              ihp=iphp(j)%hol
              jpp=levp(ipp)%j
              jhp=levp(ihp)%j
       
      amat=0.d0

      
      if (ih.eq.ihp.and.ip.eq.ipp) amat=amat+levp(ip)%en-levp(ih)%en
      
      ifz=(-1)**((levp(ipp)%j+levp(ihp)%j)/2+ijj)
      amat=amat-fp(ijj,ip,ih,ihp,ipp)*dfloat(ifz)
              amtr(i,j)=amat


            enddo
          enddo

          do i=1,idphn
            ip=iphn(i)%par
            ih=iphn(i)%hol
            jp=levn(ip)%j
            jh=levn(ih)%j

            do j=1,idphn
              ipp=iphn(j)%par
              ihp=iphn(j)%hol
              jpp=levn(ipp)%j
              jhp=levn(ihp)%j


      
      amat=0.d0

      
      if (ih.eq.ihp.and.ip.eq.ipp) amat=amat+levn(ip)%en-levn(ih)%en
          

      ifz=(-1)**((levn(ipp)%j+levn(ihp)%j)/2+ijj)
      amat=amat-fn(ijj,ip,ih,ihp,ipp)*dfloat(ifz)

      
              amtr(i+idphp,j+idphp)=amat
            enddo
          enddo

          do i=1,idphp
            ip=iphp(i)%par
            ih=iphp(i)%hol
            do j=1,idphn
              ipp=iphn(j)%par
              ihp=iphn(j)%hol

      
      amat=0.d0

      
      ifz=(-1)**((levp(ip)%j+levp(ih)%j)/2+ijj)
      amat=amat-fpn(ijj,ih,ip,ipp,ihp)*dfloat(ifz)
 
              amtr(i,j+idphp)=amat
            enddo
          enddo          
          
          do i=1,idphn
            ip=iphn(i)%par
            ih=iphn(i)%hol
            do j=1,idphp
              ipp=iphp(j)%par
              ihp=iphp(j)%hol

      
      amat=0.d0


      ifz=(-1)**((levn(ip)%j+levn(ih)%j)/2+ijj)
      amat=amat-fpn(ijj,ipp,ihp,ih,ip)*dfloat(ifz)


              amtr(i+idphp,j)=amat
            enddo
          enddo 



       
          do i=1,idphp+idphn
            do j=1,idphp+idphn
            if (dabs(amtr(i,j)-amtr(j,i)).gt.1.d-10) write(888,*)'nsym   ',i,j,amtr(i,j),amtr(j,i)
            enddo
          enddo

!          do i=1,idphp
!            do j=1,idphp
!            if (dabs(amtr(i,j)-amtr(i+idphp,j+idphp)).gt.1.d-10) write(888,*)'*********',i,j,amtr(i,j),amtr(i+idphp,j+idphp)
!            enddo
!          enddo

          ndmx=idphp+idphn

!c   control output of A matrix
 33       format(100f15.5)
!          write(999,33)          
!          write(999,33)
!          do i=1,ndmx
!          write(999,33)(amtr(i,j),j=1,ndmx)
!          enddo

!          do j=1,ndmx
!            do i=1,j-1
!              amtr(i,j)=0.d0
!            enddo
!          enddo

!c          do i=1,ndmx
!c            amtr(i,i)=amtr(i,i)!+ezerocm
!c          enddo

          lwork=26*ndmx 
          
          allocate(work(26*ndmx))
          allocate(e(ndmx))

           
          write(*,*)' ndmx = ',ndmx 
         
 
          if (ipar.eq.-1.and.ijj.eq.1.and.tcm.eq.'y') then          
          
          allocate(hami(ndmx,ndmx))
          hami=0.0d0
          
          allocate(amtrold(ndmx,ndmx))
          amtrold=0.d0
          amtrold=amtr

          call cmcreate(idphp,idphn,iphp,iphn,tbase)

          call ortog(ndmx,tbase)
          
!  control output od tbase          
!          write(777,33)          
!          write(777,33)
!          do i=1,ndmx
!          write(777,33)(tbase(i,j),j=1,ndmx)
!          enddo


          call dgemm('N','N',ndmx,ndmx,ndmx,1.0d0,amtr,ndmx,tbase,ndmx,0.0d0,hami,ndmx)
!          call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
!           call gemm(amtr,tbase,hami,'N','N',1.0d0,0.0d0)

          amtr=0.0d0

          call dgemm('T','N',ndmx,ndmx,ndmx,1.0d0,tbase,ndmx,hami,ndmx,0.0d0,amtr,ndmx)


 321      format(300f15.10)
!          write(777,*)' A matica'
!          do i=1,ndmx
!            write(777,321)(amtr(i,j),j=1,ndmx)
!          enddo

!   R decouple
          do i=2,ndmx
           amtr(i,1)=0.0d0
           amtr(1,i)=0.0d0
          enddo


          do j=1,ndmx
            do i=1,j-1
              amtr(i,j)=0.d0
            enddo
          enddo

          amtr(1,1)=10000000.0d0


          endif

!         
!          write(777,33)          
!          write(777,*)'A matr '
!          do i=1,ndmx
!          write(777,33)(amtr(i,j),j=1,ndmx)
!          enddo

 
 
          call DSYEV('V','l',ndmx,amtr,ndmx,e,WORK,LWORK,INFO )

          write(777,33)          
          write(777,*)'vl. vektory A'
          do i=1,ndmx
          write(777,33)(amtr(i,j),j=1,ndmx)
          enddo
          
 
        if (ipar.eq.-1.and.ijj.eq.1.and.tcm.eq.'y') then
          hami=0.0d0
          call dgemm('N','N',ndmx,ndmx,ndmx,1.d0,tbase,ndmx,amtr,ndmx,0.d0,hami,ndmx)

          do i=1,ndmx
           do j=1,ndmx
           amtr(i,j)=hami(i,j)
           enddo
          enddo
          
          write(777,33)          
          write(777,*)'vl. vektory v povodnej baze'
          do i=1,ndmx
          write(777,33)(amtr(i,j),j=1,ndmx)
          enddo
          

          hami=0.0d0
!         prepocet matice A do novej bazy

          call dgemm('N','N',ndmx,ndmx,ndmx,1.0d0,amtrold,ndmx,amtr,ndmx,0.0d0,hami,ndmx)

          amtrold=0.0d0

          call dgemm('T','N',ndmx,ndmx,ndmx,1.0d0,amtr,ndmx,hami,ndmx,0.0d0,amtrold,ndmx)
          

          write(777,33)          
          write(777,*)'A matr po diagonalizacii '
          do i=1,ndmx
          write(777,33)(amtrold(i,j),j=1,ndmx)
          enddo


!    test diagonalizacie s CM
!          do j=1,ndmx
!            do i=1,j-1
!              amtrold(i,j)=0.d0
!            enddo
!          enddo

!          e=0.0d0


!          call DSYEV('V','l',ndmx,amtrold,ndmx,e,WORK,LWORK,INFO )
!          write(777,*)' test diagonalizacia s CM'
!          write(777,33)(e(i),i=1,ndmx)


          
          deallocate(hami)
!          deallocate(amtrold)
          

          endif

   


!c     test ortogonality
!c        do it=1,ndmx
!c          do jt=1,ndmx
!c           xst=0.d0
!c           do ik=1,ndmx
!c             xst=xst+amtr(ik,it)*amtr(ik,jt)
!c           enddo
!c           if (dabs(xst).gt.1.d-8) write(971,*)it,jt,xst
!c          enddo
!c        enddo


!c          if (alfa.eq.0.d0) then

!c                  amtre=0.d0
!c                  do i=1,ndmx
!c                   do j=1,ndmx

!c                    amtre(i,j)=amtr(i,j)

!c                    enddo
!c                  enddo
!c                  ndimeor=ndmx

!c                  call ortog

          
!c                  do i=1,ndmx
!c                   do j=1,ndmx

!c                    amtr(i,j)=amtre(i,j)

!c                    enddo
!c                  enddo



!c          endif 


          
          
          write(77,*)'  Parity = ',ipar,'    J = ',ijj
          write(77,792)(e(ii),ii=1,ndmx)  !-135.89613000d0
          write(77,*)


      xsumph=0.0d0
      allocate(ipozi(idphp))
      ipozi=0

      do i=1,ndmx
        write(tnf,101)nlam+i
        write(33)nlam+i,ipar,ijj,e(i)
        if (ipar.eq.-1.and.ijj.eq.1.and.tcm.eq.'y') write(881,*)nlam+i,e(i),amtrold(i,ndmx)
        write(94,*)
        write(94,*)'________________________________________'
        write(94,*)
        write(94,*)'i     Parity    J        energy '
        write(94,'(3i5,f10.5)')nlam+i,ipar,ijj,e(i)
        write(94,*)
!c        open(2,file='1phonon/1f_cp.'//tnf//''
!c     *,status='unknown',form='unformatted')

        iii=0

        do ii=1,idphp
         if (dabs(amtr(ii,i)).gt.xtrunc) then 
                iii=iii+1
                ipozi(iii)=ii
         endif
        enddo

        write(2)ipar,ijj,iii


        write(2)(iphp(ipozi(j))%par,iphp(ipozi(j))%hol,amtr(ipozi(j),i),j=1,iii)
!c        close(2)

        write(94,*)'- proton p(h)^-1 -'

        do ill=1,iii
         if (dabs(amtr(ipozi(ill),i)).gt.5.d-2) then 
!         write(94,'(2i5,f10.5)')iphp(ipozi(ill))%par,
!     *iphp(ipozi(ill))%hol,
!     *amtr(ipozi(ill),i)


       write(94,'(i2,a2,i2,5x,i2,a2,i2,5x,f10.5)')levp(iphp(ipozi(ill))%par)%nr,orbit(levp(iphp(ipozi(ill))%par)%l),levp(iphp(ipozi(ill))%par)%j,levp(iphp(ipozi(ill))%hol)%nr,orbit(levp(iphp(ipozi(ill))%hol)%l),levp(iphp(ipozi(ill))%hol)%j,amtr(ipozi(ill),i)


         endif
        enddo

!        write(94,*)'- neutron p(h)^-1 -'

!        do ill=1,iii
!         if (dabs(amtr(ipozi(ill)+idphp,i)).gt.5.d-2) then
!         write(94,'(2i5,f10.5)')iphn(ipozi(ill))%par,
!     *iphn(ipozi(ill))%hol,
!     *amtr(ipozi(ill)+idphp,i)


!       write(94,'(i2,a2,i2,5x,i2,a2,i2,5x,f10.5)')levn(iphn(ipozi(ill))%par)%nr,orbit(levn(iphn(ipozi(ill))%par)%l),levn(iphn(ipozi(ill))%par)%j,levn(iphn(ipozi(ill))%hol)%nr,orbit(levn(iphn(ipozi(ill))%hol)%l),levn(iphn(ipozi(ill))%hol)%j,amtr(ipozi(ill)+idphp,i)


!         endif
!        enddo


!c       ph content

        do iph=1,iii

          ippp=iphp(ipozi(iph))%par
          ihhh=iphp(ipozi(iph))%hol
          nroz=levp(ippp)%n-levp(ihhh)%n
          xsumph(nlam+i,nroz)=xsumph(nlam+i,nroz)+amtr(ipozi(iph),i)**2.d0
        enddo

!c        write(743,*)xsum

      enddo
       
      deallocate(ipozi)


      allocate(ipozi(idphn))
      ipozi=0 

      do i=1,ndmx
        write(tnf,101)nlam+i

!c        open(2,file='1phonon/1f_cn.'//tnf//''
!c     *,status='unknown',form='unformatted')

        iii=0

        do ii=1,idphn
         if (dabs(amtr(ii+idphp,i)).gt.xtrunc) then 
                iii=iii+1
                ipozi(iii)=ii
         endif
        enddo
        
        write(94,*)'- neutron p(h)^-1 -'

        do ill=1,iii
         if (dabs(amtr(ipozi(ill)+idphp,i)).gt.5.d-2) then
!         write(94,'(2i5,f10.5)')iphn(ipozi(ill))%par,
!     *iphn(ipozi(ill))%hol,
!     *amtr(ipozi(ill)+idphp,i)


       write(94,'(i2,a2,i2,5x,i2,a2,i2,5x,f10.5)')levn(iphn(ipozi(ill))%par)%nr,orbit(levn(iphn(ipozi(ill))%par)%l),levn(iphn(ipozi(ill))%par)%j,levn(iphn(ipozi(ill))%hol)%nr,orbit(levn(iphn(ipozi(ill))%hol)%l),levn(iphn(ipozi(ill))%hol)%j,amtr(ipozi(ill)+idphp,i)


         endif
        enddo

                

        write(4)ipar,ijj,iii
        write(4)(iphn(ipozi(j))%par,iphn(ipozi(j))%hol,amtr(ipozi(j)+idphp,i),j=1,iii)
!c        close(2)


        do iph=1,iii

          ippp=iphn(ipozi(iph))%par
          ihhh=iphn(ipozi(iph))%hol
          nroz=levn(ippp)%n-levn(ihhh)%n

          xsumph(nlam+i,nroz)=xsumph(nlam+i,nroz)+amtr(ipozi(iph)+idphp,i)**2.d0
        enddo

!        write(743,742)nlam+i,ipar,ijj,e(i),
!     *xsumph(nlam+i,1),xsumph(nlam+i,2),xsumph(nlam+i,3)
!     *,xsumph(nlam+i,4),xsumph(nlam+i,5)
 
!        if (xsumph(nlam+i,1).gt.0.3d0) 
!     *write(99,*)nlam+i

!c,xsumph(nlam+i,6)
!c       *,xsumph(nlam+i,7),xsumph(nlam+i,8),xsumph(nlam+i,9)


      enddo

      deallocate(ipozi) 
     
                      
!          deallocate(amtr,work,e)

          do i=1,ndmx           
          write(743,742)nlam+i,ipar,ijj,e(i),xsumph(nlam+i,0),xsumph(nlam+i,1),xsumph(nlam+i,2),xsumph(nlam+i,3),xsumph(nlam+i,4),xsumph(nlam+i,5),xsumph(nlam+i,6),xsumph(nlam+i,1)+xsumph(nlam+i,2)+xsumph(nlam+i,3)+xsumph(nlam+i,4)+xsumph(nlam+i,0)+xsumph(nlam+i,5)+xsumph(nlam+i,6)
      
         if (xsumph(nlam+i,1).gt.0.1d0) write(91,*)nlam+i,e(i)
         if (xsumph(nlam+i,2).gt.0.1d0) write(92,*)nlam+i,e(i)
         if (xsumph(nlam+i,3).gt.0.1d0) write(93,*)nlam+i,e(i)
         if (xsumph(nlam+i,4).gt.0.1d0) write(944,*)nlam+i,e(i)
         if (xsumph(nlam+i,0).gt.0.1d0) write(90,*)nlam+i,e(i)



          enddo


          deallocate(amtr,work,e)

         
          nlam=nlam+ndmx
          
          endif
          
          
          
        enddo
      enddo

      close(33)
      close(2)
      close(4)

!      close(11)
!      close(12)
!      close(13)
!      close(14)
      close(99)

      write(*,*)'number of proton and neutron conf. ', idimn,idimp
      write(*,*)'dimension of 1 phonon space  ',  idimn+idimp

      deallocate(xsumph)



      
      end 
!*     
!*     END of the main program 
!* 


      
!*     subrouutines and functions       


      
!**************************************************************************
!*     phbase generates particle-hole combinations with parity 
!*     ipar and J and stores indices ip and ih to arrays
!*     iphp, iphn for protons, neutrons respectively. 
!************************************************************************** 
      subroutine phbase(ipar,ijj,nhomx,levne,levpe,iphne,iphpe,idphp,idphn)

      implicit double precision (a-h,o-z)

      include 'types_tda_cp.inc'
      include 'input_tda_cp.inc'

      type(level_typ),dimension( * ) :: levne,levpe
      type(ph_typ),dimension( * ) :: iphne,iphpe
      
      
!*     generating p-h combinations for neutrons 
 8    format(10i4)      
      i=0
      do ip=ipnmn,ipnmx

        iparp=(-1)**levne(ip)%l
        
        do ih=ihnmn,ihnmx
          iparh=(-1)**levne(ih)%l
          ipart=iparh*iparp
          if ((levne(ip)%n-levne(ih)%n).gt.nhomx) goto 6
          if (ipart.ne.ipar) goto 6
          jjmn=iabs(levne(ip)%j-levne(ih)%j)/2
          jjmx=(levne(ip)%j+levne(ih)%j)/2
          if (ijj.lt.jjmn) goto 6
          if (ijj.gt.jjmx) goto 6
          i=i+1
          iphne(i)%par=ip
          iphne(i)%hol=ih
  6       continue         
        end do
        
      end do
      idphn=i

!*     generating p-h combinations for protons 
      
      i=0 
      do ip=ippmn,ippmx

        iparp=(-1)**levpe(ip)%l
        do ih=ihpmn,ihpmx
          iparh=(-1)**levpe(ih)%l
          ipart=iparh*iparp
          if ((levpe(ip)%n-levpe(ih)%n).gt.nhomx) goto 7
          if (ipart.ne.ipar) goto 7
          jjmn=iabs(levpe(ip)%j-levpe(ih)%j)/2
          jjmx=(levpe(ip)%j+levpe(ih)%j)/2
          if (ijj.lt.jjmn) goto 7
          if (ijj.gt.jjmx) goto 7

          i=i+1 
          iphpe(i)%par=ip
          iphpe(i)%hol=ih
 7        continue         
        end do
        
      end do
      idphp=i

      return
      end 

!*********************************************************************
      
      

      
