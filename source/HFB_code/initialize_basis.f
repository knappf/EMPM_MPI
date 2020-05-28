      subroutine initialize_basis

       USE technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       allocate(levp(id),levn(id))
       allocate(lhfp(id),lhfn(id))

       im = 0

       do N = 0,noscmax
        do l = N,0,-1
         if(mod(N,2).eq.mod(l,2)) then
          nn = (N - l)/2
           do jj = 2*l+1, 2*l-1, -2
            if(jj.gt.0) then
             im = im + 1
             levp(im)%index = im
             levn(im)%index = im
             levp(im)%ipar = (-1)**l
             levn(im)%ipar = (-1)**l
             levp(im)%N = N
             levn(im)%N = N
             levp(im)%nn = nn
             levn(im)%nn = nn
             levp(im)%l = l
             levn(im)%l = l
             levp(im)%j2 = jj
             levn(im)%j2 = jj
             levp(im)%spenrg = hbarom * (dble(N)+1.5d0)
             levn(im)%spenrg = hbarom * (dble(N)+1.5d0)
             levp(im)%ei = levp(im)%spenrg
             levn(im)%ei = levn(im)%spenrg
             levp(im)%qei = levp(im)%spenrg
             levn(im)%qei = levn(im)%spenrg
            endif
           end do
         endif
        end do
       end do

      return
      end
