        Real Function Ran1(idum)

c Function to do a better job of random number generation than the
c system call 'Ran'. Taken from Numerical Recipes vol 2, p.271.

        implicit none

        Integer, INTENT(INOUT) :: idum

        Integer ia,im,iq,ir,ntab,ndiv
        Real am,eps,rnmx
        Parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     &             ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
        Integer j,k,iv(ntab),iy
        Save iv,iy
        Data iv /ntab*0/, iy /0/

        If (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
          enddo
          iy=iv(1)
        Endif
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum.lt.0) idum=idum+im
        j=1+iy/ndiv
        iy=iv(j)
        iv(j)=idum
        ran1=min(am*iy,rnmx)

        return

        end function Ran1

