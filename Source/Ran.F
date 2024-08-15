        Function Ran(idum)

c "Minimal" random number generator of Park and Miller combined with a
c Marsaglia shift sequence. Returns a uniform random deviate between 
c 0.0 and 1.0 (exclusive of the endpoint values). This fully portable, 
c scalar generator has the “traditional” (not Fortran 90) calling
c sequence with a random deviate as the returned function value: call with
c idum a negative integer to initialize; thereafter, do not alter idum 
c except to reinitialize. The period of this generator is about 3.1 × 1018.

        implicit none

        Integer, INTENT(INOUT) :: idum
        Real :: Ran

        Integer, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
        Real, SAVE :: am
        Integer, SAVE :: ix=-1,iy=-1,k

        if (idum <= 0 .or. iy < 0) then 
c Initialize
         am=nearest(1.0,-1.0)/IM
         iy=ior(ieor(888889999,abs(idum)),1)
         ix=ieor(777755555,abs(idum))
c Set idum positive
         idum=abs(idum)+1 
        end if

c Marsaglia shift sequence with period 232 − 1
        ix=ieor(ix,ishft(ix,13)) 
        ix=ieor(ix,ishft(ix,-17))
        ix=ieor(ix,ishft(ix,5))
        k=iy/IQ 
c Park-Miller sequence by Schrage’s method, period 231 − 2
        iy=IA*(iy-k*IQ)-IR*k
        if (iy < 0) iy=iy+IM
         ran=am*ior(iand(IM,ieor(ix,iy)),1) 
c Combine the two generators with masking to ensure nonzero value

        end Function Ran

