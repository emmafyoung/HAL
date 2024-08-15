        Subroutine Calc_dens(i,rho)

c This subroutine calculates the local water density at the particle
c position using the same formulation as the Princeton Ocean Model
c The effect of pressure is neglected i.e. it is assumed that egg
c density and surrounding water will react to pressure in the same way

        use param

        implicit none

        Integer, INTENT(IN) :: i
        Real, INTENT(OUT) :: rho

c Declare local variables

        Real :: tr,tr2,tr3,tr4,sr,cr,p,grav
        Data grav/9.806/

        tr=temper(i)
        sr=salin(i)
        tr2=tr*tr
        tr3=tr2*tr
        tr4=tr3*tr

        rho = 999.842594  + 6.793952E-2*tr
     $        - 9.095290E-3*tr2 + 1.001685E-4*tr3
     $        - 1.120083E-6*tr4 + 6.536332E-9*tr4*tr

        rho = rho + (0.824493 - 4.0899E-3*tr
     $       + 7.6438E-5*tr2 - 8.2467E-7*tr3
     $       + 5.3875E-9*tr4) * sr
     $       + (-5.72466E-3 + 1.0227E-4*tr
     $       - 1.6546E-6*tr2) * ABS(sr)**1.5
     $       + 4.8314E-4 * sr*sr


        Return

        End subroutine

