       Subroutine Input_bio (Nrel_p,Start_T,Rstart_T,Rend_T,Deppar)

c Routine to read particle information from input file particle_spec.dat

        use param

        implicit none

        Integer, INTENT(IN) :: Nrel_p
        Real, INTENT(IN) :: Start_T, Deppar
        Real, INTENT(OUT) :: Rstart_T, Rend_T

        Real :: Pstart_h, Pend_h
        Real :: th, r, xdsh, ydsh, xtemp, ytemp
        Real :: clatd, clatm, clond, clonm
        Real :: strch,ran1

        Integer :: Pstart_d, Pend_d, Npart_tot
        Integer :: Np1, Np2, dum, T, idum,i,l
        Integer :: uniti,uniti1,iseed2

        T=0
        Rstart_T=31536000.
        Rend_T=0.
        uniti=luinput(2)
        uniti1=luinput(1)
        iseed2=5

c Loop round all releases

        Npart_tot=0
        Do l=1,Nrel_p

c Input number of particles

         Read (uniti,*) dum
         Read (uniti,*) Npart
         Np1=Npart_tot+1
         Np2=Npart_tot+Npart

c Read position and times of release and removal of first particle
c If named=n this is the centre of release site (lat,long)
c If namet=n these are the release and end times for all particles in 
c this release period

         Read (uniti,*) Clatd,Clatm
         Read (uniti,*) Clond,Clonm
         If (Clond.lt.0.) Clonm=Clonm*(-1.)
         If (Clatd.lt.0.) Clatm=Clatm*(-1.)
         Read (uniti,*) Pstart_d, Pstart_h
         Read (uniti,*) Pend_d, Pend_h
         Pstart_T(Np1)=((Pstart_d-1)*86400.)+(Pstart_h*3600.)-Start_T
         Pend_T(Np1)=((Pend_d-1)*86400.)+(Pend_h*3600.)-Start_T
         If (Namet.eq.'n'.or.Namet.eq.'N') then
          Do i=Np1+1,Np2
            Pstart_T(i)=Pstart_T(Np1)
            Pend_T(i)=Pend_T(Np1)
          Enddo
         Endif
c If Nameh='n', set all start depths to Deppar
         if (Nameh.eq.'y'.or.Nameh.eq.'Y') then
          Read (uniti,*) Deppart(Np1)
         else
          Do i=Np1,Np2
            Deppart(i)=Deppar
          Enddo
         endif

c Particle weighting
         Read (uniti,*) Partwt(Np1)

c If only one particle, set it at this location

         If (Npart.eq.1) then
          Xdeg(Np1)=Clond+(Clonm/60.)
          Ydeg(Np1)=Clatd+(Clatm/60.)

c If specifying all particle locations, read all sites

         Else If (Named.eq.'Y'.or.Named.eq.'y') then
          Xdeg(Np1)=Clond+(Clonm/60.)
          Ydeg(Np1)=Clatd+(Clatm/60.)
          Do I=Np1+1,Np2
            Read (uniti,*) Clatd,Clatm
            Read (uniti,*) Clond,Clonm
            If (Clond.lt.0.) Clonm=Clonm*(-1.)
            If (Clatd.lt.0.) Clatm=Clatm*(-1.)
            Xdeg(I)=Clond+(Clonm/60.)
            Ydeg(I)=Clatd+(Clatm/60.)
            If (Namet.eq.'Y'.or.Namet.eq.'y') then
              Read (uniti,*) Pstart_d, Pstart_h
              Read (uniti,*) Pend_d, Pend_h
              Pstart_T(I)=((Pstart_d-1)*86400.)+(Pstart_h*3600.)-Start_T
              Pend_T(I)=((Pend_d-1)*86400.)+(Pend_h*3600.)-Start_T
            Endif
            If (Nameh.eq.'y'.or.Nameh.eq.'Y') Read (uniti,*) Deppart(i)
c Read particle weightings
            Read (uniti,*) Partwt(I)
          Enddo

c Else use ellipses to calculate particle locations.

         Else

          If (dum.eq.1) then
            Read (uniti1,70) aell,bell,dell
   70       Format (10F6.2)
            aell=aell*1000.
            bell=bell*1000.
          Endif

c Calculate a distribution of points about this site

          do i=Np1,Np2
           th=Ran(iseed)*2.*Pi
           r=Ran(iseed)
           xdsh=r*aell*cos(th)
           ydsh=r*bell*sin(th)
           xtemp=xdsh*cos(dell*Deg2Rad) - ydsh*sin(dell*Deg2Rad)
           ytemp=xdsh*sin(dell*Deg2Rad) + ydsh*cos(dell*Deg2Rad)
           ydeg(i)=Clatd+(Clatm/60.)+(ytemp/1.11e5)
           strch=1./cos(ydeg(i)*Deg2Rad)
           xdeg(i)=Clond+(Clonm/60.)+(strch*xtemp/1.11e5)
          enddo

c Assign particle weighting of particle 1 to all 
c particles in the ellipse

          Do i=Np1+1,Np2
           Partwt(i) = Partwt(Np1)
          Enddo

         Endif

         Npart_tot = Npart_tot + Npart

        Enddo

        Npart=Npart_tot

c Ensure particle start and end times are multiples of Dt
c Also store start depths

        Do I=1,Npart
          Pstart_T(I)=NINT(Pstart_T(I)/Dt) * Dt
          Pend_T(I)=NINT(Pend_T(I)/Dt) * Dt
          Rstart_T=MIN(Rstart_T,Pstart_T(I))
          Rend_T=MAX(Rend_T,Pend_T(I))
          Deppart_st(I)=Deppart(I)
        Enddo

        close (uniti1)
        close (uniti)

        Return

        End subroutine

