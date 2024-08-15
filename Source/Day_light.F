       Subroutine Day_light(T,Start_T,sflag)

c Sets a particle-specific flag for day or night time.

        use param

        implicit none

        Integer, INTENT(IN) :: T
        Real, INTENT(IN) :: Start_T
        Integer, INTENT(INOUT) :: sflag

        Integer, SAVE :: sunrise, sunset
        Integer :: actualtime, day, Tlx,Tly,I
        Real :: P

c Calculate positions in degrees

        If (curvi.eq.0) then
          Do I=1,Npart
            Xdeg(i)=((Xelev(i)-1.)*dch) + chst + (dch/2.)
            if (Xdeg(i).gt.180) Xdeg(i)=Xdeg(i)-360.
            Ydeg(i)=fist - (dfi/2.) - ((Yelev(i)-1.)*dfi)
          Enddo
        Elseif (curvi.eq.1) then
          Do I=1,Npart
            TLx=INT(Xelev(i))
            TLy=INT(Yelev(i))
            If (pom) then
              Xdeg(i)=Lons(TLx,TLy) + (Xelev(i)-TLx)*dchcu(Tlx,Tly)
     &                  /Deg2Rad
              Ydeg(i)=Lats(Tlx,Tly) - (Yelev(i)-Tly)*dfic(Tlx,Tly)
     &                  /Deg2Rad
            Else
              Xdeg(i)=Lons(TLx,TLy) + (Xelev(i)-TLx)*dchcu(Tlx,Tly)
              Ydeg(i)=Lats(Tlx,Tly) - (Yelev(i)-Tly)*dfic(Tlx,Tly)
            Endif
          Enddo
        Elseif (curvi.eq.2) then
          Do I=1,Npart
            Xdeg(i)=((Xelev(i)-1.)*dch) + chst
            if (Xdeg(i).gt.180) Xdeg(i)=Xdeg(i)-360.
            Tly=INT(Yelev(i))
            Ydeg(i)=lat(Tly) - (Yelev(i)-Tly)*dlat(Tly)
          Enddo
        Endif

c Calculate day length characteristics, just once per day, based on GMT

        If (mod(INT(Start_T+T),86400).eq.0.or.sflag.eq.0) then
         day=(Start_T+T)/86400
         if (day.lt.0) then
          day=day+365
         endif
         P=ASIN(0.39795*COS(0.2163108+2*
     &     ATAN(0.9671396*TAN(0.00860*(day-186.0)))))
         Do i=1,Npart
          nhours(i)=REAL(24.-(24./Pi)*ACOS((SIN(0.0145438)+
     & SIN(Pi*Ydeg(i)/180.)*SIN(P))/(COS(Pi*Ydeg(i)/180.)*COS(P))))
          sunrise=(12.*3600.)-(nhours(i)*3600./2.)
          sunset=(12.*3600.)+(nhours(i)*3600./2.)
c          write(52,*) 'In loop, i, nhours :',i,nhours(i)
         Enddo
        Endif

c Loop around particles
c Actualtime is the local time at the particle location

        Do i=1,Npart
         actualtime=Start_T+T+3600*NINT((Xdeg(i)/15.))
         If (MOD(actualtime,86400).gt.sunrise.and.
     &       MOD(actualtime,86400).lt.sunset) then
          light(i)='y'
         Else
          light(i)='n'
         Endif
        Enddo

        sflag=1

        Return

        End subroutine

