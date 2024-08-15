        Subroutine Convert

        use param

        implicit none

        Integer :: rguess,cguess,i,j,l
        Real :: refg,g,gx,gy
        Integer :: intc(jm)

c Initialise intc
        forall (j=1:jm) intc(j)=0

        If (curvi.eq.0) then
          do i=1,Npart
            if (circumpolar.and.xdeg(i).le.chst) then
             Xelev(i)=(xdeg(i)-(chst-360.)-(dch/2.))*ndch + 1.
            else
             Xelev(i)=(xdeg(i)-chst-(dch/2.))*ndch + 1.
            endif
            Yelev(i)=(fist-(dfi/2.)-ydeg(i))*ndfi + 1.
            if (polcoms) then
c U and V-pts are at SW corner of grid box in POLCOMS
              Xu(i)=Xelev(i) + 0.5
              Yu(i)=Yelev(i) - 0.5
              Xv(i)=Xelev(i) + 0.5
              Yv(i)=Yelev(i) - 0.5
            elseif (occam) then
c U and V-pts are at NE corner of grid box in OCCAM
              Xu(i)=Xelev(i) - 0.5
              Yu(i)=Yelev(i) + 0.5
              Xv(i)=Xelev(i) - 0.5
              Yv(i)=Yelev(i) + 0.5
            elseif (nemo) then
c U(V)-pt is on East(North) face of grid box in NEMO
              Xu(i)=Xelev(i) - 0.5
              Yu(i)=Yelev(i)
              Xv(i)=Xelev(i)
              Yv(i)=Yelev(i) + 0.5
            else
c U(V)-pt is on West(South) face of grid box in POM
              Xu(i)=Xelev(i) + 0.5
              Yu(i)=Yelev(i)
              Xv(i)=Xelev(i)
              Yv(i)=Yelev(i) - 0.5
            endif
          enddo
          print *,'Curvi = ',curvi
        Elseif (curvi.eq.1) then
         if (pom) then
          do i=1,Npart
            rguess=0
            cguess=0
            g=0.
c First guess at bottom right corner of grid box containing particle
            do l=1,jm
              do j=1,im
                if (lons(j,l).ge.Xdeg(i).and.lats(j,l).le.Ydeg(i)) then
                  rguess=l
                  cguess=j
                  goto 50
                endif
              enddo
            enddo
  50        l=rguess
            j=cguess
c Now test that have correct bottom right corner
            If (lons(j,l).gt.lons(j,l-1)) then
              refg=(lats(j,l-1)-lats(j,l))/(lons(j,l-1)-lons(j,l))
              If (lons(j,l-1).ne.Xdeg(i)) then
                g=(lats(j,l-1)-Ydeg(i))/(lons(j,l-1)-Xdeg(i))
                If (g.gt.refg) j=j+1
              Endif
            Elseif (lons(j-1,l).lt.lons(j-1,l-1)) then
              refg=(lats(j-1,l-1)-lats(j-1,l))/(lons(j-1,l-1)-
     &              lons(j-1,l))
              If (lons(j-1,l).ne.Xdeg(i)) then
                g=(Ydeg(i)-lats(j-1,l))/(Xdeg(i)-lons(j-1,l))
              Endif
              If (g.gt.refg.or.lons(j-1,l).eq.Xdeg(i)) j=j-1
            Endif
c Now calculate position in terms of grid elements
            Yelev(i) = l - (Ydeg(i)-lats(j,l))*Deg2Rad/dfic(j,l)
            Xelev(i) = j - (lons(j,l)-Xdeg(i))*Deg2Rad/dchcu(j,l)

c On regular space co-ordinates of transformed grid, u and v points are
c just half a grid element distance from elevation points...
c NB U-pt on West face of grid box in POM
            Xu(i) = Xelev(i) + 0.5
            Yu(i) = Yelev(i)
            Xv(i) = Xelev(i)
            Yv(i) = Yelev(i) - 0.5
          Enddo

         elseif (nemo) then
          do i=1,Npart
c Option for circumpolar
c Assumes loops from -179.75 to +180
            if (circumpolar.and.Xdeg(i).le.chst) then
             Xdeg(i)=Xdeg(i)+360.
            elseif (circumpolar.and.Xdeg(i).gt.180) then
             Xdeg(i)=Xdeg(i)-360.
            Endif

            rguess=0
            cguess=0
c Find closest point in longitude along each line of latitude
            do l=1,jm
              gx=1000.
              do j=1,im
                if (abs(lons(j,l)-Xdeg(i)).lt.gx) then
                  intc(l)=j
                  gx=abs(lons(j,l)-Xdeg(i))
                endif
              enddo
            enddo

c Find closest point in latitude from line of closest longitudes
            gy=1000.
            do l=1,jm
              if (abs(lats(intc(l),l)-Ydeg(i)).lt.gy) then
                  cguess=intc(l)
                  rguess=l
                  gy=abs(lats(intc(l),l)-Ydeg(i))
              endif
            enddo

c Now calculate position in terms of grid elements
            l=rguess
            j=cguess
            Xelev(i) = j - ((lons(j,l)-Xdeg(i))/dchcu(j,l))
            Yelev(i) = l - ((Ydeg(i)-lats(j,l))/dfic(j,l))

c On regular space co-ordinates of transformed grid, u and v points are
c just half a grid element distance from elevation points...
c NB U/V-pt on East/North face of grid box in NEMO
            Xu(i) = Xelev(i) - 0.5
            Yu(i) = Yelev(i)
            Xv(i) = Xelev(i)
            Yv(i) = Yelev(i) + 0.5
          enddo
         endif
        Elseif (curvi.eq.2) then
c Option for NEMO - regular in longitude, not in latitude
c Note, chst is top left t-point longitude!
c Assumes loops from -179.75 to +180
          Do i=1,Npart
            if (circumpolar.and.xdeg(i).le.chst) then
             Xelev(i)=(xdeg(i)-(chst-360.))*ndch + 1.
            elseif (circumpolar.and.xdeg(i).gt.180) then
             Xelev(i)=((xdeg(i)-360.)-chst)*ndch + 1.
            else
             Xelev(i)=(xdeg(i)-chst)*ndch + 1.
            endif
            rguess=0
            Do j=1,jm
c Find t-point immediately north of particle
              if (lat(j).gt.ydeg(i)) then
                rguess=j
              endif
            Enddo
            Yelev(i)=rguess-((ydeg(i)-lat(rguess))/(dlat(rguess)))
c On regular space co-ordinates of transformed grid, u and v points are
c just half a grid element distance from elevation points...
c NB U/V-pt on East/North face of grid box in NEMO
            Xu(i) = Xelev(i) - 0.5
            Yu(i) = Yelev(i)
            Xv(i) = Xelev(i)
            Yv(i) = Yelev(i) + 0.5

          Enddo
        Endif

        Return

        End subroutine

