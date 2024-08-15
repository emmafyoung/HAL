        Subroutine Output_netcdf(T,nrec,x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid,
     &endflag)

        use param
        use netcdf

        implicit none

        Integer, INTENT(IN) :: T,nrec,x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid
        Integer, INTENT(INOUT) :: endflag

        Integer :: unito,status
        Integer :: I
        Integer :: Tlx,Tly

c Data dimensions
        Integer, parameter :: ndims=2

c Start and count arrays for writing to netCDF output file
        Integer :: startp(ndims), countp(ndims)

c Program variables to hold the data to be written out
        Real :: Xout(Npart),Yout(Npart),Depout(Npart),Wtout(Npart)
        Real :: Swout(Npart),Xgout(Npart),Ygout(Npart),Dflag(Npart)
        Real :: Tavout(Npart),Savout(Npart)

        unito=luoutput(1)

        countp=(/ Npart, 1 /)
        startp=(/ 1, nrec /)

c Bypass this bit if first call to output initial positions
        If (T.gt.0) then
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
        Endif

c       Fill output arrays

        Xout=Xdeg(1:Npart)
        Yout=Ydeg(1:Npart)
        Depout=Deppart(1:Npart)
        Wtout=Partwt(1:Npart)
        Swout=Part_sw(1:Npart)
        Xgout=Xelev(1:Npart)
        Ygout=Yelev(1:Npart)
        Dflag=bathy_flag(1:Npart)
        If (TSopt) then
         Tavout=Tempav(1:Npart)
         Savout=Salinav(1:Npart)
        Endif

c First write time stamp
        status=nf90_put_var(unito,t_varid,T,start=(/ nrec /))
        call Check(status,endflag,'time stamp write')

c Now write particle information
        status=nf90_put_var(unito,x_varid,Xout,start=startp,
     &count=countp)
        call Check(status,endflag,'longitude write')
        status=nf90_put_var(unito,y_varid,Yout,start=startp,
     &count=countp)
        call Check(status,endflag,'latitude write')
        status=nf90_put_var(unito,d_varid,Depout,start=startp,
     &count=countp)
        call Check(status,endflag,'depth write')
        status=nf90_put_var(unito,w_varid,Wtout,start=startp,
     &count=countp)
        call Check(status,endflag,'weight write')
        status=nf90_put_var(unito,s_varid,Swout,start=startp,
     &count=countp)
        call Check(status,endflag,'switch write')
        status=nf90_put_var(unito,f_varid,Dflag,start=startp,
     &count=countp)
        call Check(status,endflag,'depth flag write')
        status=nf90_put_var(unito,xg_varid,Xgout,start=startp,
     &count=countp)
        call Check(status,endflag,'xelev write')
        status=nf90_put_var(unito,yg_varid,Ygout,start=startp,
     &count=countp)
        call Check(status,endflag,'yelev write')

        if (TSopt) then
         status=nf90_put_var(unito,tav_varid,Tavout,start=startp,
     &count=countp)
         call Check(status,endflag,'Tempav write')
         status=nf90_put_var(unito,sav_varid,Savout,start=startp,
     &count=countp)
         call Check(status,endflag,'Salinav write')
        endif

        Return

        End subroutine

