        Subroutine Setup_output_netcdf(x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid,
     &compflag,endflag)

        use param
        use netcdf

        implicit none

        Integer, INTENT(INOUT) :: compflag,endflag
        Integer, INTENT(OUT) :: x_varid,y_varid,d_varid,w_varid,s_varid,
     &t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid

        Integer :: unito

c Define the dimensions of the data
        Integer, parameter :: ndims=2
        Character (len=*), parameter :: part_name = "particle"
        Character (len=*), parameter :: rec_name = "timestep"
        Integer :: part_dimid, rec_dimid

c Create ten netCDF variables for particle characteristics and a time
c stamp
        Character (len=*), parameter :: x_name="longitude"
        Character (len=*), parameter :: y_name="latitude"
        Character (len=*), parameter :: d_name="depth"
        Character (len=*), parameter :: w_name="weight"
        Character (len=*), parameter :: s_name="switch"
        Character (len=*), parameter :: t_name="time"
        Character (len=*), parameter :: xg_name="xelev"
        Character (len=*), parameter :: yg_name="yelev"
        Character (len=*), parameter :: f_name="depth flag"
        Character (len=*), parameter :: tav_name="temperature"
        Character (len=*), parameter :: sav_name="salinity"
        Integer :: dimids(ndims)

c Create units attributes
        Character (len=*), parameter :: units="units"
        Character (len=*), parameter :: x_units="degrees_east"
        Character (len=*), parameter :: y_units="degrees_north"
        Character (len=*), parameter :: d_units="meters_below_sea_level"
        Character (len=*), parameter :: w_units="none"
        Character (len=*), parameter :: s_units="none"
        Character (len=*), parameter :: f_units="none"
        Character (len=*), parameter ::
     &t_units="seconds_from_model_start"
        Character (len=*), parameter :: xg_units="grid cells east"
        Character (len=*), parameter :: yg_units="grid cells north"
        Character (len=*), parameter :: tav_units="degrees C"
        Character (len=*), parameter :: sav_units="none"

        Integer :: status

        unito=luoutput(1)

c Define the dimensions

        status=nf90_def_dim(unito,part_name,Npart,part_dimid)
        call Check(status,endflag,'defining particle dimension')
        if (endflag.lt.0) goto 500
        status=nf90_def_dim(unito,rec_name,nf90_unlimited,
     &rec_dimid)
        call Check(status,endflag,'defining time dimension')
        if (endflag.lt.0) goto 500
        dimids=(/ part_dimid, rec_dimid /)

c Define the netCDF variables for the particle data and time

        status=nf90_def_var(unito,x_name,nf90_real,dimids,x_varid)
        call Check(status,endflag,'defining longitude variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,x_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,x_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining longitude compression')
         if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,y_name,nf90_real,dimids,y_varid)
        call Check(status,endflag,'defining latitude variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,y_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,y_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
        call Check(status,endflag,'defining latitude compression')
        if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,d_name,nf90_real,dimids,d_varid)
        call Check(status,endflag,'defining depth variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,d_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,d_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining depth compression')
         if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,w_name,nf90_real,dimids,w_varid)
        call Check(status,endflag,'defining weight variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,w_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,w_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining weight compression')
         if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,s_name,nf90_real,dimids,s_varid)
        call Check(status,endflag,'defining switch variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,s_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,s_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining switch compression')
         if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,t_name,nf90_int,rec_dimid,t_varid)
        call Check(status,endflag,'defining time variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,t_varid,nf90_chunked,
     &                               [1])
         status=nf90_def_var_deflate(unito,t_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining time compression')
         if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,xg_name,nf90_real,dimids,xg_varid)
        call Check(status,endflag,'defining xelev variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,xg_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,xg_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining xelev compression')
         if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,yg_name,nf90_real,dimids,yg_varid)
        call Check(status,endflag,'defining yelev variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,yg_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,yg_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining yelev compression')
         if (endflag.lt.0) goto 500
        endif

        status=nf90_def_var(unito,f_name,nf90_real,dimids,f_varid)
        call Check(status,endflag,'defining depth flag variable')
        if (endflag.lt.0) goto 500
        if (compflag.eq.1) then
         status=nf90_def_var_chunking(unito,f_varid,nf90_chunked,
     &                               [Npart,1])
         status=nf90_def_var_deflate(unito,f_varid,shuffle=1,deflate=1,
     &                              deflate_level=5)
         call Check(status,endflag,'defining depth flag compression')
         if (endflag.lt.0) goto 500
        endif

        if (TSopt) then
         status=nf90_def_var(unito,tav_name,nf90_real,dimids,tav_varid)
         call Check(status,endflag,'defining tempav variable')
         if (endflag.lt.0) goto 500
         if (compflag.eq.1) then
          status=nf90_def_var_chunking(unito,tav_varid,nf90_chunked,
     &                                [Npart,1])
          status=nf90_def_var_deflate(unito,tav_varid,shuffle=1,
     &                               deflate=1,
     &                               deflate_level=5)
          call Check(status,endflag,'defining tempav compression')
          if (endflag.lt.0) goto 500
         endif

         status=nf90_def_var(unito,sav_name,nf90_real,dimids,sav_varid)
         call Check(status,endflag,'defining salinav variable')
         if (endflag.lt.0) goto 500
         if (compflag.eq.1) then
          status=nf90_def_var_chunking(unito,sav_varid,nf90_chunked,
     &                                [Npart,1])
          status=nf90_def_var_deflate(unito,sav_varid,shuffle=1,
     &                               deflate=1,
     &                               deflate_level=5)
          call Check(status,endflag,'defining salinav compression')
          if (endflag.lt.0) goto 500
         endif
        endif

c Assign units to the variables
        status=nf90_put_att(unito,x_varid,units,x_units)
        call Check(status,endflag,'defining longitude units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,y_varid,units,y_units)
        call Check(status,endflag,'defining latitude units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,d_varid,units,d_units)
        call Check(status,endflag,'defining depth units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,w_varid,units,w_units)
        call Check(status,endflag,'defining weight units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,s_varid,units,s_units)
        call Check(status,endflag,'defining switch units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,t_varid,units,t_units)
        call Check(status,endflag,'defining time units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,xg_varid,units,xg_units)
        call Check(status,endflag,'defining xelev units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,yg_varid,units,yg_units)
        call Check(status,endflag,'defining yelev units')
        if (endflag.lt.0) goto 500
        status=nf90_put_att(unito,f_varid,units,f_units)
        call Check(status,endflag,'defining flag units')
        if (endflag.lt.0) goto 500

        if (TSopt) then
         status=nf90_put_att(unito,tav_varid,units,tav_units)
         call Check(status,endflag,'defining tempav units')
         if (endflag.lt.0) goto 500
         status=nf90_put_att(unito,sav_varid,units,sav_units)
         call Check(status,endflag,'defining salinav units')
         if (endflag.lt.0) goto 500
        endif

c End define mode
        status=nf90_enddef(unito)
        call Check(status,endflag,'ending define mode')
        if (endflag.lt.0) goto 500

  500   return

        End subroutine

