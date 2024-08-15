        Subroutine Mod_param (Outdt,Iseed2,Nrel_p,Start_T,Deppar,
     &                        endflag,slength,res_sw)

c Routine to read model parameterisation input file model_spec.dat, including
c information on model time step, source flow fields etc.

        use param

        implicit none

        Integer, INTENT(OUT) :: Outdt, Iseed2, Nrel_p, slength
        Integer, INTENT(INOUT) :: endflag
        Real, INTENT(OUT) :: Deppar, Start_T
        Character, INTENT(OUT) :: res_sw

        Integer Start_d
        Real Start_h
        Integer :: iunit

        iunit=luinput(1)

c Output frequency, diffusion coefficients and random number seed

        Read (iunit,5) Outdt
   5    Format (I5)
        Read (iunit,6) Ah
   6    Format (F7.3)
        Read (iunit,7) Az
   7    Format (F9.5)
        Read (iunit,5) Iseed2

c Model and forcing timesteps (secs) and start time of first velocity 
c file (Julian days). Also length of source files for warm start

        Read (iunit,10) Dt
  10    Format (I4)
        Read (iunit,11) Dt2
  11    Format (I7)
        print *,'Flows timestep: ',Dt2

        if (mod(int(Dt2),Dt).ne.0) then
          print *,'Timestep must be a divisor of hydrodynamic model 
     &                  output'
          print *,'Dt2: ',Dt2
          print *,'Dt: ',Dt
          endflag=-1
          goto 999
        endif

        Read (iunit,15) Start_d, Start_h
  15    Format (I4,F8.4)
        Start_T=((Start_d-1)*86400.)+(Start_h*3600.)

        Read (iunit,20) slength
  20    Format (I3)

c Type of velocities - depth-mean (y) or specified depth (n) and option
c to stop run at the end of the source data file (for multiple runs with
c restart)

        Read (iunit,60) veltyp
        Read (iunit,60) res_sw

c Specified start positions for each particle (y) or random elliptical
c distribution about one location (n)

        Read (iunit,60) Named
   60   Format (A1)

c All released at the same time (n) or staggered input - user specified
c time (y)

        Read (iunit,60) Namet

c All released at the same depth (n) or individually specified depth (y)

        Read (iunit,60) Nameh

c N.B. For staggered releases (namet = y), ellipse formulation is not 
c      allowed, so named must be set = y

        If ((Namet.eq.'y'.or.Namet.eq.'Y').and.(Named.eq.'n'.or.
     1       Named.eq.'N')) then
          print *,'For staggered releases, Named must be set = y'
          stop
        Endif

c N.B. For particle-specific depths (nameh = y), ellipse formulation is
c      no allowed, so named must be set = y

        If ((Nameh.eq.'y'.or.Nameh.eq.'Y').and.(Named.eq.'n'.or.
     1       Named.eq.'N')) then
          print *,'For particle-specific depths, Named must be set = y'
          stop
        Endif

c Input start depths for when Nameh = n (not used if Nameh = y)
        Read (iunit,70) Deppar
   70   Format (F8.3)

c Input number of release periods

        Read (iunit,5) Nrel_p

 999    return

        End subroutine

