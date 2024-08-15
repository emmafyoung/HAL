        Program IBM

c HAL : Hydrodynamics-based Algorithm for Lagrangian simulations

c This program calculates egg/larval trajectories and environmental variables
c for a user-supplied number and location of particles.
c Adapted from PhD thesis (Young,1996) for any model grid.
c Note indices are relative to TOP LEFT!
c This version adapted for POLCOMS, POM, OCCAM and NEMO
c Only does temporal interpolation for points required by particles
c Can do restarts to bypass system memory problem - only for particle
c locations and flags
c Uses netcdf input and output files
c Type of grid set with curvi option in comblk_ibm.h file
c : 0 for regular grid
c : 1 for curvilinear grid, 
c : 2 for curvilinear only in latitude 
c *** For NEMO, use option 1 or 2 ***

c 1) User-supplied timestep
c 2) Staggered velocity input
c 3) Particles can be released at different times
c 4) Particles can be released at different depths
c 5) Successive current files can be used
c 6) Can read in model diffusivities for particle diffusion
c 7) Biological behaviour options; egg and larval growth and mortality,
c    sea ice-associated behaviour (krill), directed swimming (tidal,
c    bathymetric gradients, or related to currents), particle buoyancy (eggs
c    or oil droplets), DVM, diapause (copepods), beaching option

        use param
        use netcdf

        implicit none

c Declare external subroutines

        INTERFACE 
         SUBROUTINE Bio_Calc_Copepod(i,mod_jday)
          Integer, INTENT(IN) :: i, mod_jday
         END SUBROUTINE Bio_Calc_Copepod
         SUBROUTINE Bio_Calc_Egg(i,T)
          Integer, INTENT(IN) :: i, T
         END SUBROUTINE Bio_Calc_Egg
         SUBROUTINE Bio_Calc_Larvae(i,i1,i2,i3,i4,x1,y1,T)
          Integer, INTENT(IN) :: i,i1,i2,i3,i4,x1,y1,T
         END SUBROUTINE Bio_Calc_Larvae
         SUBROUTINE Bio_Calc_Krill(i)
          Integer, INTENT(IN) :: i
         END SUBROUTINE Bio_Calc_Krill
         SUBROUTINE Calc_dens(i,rho)
          Integer, INTENT(IN) :: i
          Real, INTENT(OUT) :: rho
         END SUBROUTINE Calc_dens
         SUBROUTINE Check(status,endflag,operation)
          Integer, INTENT(IN) :: status
          Integer, INTENT(INOUT) :: endflag
          Character(len=*),INTENT(IN) :: operation
         END SUBROUTINE Check
         SUBROUTINE Convert
         END SUBROUTINE Convert
         SUBROUTINE Day_light(T,Start_T,sflag)
          Integer, INTENT(IN) :: T
          Real, INTENT(IN) :: Start_T
          Integer, INTENT(INOUT) :: sflag
         END SUBROUTINE Day_light
         SUBROUTINE Diapause_depth(dia_distr)
          Character, INTENT(IN) :: dia_distr
         END SUBROUTINE Diapause_depth
         SUBROUTINE Filope(warmst,filenm,numnamesi)
          Character, INTENT(INOUT) :: warmst
          Character*50, INTENT(IN) :: filenm
          Integer, INTENT(OUT) :: numnamesi
         END SUBROUTINE Filope
         SUBROUTINE Input_bio(Nrel_p,Start_T,Rstart_T,Rend_T,Deppar)
          Integer, INTENT(IN) :: Nrel_p
          Real, INTENT(IN) :: Start_T, Deppar
          Real, INTENT(OUT) :: Rstart_T, Rend_T
         END SUBROUTINE Input_bio
         SUBROUTINE Input_bio2
         END SUBROUTINE Input_bio2
         SUBROUTINE Input_NEMO_netcdf(endflag)
          Integer, INTENT(INOUT) :: endflag
         END SUBROUTINE Input_NEMO_netcdf
         SUBROUTINE Input_phys(endflag)
          Integer, INTENT(INOUT) :: endflag
         END SUBROUTINE Input_phys
         SUBROUTINE Input_POM(endflag)
          Integer, INTENT(INOUT) :: endflag
         END SUBROUTINE Input_POM
         SUBROUTINE Mod_param(Outdt,Iseed2,Nrel_p,Start_T,Deppar,
     &                        endflag,slength,res_sw)
          Integer, INTENT(OUT) :: Outdt, Iseed2, Nrel_p, slength
          Integer, INTENT(INOUT) :: endflag
          Real, INTENT(OUT) :: Deppar, Start_T
          Character, INTENT(OUT) :: res_sw
         END SUBROUTINE Mod_param
         SUBROUTINE Output_netcdf(T,nrec,x_varid,y_varid,d_varid,
     &w_varid,s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,
     &sav_varid,endflag)
          Integer, INTENT(IN) :: T,nrec,x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid
          Integer, INTENT(INOUT) :: endflag
         END SUBROUTINE Output_netcdf
         REAL FUNCTION Ran(idum)
          Integer, INTENT(INOUT) :: idum
         END FUNCTION Ran
         SUBROUTINE Read_NEMO_netcdf(velcount,timecount,endflag,
     &switch,fcount,numdt,ncidu,ncidv,ncidw,ncidts,ncidi,varidu,varidv,
     &varidw,varidt,varids,varidi,varidiu,varidiv,numlons,numlats,
     &numlevs)
          Integer, INTENT(IN) :: timecount,switch
          Integer, INTENT(INOUT) :: velcount,endflag,fcount,numdt,ncidu,
     &ncidv,ncidw,ncidts,ncidi,varidu,varidv,varidw,varidt,varids,
     &varidi,varidiu,varidiv,numlons,numlats,numlevs
         END SUBROUTINE Read_NEMO_netcdf
         SUBROUTINE Read_OCCAM(velcount,timecount,endflag,switch)
          Integer, INTENT(IN) :: timecount,switch
          Integer, INTENT(INOUT) :: velcount,endflag
         END SUBROUTINE Read_OCCAM
         SUBROUTINE Read_POLCOMS(velcount,timecount,endflag,switch)
          Integer, INTENT(IN) :: timecount,switch
          Integer, INTENT(INOUT) :: velcount,endflag
         END SUBROUTINE Read_POLCOMS
         SUBROUTINE Read_POM(velcount,timecount,endflag,switch)
          Integer, INTENT(IN) :: timecount,switch
          Integer, INTENT(INOUT) :: velcount,endflag
         END SUBROUTINE Read_POM
         SUBROUTINE Scoordset
         END SUBROUTINE Scoordset
         SUBROUTINE Sigmaset
         END SUBROUTINE Sigmaset
         SUBROUTINE Setup_output_netcdf(x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid,
     &compflag,endflag)
          Integer, INTENT(INOUT) :: compflag,endflag
          Integer, INTENT(OUT) :: x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid
         END SUBROUTINE Setup_output_netcdf
        END INTERFACE

c Declare local variables

        Real :: Start_T
        Real :: Rstart_T, Rend_T
        Real :: th, r, xdsh, ydsh, xtemp, ytemp, strch
        Real :: Deppar, adj_res, end_res

        Integer :: T, Outdt, timecount, velcount, Nrel_p, Iseed2, I
        Integer :: numf,icn,icn2,itmp,endflag,switch,nfstart,ierror
        Integer :: slength,status,mod_jday,sflag,compflag

        Integer :: x_varid,y_varid,d_varid,w_varid,s_varid,t_varid
        Integer :: xg_varid,yg_varid,f_varid,tav_varid,sav_varid
        Integer :: nrec
        Integer :: iunit
        Integer :: fcount,numdt,ncidu,ncidv,ncidw
        Integer :: ncidts,ncidi
        Integer :: varidu,varidv,varidw
        Integer :: varidt,varids,varidi,varidiu,varidiv
        Integer :: numlons,numlats,numlevs,numnamesi

        Character :: warmst,res_sw,dia_distr
        Character*50 :: filenm

c Allocate arrays

        allocate(nemolev(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for nemolev***"
        allocate(nemowlev(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for nemowlev***"
        allocate(Drift1(Mx2,3,kb-1),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for drift1***"
        allocate(Drift2(Mx2,3,kb-1),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for drift2***"
        allocate(z(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for z***"
        allocate(zz(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for zz***"
        allocate(zzv(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for zzv***"
        allocate(U(Mx2u,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for U***"
        allocate(V(Mx2u,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for V***"
        allocate(W(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for W***"
        allocate(Tem(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Tem***"
        allocate(Sal(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Sal***"
        allocate(Ice(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ice***"
        allocate(Uice(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Uice***"
        allocate(Vice(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Vice***"
        allocate(PAh(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PAh***"
        allocate(PKh(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PKh***"
        allocate(U1(Mx2u,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for U1***"
        allocate(U2(Mx2u,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for U2***"
        allocate(V1(Mx2u,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for V1***"
        allocate(V2(Mx2u,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for V2***"
        allocate(W1(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for W1***"
        allocate(W2(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for W2***"
        allocate(T1(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for T1***"
        allocate(T2(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for T2***"
        allocate(Sa1(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Sa1***"
        allocate(Sa2(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Sa2***"
        allocate(Ice1(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ice1***"
        allocate(Ice2(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ice2***"
        allocate(Uice1(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Uice1***"
        allocate(Uice2(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Uice2***"
        allocate(Vice1(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Vice1***"
        allocate(Vice2(Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Vice2***"
        allocate(PAh1(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PAh***"
        allocate(PAh2(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PAh2***"
        allocate(PKh1(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PKh1***"
        allocate(PKh2(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PKh2***"

c Initialise variables

        timecount=0
        iseed=13
        T=0
        velcount=0
        velnum=1
        endflag=0
        switch=0
        warmst='n'
        adj_res=0.
        end_res=100000000.
        nrec=1
        slength=0
        res_sw='n'
        fcount=1
        sflag=0
        compflag=0

        Do I=1,Maxp
          Part_sw(I)=0
          dia_flag(I)=0
          bathy_flag(I)=0
          Tempav(I)=0.
          Salinav(I)=0.
          Tempav_bio(I)=0.
        Enddo

        print *,'Iopt, TSopt : ',Iopt,TSopt

c Only allow certain combinations of model choices

        if (TSopt.and.(.not.(polcoms.or.nemo))) then
          print *,'T and S option only available for POLCOMS or NEMO'
          stop
        endif

        if (Iopt.and.(.not.nemo)) then
          print *,'Sea ice option only available for NEMO'
          stop
        endif

        if (nemo.and.(curvi.eq.0)) then
          print *,'Choose curvi option 1 or 2 for NEMO'
          stop
        endif

c Get argument (filename) from command line

        call getarg (1,filenm)
        print *,'filename is : ',trim(filenm)

c*****START OF PROGRAM!!!*****

c Call subroutine to open input and output files

        Call Filope(warmst,filenm,numnamesi)

c Input run information

        call Mod_param (Outdt,Iseed2,Nrel_p,Start_T,Deppar,endflag,
     &slength,res_sw)
        If (endflag.lt.0) goto 2000

c Read particle information file

        Call Input_bio (Nrel_p,Start_T,Rstart_T,Rend_T,Deppar)
        If (Npart.ge.1000) compflag=1

c Set up the output netcdf file
        if (nemo) then
          Call Setup_output_netcdf(x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid,
     &compflag,endflag)
        endif
        If (endflag.lt.0) goto 2000

c Input physical and model setup information

        if (polcoms.or.occam) then
         Call Input_phys(endflag)
        elseif (nemo) then
          Call Input_NEMO_netcdf(endflag)
        else
          Call Input_POM(endflag)
        endif
        If (endflag.lt.0) goto 2000

c Calculate vertical levels for POLCOMS model

        if (polcoms) then
          if (slevel) then
            call Scoordset
          else
            call Sigmaset
          endif
        endif

c If warm start, overwrite position and flag info. from restart file

        if (warmst.eq.'y') then
c Assumes restart file is the last in the input file list
         iunit=luinput(numnamesi)
         do i=1,Npart
          read (iunit,*) Xdeg(i),Ydeg(i),Deppart(i),Partwt(i),Part_sw(i)
         enddo
        endif
        close (iunit)

c Possible correction for circumpolar NEMO
         if (nemo.and.circumpolar) then
          do i=1,Npart
           if (chst.gt.0.and.xdeg(i).lt.0) then
            xdeg(i)=xdeg(i)+360
           endif
          enddo
         endif

c Call a subroutine to convert from lat/long to grid co-ordinates

        Call Convert

c Write initial positions to output file - only if not warm start

        if (warmst.ne.'y') then
          Call Output_netcdf(T,nrec,x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid,
     &endflag)
          nrec=nrec+1
        endif
        If (endflag.lt.0) goto 2000

c Read biological parameter file

        Call Input_bio2

c Extract flow field file information
c Assumes this file is the 4th in the input file list

        read(luinput(4),*) numf
        if (numf.ne.0) then
         do i=1,numf
          read(luinput(4),'(a)') velfiles(i)
          print *,'Name of velfile ',i,' : ',velfiles(i)
         enddo
        endif
        read(luinput(4),*) nfstart

c Calculate time adjustment for warm start
c Also for adjusting time origin for larval runs subsequent to egg runs
c Length of velocity field included in ibm.dat (slength)

        adj_res=(nfstart-1)*slength*86400.
        if (res_sw.eq.'y') then
          end_res=(nfstart)*slength*86400.
        endif

        Close (luinput(4))

c *** MAIN LOOP *** 

  150   Continue

c Read forcing data from ocean model output file

        if (polcoms) then
          Call Read_POLCOMS(velcount,timecount,endflag,switch)
        elseif (occam) then
          Call Read_OCCAM(velcount,timecount,endflag,switch)
        elseif (nemo) then
          Call Read_NEMO_netcdf(velcount,timecount,endflag,
     &switch,fcount,numdt,ncidu,ncidv,ncidw,ncidts,ncidi,varidu,
     &varidv,varidw,varidt,varids,varidi,varidiu,varidiv,numlons,
     &numlats,numlevs)
        else
          Call Read_POM(velcount,timecount,endflag,switch)
        endif
        If (endflag.lt.0) goto 2000

        timecount=timecount+1
        T=timecount*Dt + adj_res
        mod_jday=INT((Start_T+T)/86400)+1


c Update particle switches and switch for input field interpolation
c The latter requires a consideration of the next timestep

        Do I=1,Npart
          if (T.ge.Pstart_T(I)) then
            if (Part_sw(I).ne.-2) then
              Part_sw(I)=1
            endif
            switch=1
          endif
          if (T+Dt.ge.Pstart_T(I)) switch=1
          If (T.ge.Pend_T(I)) Part_sw(I)=-1
        Enddo

        if (switch.eq.1.and.species.ne.'passive')
     &    Call Day_light(T,Start_T,sflag)

c Call main calculation routine

        If (T.ge.Rstart_T) then
          Call Maincalc(T,Iseed2,velcount,endflag,mod_jday,Outdt,
     &Rstart_T)
        Else
          Goto 150
        Endif

c Exit program if all particles are at their end time or an error was generated

        If (T.eq.Rend_T.or.T.eq.end_res) Goto 2000
        If (endflag.lt.0) goto 2000

c Write new larval positions and temperatures

        If ((T-NINT(Rstart_T)).gt.0.and.mod((T-NINT(Rstart_T)+Dt),Outdt)
     &       .eq.0) then
          Call Output_netcdf(T,nrec,x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid,
     &endflag)
          nrec=nrec+1
        Endif
        If (endflag.lt.0) goto 2000

        goto 150

 2000   Call Output_netcdf(T,nrec,x_varid,y_varid,d_varid,w_varid,
     &s_varid,t_varid,xg_varid,yg_varid,f_varid,tav_varid,sav_varid,
     &endflag)

c Create restart file
        do i=1,Npart
         write (luoutput(2),*) Xdeg(i),Ydeg(i),Deppart(i),Partwt(i),
     &          Part_sw(i)
        enddo
        close (luoutput(2))

c Close output netcdf file and NEMO input netcdf files
        status=nf90_close(luoutput(1))
        if (nemo) then
         status=nf90_close(ncidu)
         status=nf90_close(ncidv)
         status=nf90_close(ncidw)
         if (TSopt) then
          status=nf90_close(ncidts)
         endif
         if (Iopt) then
          status=nf90_close(ncidi)
         endif
        endif

c Deallocate arrays

        deallocate(nemolev)
        deallocate(nemowlev)
        deallocate(Drift1)
        deallocate(Drift2)
        deallocate(z)
        deallocate(zz)
        deallocate(zzv)
        deallocate(U)
        deallocate(V)
        deallocate(W)
        deallocate(Tem)
        deallocate(Sal)
        deallocate(Ice)
        deallocate(Uice)
        deallocate(Vice)
        deallocate(PKh)
        deallocate(PAh)
        deallocate(U1)
        deallocate(U2)
        deallocate(V1)
        deallocate(V2)
        deallocate(W1)
        deallocate(W2)
        deallocate(T1)
        deallocate(T2)
        deallocate(Sa1)
        deallocate(Sa2)
        deallocate(Ice1)
        deallocate(Ice2)
        deallocate(Uice1)
        deallocate(Uice2)
        deallocate(Vice1)
        deallocate(Vice2)
        deallocate(PKh1)
        deallocate(PKh2)
        deallocate(PAh1)
        deallocate(PAh2)

        end

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Maincalc(T,Iseed2,velcount,endflag,mod_jday,Outdt,
     &Rstart_T)

        use param

        implicit none

        Real*8 :: Xdgrid, Ydgrid, OldXelev(Maxp), OldYelev(Maxp),
     &          OldXu(Maxp), OldYu(Maxp), OldXv(Maxp), OldYv(Maxp)
        Real :: Plat, ran1,Rstart_T
        Real :: Hfull(Mx2)
        Real :: Newz, vsp1, vsp2, vsp3, vsp4
        Real :: pt, qt, alphat, betat, gammat, deltat, mdch, mdfi
        Real :: pu, qu, alphau, betau, gammau, deltau
        Real :: pv, qv, alphav, betav, gammav, deltav
        Real :: ptu, ptv, ptw, Xdist, Ydist, Zdist, pth, ptkh, ptah
        Real :: Theta, Dah, Daz
        Real :: layerdep, grav, visc, rho, layerdepb
        Real :: zt(kb),s1(kb),s2,vz(kb),dzt(kb)
        Real :: tsp1, tsp2, tsp3, tsp4, ssp1, ssp2, ssp3, ssp4
        Real :: ptt,qtt,alphatt,betatt,gammatt,deltatt

        Integer :: count, Iseed2, i, j, k, ierror, Outdt
        Integer :: flag, velcount, idum, ibox, flagi, iboxu, iboxv
        Integer :: T, Tflag, Back, Landflag, endflag
        Integer :: ilev, timerem, nlev
        Integer :: i1t, i2t, i3t, i4t, x1t, y1t
        Integer :: x2t, y2t, x3t, y3t, x4t, y4t
        Integer :: i1u, i2u, i3u, i4u, x1u, y1u
        Integer :: x2u, y2u, x3u, y3u, x4u, y4u
        Integer :: i1v, i2v, i3v, i4v, x1v, y1v
        Integer :: x2v, y2v, x3v, y3v, x4v, y4v
        Integer :: mod_jday
        Integer :: bflag(Maxp)

        Character :: btyp

        Save idum

        Data Tflag /0/
        Data grav/9.81/,visc/0.001401/

c Initialise variables

        write(6,*) 'In maincalc',T

        Tflag=Tflag+1
        If (Tflag.eq.1) Idum=Iseed2

        do 5 i=1,Mx2
         Hfull(i) = 0.0
    5 	continue

        do 10 j=1,nz
         if (polcoms) then
          i=iz(j)
          Hfull(i) = Hzeta(i)
         elseif (pom) then
          do 20 i=iz1(j),iz2(j)
            Hfull(i) = Hzeta(i) + H(i)
   20     continue
         endif
   10   continue
        if (occam.or.nemo) then
         do 15 i=1,Mx2
          Hfull(i) = H(i)
   15    continue
        endif

c Set temporary arrays to zero

        if (occam.or.polcoms.or.nemo) then
         W=0.0
         Tem=0.0
         Sal=0.0
         PKh=0.0
         PAh=0.0
         U=0.0
         V=0.0
         Ice=0.0
         Uice=0.0
         Vice=0.0
        endif

c Calculate drift terms

        if (Ah.lt.0.) Call Calc_drift(Hfull)

        do i=1,Npart
          OldXelev(i)=Xelev(i)
          OldYelev(i)=Yelev(i)
          OldXu(i)=Xu(i)
          OldYu(i)=Yu(i)
          OldXv(i)=Xv(i)
          OldYv(i)=Yv(i)
          bflag(i)=0
        enddo
        If (mod((T-NINT(Rstart_T)),Outdt).eq.0) then
          do i=1,Npart
            bathy_flag(i)=0
            if (TSopt) then
             Tempav(i)=0.
             Salinav(i)=0.
            endif
          enddo
        Endif

        do 150 i=1,npart

         If (Part_sw(i).ne.1) goto 150

         count=0
         Back=0
         newz=0.
         ptu=0.
         ptv=0.
         ptw=0.
         flagi=1

c Find 4 surrounding T-points

         x1t=INT(Xelev(i))
         y1t=INT(Yelev(i))
 
         Call Findpts(x1t,x2t,x3t,x4t,y1t,y2t,y3t,y4t,i1t,i2t,i3t,i4t)
         If ((circumpolar.and.y1t.le.1.or.y3t.ge.jm).or.
     &      ((.not.circumpolar).and.(x2t.ge.im.or.y1t.le.1.or
     &        .x1t.le.1.or.y3t.ge.jm))) then
           part_sw(i)=-1
           goto 150
         Endif

c Find weighted distances relative to T-points

         pt=Xelev(i)-x1t
         qt=Yelev(i)-y1t
         alphat=pt*qt
         betat=pt*(1-qt)
         gammat=(1-pt)*(1-qt)
         deltat=(1-pt)*qt
 
c Find 4 surrounding U-points

         x1u=INT(Xu(i))
         y1u=INT(Yu(i))

         Call Findpts(x1u,x2u,x3u,x4u,y1u,y2u,y3u,y4u,i1u,i2u,i3u,i4u)
         If ((circumpolar.and.y1u.le.1.or.y3u.ge.jm).or.
     &      ((.not.circumpolar).and.(x2u.ge.im.or.y1u.le.1.or
     &        .x1u.le.1.or.y3u.ge.jm))) then
           part_sw(i)=-1
           goto 150
         Endif

c Find weighted distances relative to U-points

          pu=Xu(i)-x1u
          qu=Yu(i)-y1u
          alphau=pu*qu
          betau=pu*(1-qu)
          gammau=(1-pu)*(1-qu)
          deltau=(1-pu)*qu

c Find 4 surrounding V-points

         x1v=INT(Xv(i))
         y1v=INT(Yv(i))

         Call Findpts(x1v,x2v,x3v,x4v,y1v,y2v,y3v,y4v,i1v,i2v,i3v,i4v)
         If ((circumpolar.and.y1v.le.1.or.y3v.ge.jm).or.
     &      ((.not.circumpolar).and.(x2v.ge.im.or.y1v.le.1.or
     &        .x1v.le.1.or.y3v.ge.jm))) then
           part_sw(i)=-1
           goto 150
         Endif

c Find weighted distances relative to V-points

         pv=Xv(i)-x1v
         qv=Yv(i)-y1v
         alphav=pv*qv
         betav=pv*(1-qv)
         gammav=(1-pv)*(1-qv)
         deltav=(1-pv)*qv

c Find depth at particle position (depth of cell)

         ibox=(IDNINT(Yelev(i))-1)*im + IDNINT(Xelev(i))
         pth=Hfull(ibox)

c If point has become a tidal flat point since the last calculation,
c leave it where it is.

         if (Hfull(ibox).lt.1E-20) goto 150

c Find cos(latitude) at particle position

         plat=((1-qt)*cspz(i1t) + qt*cspz(i3t))

c Find depth layer of particle

         ilev=1
         do 200 k=1,kb-1
          if (polcoms.or.pom) then
           layerdep=(-1.)*z(ibox,k)*Hfull(ibox)
          elseif (occam) then
           layerdep=hlev(k)
          elseif (nemo) then
           layerdep=nemowlev(ibox,k+1)
          endif
c Note should be k+1 for k=1 at surface
          if (deppart(i).gt.layerdep.and.layerdep.gt.0) ilev=k+1
  200    continue

c End if having trouble finding model layer
         if (nemo.and.nemolev(ibox,ilev).lt.0) then
          print *,'Particle below nemo sea bed!'
          print *,'Particle number : ',i
          print *,'Particle depth : ',deppart(i)
          print *,'i,j,ibox : ',IDNINT(Xelev(i)),IDNINT(Yelev(i)),ibox
          print *,'Hfull : ',Hfull(ibox)
          print *,'Vertical level and depth : ',ilev,layerdep
          endflag=-1
          goto 150
         endif

c Calculate diffusion terms at particle position
c First for vertical diffusion, Kh
c If Az < 0, use model fields, else fixed in space and time.

         if (Az.lt.0.) then

c Temporally interpolate diffusivities to local points
c POLCOMS only at present...
          if (polcoms) then
           If (PKh(i1t,2).eq.0.or.PKh(i2t,2).eq.0.or.PKh(i3t,2).eq.0.or.
     &        PKh(i4t,2).eq.0) then
            Call Tinterpvel(5,i1t,i2t,i3t,i4t,velcount)
           Endif
          Endif

          flag=5
          Landflag=0
          If (veltyp.eq.'N'.or.veltyp.eq.'n') then
           Call Calcvel1(i,flag,x1t,y1t,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                   Landflag)
          Else
           Call Calcvel2(flag,x1t,y1t,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                   Landflag)
          Endif
          If (Landflag.eq.0) then
           ptkh=(alphat*vsp4 + betat*vsp2 + gammat*vsp1 + deltat*vsp3)
          Else
           Do k=2,kb-1
            zt(k)=z(ibox,k)*Hfull(ibox)
            s1(k)=PKh(ibox,k)
           Enddo
           s1(1)=s1(2)
           If (veltyp.eq.'N'.or.veltyp.eq.'n') then
            Call Intpol(zt,s1,deppart(i),s2,kb-1)
           Else
            zt(kb)=-1.*Hfull(ibox)
            s2=0.
            Do k=1,kb-1
             dzt(k)=zt(k)-zt(k+1)
             vz(k)=s1(k)*dzt(k)
             s2=s2+vz(k)
            Enddo
            s2=s2/Hfull(ibox)
           Endif
           ptkh=s2
          Endif
         else
          ptkh=Az
         endif

c Now for horizontal diffusion, Ah
c If Ah < 0, use model fields, else fixed in space and time.

         if (Ah.lt.0.) then

c Temporally interpolate diffusivities to local points
c POLCOMS only at present...
          if (polcoms) then
           If (PAh(i1t,2).eq.0.or.PAh(i2t,2).eq.0.or.PAh(i3t,2).eq.0.or.
     &        PAh(i4t,2).eq.0) then
            Call Tinterpvel(6,i1t,i2t,i3t,i4t,velcount)
           Endif
          Endif

          flag=6
          Landflag=0
          If (veltyp.eq.'N'.or.veltyp.eq.'n') then
           Call Calcvel1(i,flag,x1t,y1t,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                   Landflag)
          Else
           Call Calcvel2(flag,x1t,y1t,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                   Landflag)
          Endif
          If (Landflag.eq.0) then
           ptah=(alphat*vsp4 + betat*vsp2 + gammat*vsp1 + deltat*vsp3)
          Else
           Do k=1,kb-1
            zt(k)=z(ibox,k)*Hfull(ibox)
            s1(k)=PAh(ibox,k)
           Enddo
           If (veltyp.eq.'N'.or.veltyp.eq.'n') then
            Call Intpol(zt,s1,deppart(i),s2,kb-1)
           Else
            zt(kb)=-1.*Hfull(ibox)
            s2=0.
            Do k=1,kb-1
             dzt(k)=zt(k)-zt(k+1)
             vz(k)=s1(k)*dzt(k)
             s2=s2+vz(k)
            Enddo
            s2=s2/Hfull(ibox)
           Endif
           ptah=s2
          Endif
         else
          ptah=Ah
         endif

c Update stage and biological parameter data
c This is now done before the calculation of new particle positions
c to allow for possibility of active horizontal swimming
c Also before ice calculation so performs DVM under ice
c Note, egg and larvae calculations involving growth and mortality are
c not done at every physics timestep unless biocalcfreq is
c set equal to Dt

         if (species.eq.'passive') goto 300
         timerem=mod((T-NINT(Pstart_T(i))),biocalcfreq(part_st(i)))
         if (timerem.eq.0) then
          if (stage_spec(part_st(i)).eq.'e')
     &     Call Bio_Calc_Egg(i,T)
          if (stage_spec(part_st(i)).eq.'l')
     &     Call Bio_Calc_Larvae(i,i1t,i2t,i3t,i4t,x1t,y1t,pth,T)
         endif
         if (stage_spec(part_st(i)).eq.'k')
     &    Call Bio_Calc_Krill(i)
         if (stage_spec(part_st(i)).eq.'c')
     &    Call Bio_Calc_Copepod(i,mod_jday)

  300    if (stage_spec(part_st(i)).eq.'n') then
          if (T.gt.(Pstart_T(i)+stagelen(part_st(i)))) 
     &     part_st(i)=part_st(i)+1
         endif
         if (part_st(i).gt.true_no_stages) then
          part_sw(i)=-1
          goto 150
         endif
         btyp=buoytyp(part_st(i))

c Find ice at particle position

         if (Iopt.and.nemo) then

          If (Ice(i1t).eq.0.or.Ice(i2t).eq.0.or.Ice(i3t).eq.0.or.
     &        Ice(i4t).eq.0) then
           Call Tinterpvel(8,i1t,i2t,i3t,i4t,velcount)
          Endif
          flag=4
          Landflag=0
          Call Calcvel4(i,flag,x1t,y1t,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                  Landflag,1)
          Icep(i)=(alphat*vsp4 + betat*vsp2 + gammat*vsp1 + deltat*vsp3)

c Reset particle depth if below maximum depth under ice
c Also moves particles if moving into an ice area

          if (Icep(i).ge.0.15.and.Deppart(i).gt.maxdi) then
           Deppart(i)=maxdi
           if (Tflag.eq.1)
     &      print *,'Moving particle ',i,' to maximum depth under ice'
          endif

c If Ice concentration >= 15%, calculate ice velocities
c Extra flag for krill performing DVM under ice
c Also set particle depth to base of seasonal sea ice (2m) for
c ice-associated larvae

          if (Icep(i).ge.0.15.and.((migrationYN.eq.'n').or.
     &        (migrationYN.eq.'y'.and.bwlarvae(i).gt.0.and.
     &         deppart(i).lt.5.))) then

           flagi=0
           deppart(i)=2.0
           if (Uice(i1u).eq.0.or.Uice(i2u).eq.0.or.Uice(i3u).eq.0.or.
     &         Uice(i4u).eq.0) then
            Call Tinterpvel(9,i1u,i2u,i3u,i4u,velcount)
           endif
           if (Vice(i1v).eq.0.or.Vice(i2v).eq.0.or.Vice(i3v).eq.0.or.
     &         Vice(i4v).eq.0) then
            Call Tinterpvel(10,i1v,i2v,i3v,i4v,velcount)
           endif
           flag=5
           Landflag=0
           Call Calcvel4(i,flag,x1u,y1u,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                   Landflag,1)
           ptu=(alphau*vsp4 + betau*vsp2 + gammau*vsp1 + deltau*vsp3)
           flag=6
           Landflag=0
           Call Calcvel4(i,flag,x1v,y1v,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                   Landflag,1)
           ptv=(alphav*vsp4 + betav*vsp2 + gammav*vsp1 + deltav*vsp3)

          endif

         endif

c Bypass ocean velocity calculations if Iopt set and ice concentration >= 15%
c But don't bypass if particles performing DVM under ice

         if (flagi.eq.1) then

c Reset initial particle depth if below maximum particle depth
          if (Tflag.eq.1.and.Deppart(i).gt.maxd) then
           Deppart(i)=maxd
           print *,'Moving particle ',i,' to maximum particle depth'
          endif

c *** U ***

c Temporally interpolate velocities to 4 surrounding U-points
c OCCAM, POLCOMS and NEMO only...

          if (occam.or.polcoms.or.nemo) then
           If (U(i1u,2).eq.0.or.U(i2u,2).eq.0.or.U(i3u,2).eq.0.or.
     &        U(i4u,2).eq.0) then
            Call Tinterpvel(1,i1u,i2u,i3u,i4u,velcount)
           Endif
          Endif

c Potentially find velocity at particle depth, else use depth-mean.
c For OCCAM just extract velocity of the depth slice.

          Landflag=0
          if (occam) then
           Call Calcvel3(i,1,x1u,y1u,vsp1,vsp2,vsp3,vsp4,Hfull,Landflag,
     &                   ilev)
          elseif (nemo) then
           Call Calcvel5(i,1,x1u,y1u,vsp1,vsp2,vsp3,vsp4,Hfull,Landflag)
          else
           if (polcoms) then
            flag=8
           else
            flag=1
           endif
           If (veltyp.eq.'N'.or.veltyp.eq.'n') then
            Call Calcvel1(i,flag,x1u,y1u,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                    Landflag)
           Else
            Call Calcvel2(flag,x1u,y1u,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                    Landflag)
           Endif
          endif

c Calculate velocity at particle position

          ptu=(alphau*vsp4 + betau*vsp2 + gammau*vsp1 + deltau*vsp3)

c *** V ***

c Temporally interpolate velocities to 4 surrounding V-points

          if (occam.or.polcoms.or.nemo) then
           If (V(i1v,2).eq.0.or.V(i2v,2).eq.0.or.V(i3v,2).eq.0.or.
     &       V(i4v,2).eq.0) then
            Call Tinterpvel(2,i1v,i2v,i3v,i4v,velcount)
           Endif
          endif

c Potentially find velocity at particle depth, else use depth-mean.

          Landflag=0
          if (occam) then
           call Calcvel3(i,2,x1v,y1v,vsp1,vsp2,vsp3,vsp4,Hfull,Landflag,
     &                   ilev)
          elseif (nemo) then
           call Calcvel5(i,2,x1v,y1v,vsp1,vsp2,vsp3,vsp4,Hfull,Landflag)
          else
           if (polcoms) then
            flag=9
           else
            flag=2
           endif
           If (veltyp.eq.'N'.or.veltyp.eq.'n') then
            Call Calcvel1(i,flag,x1v,y1v,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                    Landflag)
           Else
            Call Calcvel2(flag,x1v,y1v,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                    Landflag)
           Endif
          endif

c Calculate velocity at particle position

          ptv=(alphav*vsp4 + betav*vsp2 + gammav*vsp1 + deltav*vsp3)

         Endif

c End of possible bypass for ocean velocities with Iopt

c Add on a random walk

         theta=Ran(Idum)*2.*Pi
         dah=Ran(Idum)*SQRT(12.*ptah*Dt)

c Take account of active horizontal swimming

         if (species.ne.'passive'.and.swimming.ne.'n') then
          ptu=ptu-(swimflagr*swim_speed(i)*ptu)+
     &            (swimflaga*swimu(ibox)*swim_speed(i))
          ptv=ptv-(swimflagr*swim_speed(i)*ptv)+
     &            (swimflaga*swimv(ibox)*swim_speed(i)) 
         endif       

c Hence distance moved in the x direction and new position

         if (Ah.lt.0.) then
          Xdist=ptu*Dt + Drift1(ibox,1,ilev)*Dt + Drift2(ibox,1,ilev)*Dt
     &          + dah*cos(theta)
         else
          Xdist=ptu*Dt + dah*cos(theta)
         endif
         If (curvi.eq.0.or.curvi.eq.2) then
           Xdgrid=(180.*Xdist/(Pi*Radius_E*plat))/dch
         Else
           If (pom) then
             mdch=(alphat*dchcu(x4t,y4t) + betat*dchcu(x2t,y2t) +
     &           gammat*dchcu(x1t,y1t) + deltat*dchcu(x3t,y3t))
             Xdgrid=(Xdist/(Radius_E*plat))/mdch
           Else
             Xdgrid=(180.*Xdist/(Pi*Radius_E*plat))/dchcu(x1t,y1t)
           Endif
         Endif
         Xelev(i)=Xelev(i)+Xdgrid
         Xu(i)=Xu(i)+Xdgrid
         Xv(i)=Xv(i)+Xdgrid

c Hence distance moved in the y direction and new position

         if (Ah.lt.0.) then
          Ydist=ptv*Dt + Drift1(ibox,2,ilev)*Dt + Drift2(ibox,2,ilev)*Dt
     &          + dah*sin(theta)
         else
          Ydist=ptv*Dt + dah*sin(theta)
         endif
         If (curvi.eq.0) then
           Ydgrid=(180.*Ydist/(Pi*Radius_E))/dfi
         Elseif (curvi.eq.1) then
           If (pom) then
             mdfi=(alphat*dfic(x4t,y4t) + betat*dfic(x2t,y2t) +
     &           gammat*dfic(x1t,y1t) + deltat*dfic(x3t,y3t))
             Ydgrid=(Ydist/Radius_E)/mdfi
           Else
             Ydgrid=(180.*Ydist/(Pi*Radius_E))/dfic(x1t,y1t)
           Endif
         Elseif (curvi.eq.2) then
           Ydgrid=(180.*Ydist/(Pi*Radius_E))/dlat(y1t)
         Endif
         Yelev(i)=Yelev(i)-Ydgrid
         Yu(i)=Yu(i)-Ydgrid
         Yv(i)=Yv(i)-Ydgrid

c Find 4 surrounding W-points
c Bypass this section if doing a depth-mean calculation...

         If (veltyp.eq.'N'.or.veltyp.eq.'n') then

c Also bypass parts of this section if using fixed depth particles (drifters)

          If (btyp.ne.'c'.and.btyp.ne.'C') then

c If occam, just set ptw to 0.0
c Also for POLCOMS as w requires rotating onto true vertical; add this
c in a future update

           If (occam.or.polcoms) then

            ptw=0.0
    
c Also set to 0 if using ice velocities with nemo

           elseif (flagi.eq.0) then

            ptw=0.0

           else

            if (polcoms.or.nemo) then
             If (W(i1t,2).eq.0.or.W(i2t,2).eq.0.or.W(i3t,2).eq.0.or.
     &           W(i4t,2).eq.0) then
              Call Tinterpvel(3,i1t,i2t,i3t,i4t,velcount)
             Endif
            Endif

c Find velocity at particle depth.

            flag=3
            Landflag=0
            if (nemo) then
             Call Calcvel5(i,flag,x1t,y1t,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                     Landflag)
            else
             Call Calcvel1(i,flag,x1t,y1t,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                     Landflag)
            endif

c Hence velocity at particle position

            If (Landflag.eq.0) then
             ptw=(alphat*vsp4 + betat*vsp2 + gammat*vsp1 + deltat*vsp3)
            Else
             If (polcoms) then
              Do k=1,kb-1
               zt(k)=z(ibox,k)*Hfull(ibox)
               s1(k)=w(ibox,k)
              Enddo
              Call Intpol(zt,s1,deppart(i),s2,kb-1)
              ptw=s2
             Else
              ptw=0.;
             Endif
            Endif

           Endif

c Calculate particle buoyancy (eggs/oil droplets) based on density difference 
c from previous Dt. Uses Stoke's law.
c Note viscosity is currently set constant for T~10deg, S~35
c Future updates should incorporate calculation of visc dependent on T and S

           if (stage_spec(part_st(i)).eq.'e') then
c Update only at biological times
            if (btyp.eq.'d'.or.btyp.eq.'D') then
             if (timerem.eq.0.or.T.eq.(Pstart_T(i)+Dt)) then
              if (Tflag.eq.1) then
c Don't try to calculate buoyancy at first timestep as temperatures and 
c salinities at previous timestep aren't known.
               bwlarvae(i)=0.
              else
               Call Calc_dens(i,rho)
               bwlarvae(i)=grav*(eggd(i)*eggd(i))*(rho-bwegg(i))
     &                     /(18.*visc)
              endif
             endif
            elseif (btyp.eq.'v'.or.btyp.eq.'V') then
             bwlarvae(i)=bwegg(i)
            endif
           endif

c Hence distance moved in the z direction
c No particle buoyancy term in ice

           daz=SQRT(2.*ptkh*Dt)
           if (Az.lt.0.) then
            Zdist=ptw*Dt + Drift2(ibox,3,ilev)*Dt + 
     &            ((2.*Ran(Idum)-1.)*SQRT(3.)*daz + 
     &             bwlarvae(i)*Dt)*flagi
           else
            Zdist=ptw*Dt + ((2.*Ran(Idum)-1.)*SQRT(3.)*daz + 
     &                      bwlarvae(i)*Dt)*flagi
           endif

c Else particle stays at constant depth
          else
           Zdist=0.
          endif

c Newz is now the particle depth
          Newz=Deppart(i)-Zdist

         Endif

c Test new particle position and adjust if necessary

c If circumpolar, allow particles to loop round grid

         If (circumpolar) then
          If (Xelev(i).le.0.5) Xelev(i)=Xelev(i)+im
          If (Xu(i).le.0.5) Xu(i)=Xu(i)+im
          If (Xv(i).le.0.5) Xv(i)=Xv(i)+im
          If (Xelev(i).gt.(im+0.5)) Xelev(i)=Xelev(i)-im
          If (Xu(i).gt.(im+0.5)) Xu(i)=Xu(i)-im
          If (Xv(i).gt.(im+0.5)) Xv(i)=Xv(i)-im
         Endif  

c If particle is now on a permanent land point, move back to position at
c previous timestep. Treat land as a reflecting boundary.
c If it is predicted to move onto a tidal flat point, move back to position
c at previous timestep (included in mind comparison).
c Also move back if predicted to go onto a point where depth is less
c than the minimum particle depth

         ibox=(IDNINT(Yelev(i))-1)*im + IDNINT(Xelev(i))
         iboxu=(IDNINT(OldYelev(i))-1)*im + IDNINT(Xelev(i))
         iboxv=(IDNINT(Yelev(i))-1)*im + IDNINT(OldXelev(i))
         if ((mz(ibox).eq.0).or.(Hfull(ibox).lt.mind)) Back=1
c POLCOMS requires check on umask too
         if ((polcoms.and.mu(ibox).eq.0)) Back=1
c If both new depth and old depth are deeper than new location, try
c moving with u and v components separately
c If still moving onto land, move particles back
         if (Hfull(ibox).lt.deppart(i).and.Hfull(ibox).lt.Newz) then
          Back=1
          if (bflag(i).eq.0) then 
           bathy_flag(i)=bathy_flag(i)+1
           bflag(i)=1
          endif
          if (Hfull(iboxv).gt.deppart(i).or.Hfull(iboxv).gt.Newz) then
           Back=2
          endif
          if (Hfull(iboxu).gt.deppart(i).or.Hfull(iboxu).gt.Newz) then
           Back=3
          endif
         endif

         if (Back.gt.0) then
c Put back to last coastal sea point
          if ((mz(ibox).eq.0).and.beaching) then
c If beaching set, beach particles if on land instead of reflecting
            Part_sw(i)=-2
          endif
          if (Back.eq.1.or.Back.eq.2) then
           Xelev(i)=OldXelev(i)
           Xu(i)=OldXu(i)
           Xv(i)=OldXv(i)
          endif
          if (Back.eq.1.or.Back.eq.3) then
           Yelev(i)=OldYelev(i)
           Yu(i)=OldYu(i)
           Yv(i)=OldYv(i)
          endif
         endif

         if (Part_sw(i).gt.-2) then

c Calculate new position

          x1t=INT(Xelev(i))
          y1t=INT(Yelev(i))
          Call Findpts(x1t,x2t,x3t,x4t,y1t,y2t,y3t,y4t,i1t,i2t,i3t,i4t)

c Find weighted distances

          pt=Xelev(i)-x1t
          qt=Yelev(i)-y1t
          alphat=pt*qt
          betat=pt*(1-qt)
          gammat=(1-pt)*(1-qt)
          deltat=(1-pt)*qt

          ibox=(IDNINT(Yelev(i))-1)*im + IDNINT(Xelev(i))

         endif

c Don't allow particles to go above the surface or below the bed
c Note Newz is particle depth at this point
c Bypass this section for depth-mean calculations
c Only adjust vertical position if btyp ne 'c'

         If (veltyp.eq.'N'.or.veltyp.eq.'n') then

          If (btyp.ne.'c'.and.btyp.ne.'C') then

c Option for oil spill - set particles above surface to 0.5m depth
           if (oil.and.Newz.lt.0.) Newz=0.5
c Reflect particles from sea surface
           if (Newz.lt.0.) Newz=-Newz
c Reflect particles from sea bed
           if (Newz.gt.Hfull(ibox)) then
            if (Back.gt.1) then
             Newz=Deppart(i)
            else
             Newz=2.*Hfull(ibox) - Newz
            endif
            if (bflag(i).eq.0) then
             bathy_flag(i)=bathy_flag(i)+1
             bflag(i)=1
            endif
           endif
c Correct near-surface particles that are still being a pest
           if (Newz.lt.0.) Newz=5
c Also keep particle within range defined in biology_spec.dat (put 
c back to previous depth if not)
           if (Iopt.and.nemo.and.Icep(i).ge.0.15) then
            if (Newz.lt.mind.or.Newz.gt.maxdi)
     &       Newz=Deppart(i)
           elseif (stage_spec(part_st(i)).eq.'c') then
c Option for diapausing copepods: first reset diapause flags if 
c particles have reached their target depth
            if (Newz.ge.partd_dia(i).and.dia_flag(i).eq.3)
     &          dia_flag(i)=1
            if (Newz.le.Deppart_st(i).and.dia_flag(i).eq.4)
     &          dia_flag(i)=0
c Then limit particles to either upper or lower depth range (put
c back to previous depth if not within range)
            if ((dia_flag(i).eq.0.and.(Newz.lt.mind.or.
     &           Newz.gt.maxd)).or.(dia_flag(i).eq.1.and.
     &           Newz.lt.mind_dia.or.Newz.gt.maxd_dia).or.
     &          (dia_flag(i).eq.3.and.Newz.lt.mind).or.
     &          (dia_flag(i).eq.4.and.Newz.gt.maxd_dia))
     &       Newz=Deppart(i)
           else
            if (Newz.lt.mind.or.Newz.gt.maxd) then
             Newz=Deppart(i)
            endif
           endif
c Now set Deppart to Newz
           Deppart(i)=Newz
c Do one final check for particles above the minimum depth...
           if (stage_spec(part_st(i)).eq.'c') then
            if (dia_flag(i).eq.0.and.Deppart(i).lt.mind)
     &          Deppart(i)=mind
            if (dia_flag(i).eq.1.and.Deppart(i).lt.mind_dia)
     &          Deppart(i)=mind_dia
           else
            if (Deppart(i).lt.mind) Deppart(i)=mind
           endif
c Final check at sea bed to eliminate spurious vertical jumps
           if (Deppart(i).ge.Hfull(ibox)) then
            Deppart(i)=Hfull(ibox)-5
            if (bflag(i).eq.0) then
             bathy_flag(i)=bathy_flag(i)+1
             bflag(i)=1
            endif
           endif
c Check height above seabed for later benthic interaction analyses
           if ((Hfull(ibox)-Deppart(i)).lt.20.and.bflag(i).eq.0) then 
            bathy_flag(i)=bathy_flag(i)+1
            bflag(i)=1
           endif

          Endif

c Check and write warning if particle is below bed/above surface/not 
c within mind and maxd - no check for maxdi to allow particles to
c move from ice-free to ice-covered areas

          If (stage_spec(part_st(i)).eq.'c') then
           If ((Deppart(i).gt.Hfull(ibox)).or.
     &         ((dia_flag(i).eq.0).and.
     &          (Deppart(i).lt.mind.or.Deppart(i).gt.maxd)).or.
     &         ((dia_flag(i).eq.1).and.
     &          (Deppart(i).lt.mind_dia.or.Deppart(i).gt.maxd_dia)))
     &       then
            write (55,*) 'Warning! Copepod not meeting depth criteria'
            write (55,*) 'T, particle, Deppart, Hfull, i, j : ',T,i,
     &       Deppart(i),Hfull(ibox),IDNINT(Xelev(i)),IDNINT(Yelev(i))
           Endif
          Elseif ((Deppart(i).gt.Hfull(ibox)).or.(Deppart(i).lt.mind)
     &             .or.(Deppart(i).gt.maxd)) then
           write (55,*) 'Warning! Particle not meeting depth criteria'
           write (55,*) 'T, particle, Deppart, Hfull, i, j : ',T,i,
     &       Deppart(i),Hfull(ibox),IDNINT(Xelev(i)),IDNINT(Yelev(i))
          Endif

         Endif

c Find temperature and salinity at new particle position

         if (TSopt) then
          if (polcoms.or.nemo) then 
           If (Tem(i1t,1).eq.0.or.Tem(i2t,1).eq.0.or.Tem(i3t,1).eq.0.or.
     &        Tem(i4t,1).eq.0) then
            Call Tinterpvel(4,i1t,i2t,i3t,i4t,velcount)
           Endif
          Endif

          flag=4
          Landflag=0
          if (nemo) then
           Call Calcvel5(i,flag,x1t,y1t,tsp1,tsp2,tsp3,tsp4,Hfull,
     &                   Landflag)
          else
           If (veltyp.eq.'N'.or.veltyp.eq.'n') then
            Call Calcvel1(i,flag,x1t,y1t,tsp1,tsp2,tsp3,tsp4,Hfull,
     &                    Landflag)
           Else
            Call Calcvel2(flag,x1t,y1t,tsp1,tsp2,tsp3,tsp4,Hfull,
     &                    Landflag)
           Endif
          endif
          If (Landflag.eq.0) then
           Temper(i)=(alphat*tsp4 + betat*tsp2 + gammat*tsp1 + 
     &                deltat*tsp3)
          Else
           if (polcoms) then
            Do k=1,kb-1
             zt(k)=zz(ibox,k)*Hfull(ibox)
             s1(k)=tem(ibox,k)
            Enddo
            If (veltyp.eq.'N'.or.veltyp.eq.'n') then
             Call Intpol(zt,s1,deppart(i),s2,kb-1)
            Else
             zt(kb)=-1.*Hfull(ibox)
             s2=0.
             Do k=1,kb-1
              dzt(k)=zt(k)-zt(k+1)
              vz(k)=s1(k)*dzt(k)
              s2=s2+vz(k)
             Enddo
             s2=s2/Hfull(ibox)
            Endif
            Temper(i)=s2
           elseif (nemo) then
            Do k=1,kb-1
             if (nemolev(ibox,k).lt.0.) then
              nlev=k-1
              goto 360
             endif
             zt(k)=nemolev(ibox,k)*(-1)
             s1(k)=tem(ibox,k)
            Enddo
 360        Call Intpol(zt,s1,deppart(i),s2,nlev)
            Temper(i)=s2
           endif
          Endif
          Tempav_bio(i)=Tempav_bio(i)+Temper(i)/
     &                   (biocalcfreq(Part_st(i))/Dt)
          Tempav(i)=Tempav(i)+Temper(i)/(Outdt/Dt)

c Find salinity at new particle position

          if (polcoms.or.nemo) then 
           If (Sal(i1t,1).eq.0.or.Sal(i2t,1).eq.0.or.Sal(i3t,1).eq.0.or.
     &        Sal(i4t,1).eq.0) then
            Call Tinterpvel(7,i1t,i2t,i3t,i4t,velcount)
           Endif
          Endif
 
          Landflag=0
          if (nemo) then
           flag=5
           Call Calcvel5(i,flag,x1t,y1t,ssp1,ssp2,ssp3,ssp4,Hfull,
     &                   Landflag)
          else
           flag=7
           If (veltyp.eq.'N'.or.veltyp.eq.'n') then
            Call Calcvel1(i,flag,x1t,y1t,ssp1,ssp2,ssp3,ssp4,Hfull,
     &                    Landflag)
           Else
            Call Calcvel2(flag,x1t,y1t,ssp1,ssp2,ssp3,ssp4,Hfull,
     &                    Landflag)
           Endif
          endif
          If (Landflag.eq.0) then
           Salin(i)=(alphat*ssp4 + betat*ssp2 + gammat*ssp1 + 
     &               deltat*ssp3)
          Else
           if (polcoms) then
            Do k=1,kb-1
             zt(k)=zz(ibox,k)*Hfull(ibox)
             s1(k)=sal(ibox,k)
            Enddo
            If (veltyp.eq.'N'.or.veltyp.eq.'n') then
             Call Intpol(zt,s1,deppart(i),s2,kb-1)
            Else
             zt(kb)=-1.*Hfull(ibox)
             s2=0.
             Do k=1,kb-1
              dzt(k)=zt(k)-zt(k+1)
              vz(k)=s1(k)*dzt(k)
              s2=s2+vz(k)
             Enddo
             s2=s2/Hfull(ibox)
            Endif
            Salin(i)=s2
           elseif (nemo) then
            Do k=1,kb-1
             if (nemolev(ibox,k).lt.0.) then
              nlev=k-1
              goto 370
             endif
             zt(k)=nemolev(ibox,k)*(-1)
             s1(k)=sal(ibox,k)
            Enddo
 370        Call Intpol(zt,s1,deppart(i),s2,nlev)
            Salin(i)=s2
           endif
          Endif
          Salinav(i)=Salinav(i)+Salin(i)/(Outdt/Dt)

         Endif

 150    continue

        Return
        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Calc_drift(Hfull)

c This subroutine calculates the particle drift (diffusion) based on
c hydrodynamic model output fields Ah and Kh

        use param

        implicit none

        Integer :: i,j,k,im1,ip1,imn,ipn,izm1,izp1
        Real :: delx, dely
        Real :: Km1, Kp1, Kmn, Kpn
        Real :: Hfull(Mx2)

        do 5 i=1,Mx2
         do 5 j=1,3
          do 5 k=1,kb-1
           Drift1(i,j,k)=0.
           Drift2(i,j,k)=0.
   5    continue

        dely=2.*deg2rad*dfi*Radius_E

        do 10 k=2,kb-1
         do 10 j=1,nz
          if (polcoms) then
           i=iz(j)
           im1=i-1
           ip1=i+1
           imn=i-im
           ipn=i+im
           izm1=k-1
           izp1=k+1
           delx=2.*deg2rad*dch*Radius_E*cspz(i)
           Drift2(i,1,k)=(PAh(ip1,k)-PAh(im1,k))/delx
           Drift2(i,2,k)=(PAh(ipn,k)-PAh(imn,k))/dely
           if (k.eq.kb-1) then
            Drift2(i,3,k)=(PKh(i,izm1)-PKh(i,k))/
     &                     abs(z(i,k-1)*Hfull(i)-z(i,k)*Hfull(i))
           elseif (k.eq.2) then
            Drift2(i,3,k)=(PKh(i,k)-PKh(i,izp1))/
     &                     abs(z(i,k)*Hfull(i)-z(i,k+1)*Hfull(i))
           else
            Drift2(i,3,k)=(PKh(i,izm1)-PKh(i,izp1))/
     &                     abs(z(i,k-1)*Hfull(i)-z(i,k+1)*Hfull(i))
           endif
           Drift1(i,1,k)=(PAh(i,k)/Hfull(i))*(Hfull(ip1)-Hfull(im1))
     &                          /delx
           Drift1(i,2,k)=(PAh(i,k)/Hfull(i))*(Hfull(ipn)-Hfull(imn))
     &                          /dely
          else
           do 20 i=iz1(j),iz2(j)
            im1=i-1
            ip1=i+1
            imn=i-im
            ipn=i+im
            izm1=k-1
            izp1=k+1
            delx=2.*deg2rad*dch*Radius_E*cspz(i)
            Drift2(i,1,k)=(PAh(ip1,k)-PAh(im1,k))/delx
            Drift2(i,2,k)=(PAh(ipn,k)-PAh(imn,k))/dely
            if (k.eq.kb-1) then
             Drift2(i,3,k)=(PKh(i,izm1)-PKh(i,k))/
     &                      abs(zz(i,k-1)*Hfull(i)-zz(i,k)*Hfull(i))
            elseif (k.eq.2) then
             Drift2(i,3,k)=(PKh(i,k)-PKh(i,izp1))/
     &                      abs(zz(i,k)*Hfull(i)-zz(i,k+1)*Hfull(i))
            else
             Drift2(i,3,k)=(PKh(i,izm1)-PKh(i,izp1))/
     &                      abs(zz(i,k-1)*Hfull(i)-zz(i,k+1)*Hfull(i))
            endif
            Drift1(i,1,k)=(PAh(i,k)/Hfull(i))*(Hfull(ip1)-Hfull(im1))
     &                           /delx
            Drift1(i,2,k)=(PAh(i,k)/Hfull(i))*(Hfull(ipn)-Hfull(imn))
     &                           /dely
   20      continue
          endif
   10   continue

        Return
        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Calcvel1 (ipt,flag,x1,y1,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                       Landflag)

c This subroutine interpolates velocities/temperatures/salinity/
c diffusivities at sigma levels for 4 surrounding points to a specified 
c depth
c Flag : 1=u; 2=v; 3=w; 4=Tem; 5=PKh; 6=PAh; 7=Sal; 8=POLCOMS u;
c 9=POLCOMS v

        use param

        implicit none

        Real :: zt(kb),s1(kb),s2
        Real :: Hfull(Mx2)
        Real :: vsp1,vsp2,vsp3,vsp4

        Integer :: flag,count,x1,y1,Landflag,i,j,k,ii,iim1,iipn,iipnm1,
     &             ipt

C Interpolate velocity/temperature/diffusivity at each surrounding point
C vertically to depth of particle

        count=1
        Do 10 j=y1,y1+1
        Do 10 i=x1,x1+1
         ii=(j-1)*im + i
         iim1=(j-1)*im + i-1
         iipn=j*im + i
         iipnm1=j*im + i-1
         Do 20 k=1,kb-1
          if (flag.eq.1) then
           if (Hfull(ii).lt.1.0.or.Hfull(iim1).lt.1.0.or.
     &          deppart(ipt).gt.((Hfull(ii)+Hfull(ii-1))*0.5)) then
                s2=0.0
            goto 15
           endif
           zt(k)=zz(ii,k)*0.5*(Hfull(ii)+Hfull(iim1))
           s1(k)=u(ii,k)
          elseif (flag.eq.2) then
           if (Hfull(ii).lt.1.0.or.Hfull(iipn).lt.1.0.or.
     &          deppart(ipt).gt.((Hfull(ii)+Hfull(iipn))*0.5)) then
                s2=0.0
            goto 15
           endif
           zt(k)=zz(ii,k)*0.5*(Hfull(ii)+Hfull(iipn))
           s1(k)=v(ii,k)
          elseif (flag.eq.3) then
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           zt(k)=z(ii,k)*Hfull(ii)
           s1(k)=w(ii,k)
          elseif (flag.eq.4) then
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=tem(ii,kb-1)
                if (Hfull(ii).lt.1.0) then
                  Landflag=1
                  goto 999
                endif
            goto 15
           endif
           zt(k)=zz(ii,k)*Hfull(ii)
           s1(k)=tem(ii,k)
          elseif (flag.eq.5) then
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=PKh(ii,kb-1)
                if (Hfull(ii).lt.1.0) then
                  Landflag=1
                  goto 999
                endif
            goto 15
           endif
           zt(k)=zz(ii,k)*Hfull(ii)
           s1(k)=PKh(ii,k)
           if (k.eq.1) s1(k)=PKh(ii,2)
          elseif (flag.eq.6) then
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=PAh(ii,kb-1)
                if (Hfull(ii).lt.1.0) then
                  Landflag=1
                  goto 999
                endif
            goto 15
           endif
           zt(k)=zz(ii,k)*Hfull(ii)
           s1(k)=PAh(ii,k)
          elseif (flag.eq.7) then
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=sal(ii,kb-1)
                if (Hfull(ii).lt.1.0) then
                  Landflag=1
                  goto 999
                endif
            goto 15
           endif
           zt(k)=zz(ii,k)*Hfull(ii)
           s1(k)=sal(ii,k)
          elseif (flag.eq.8) then
c Option for POLCOMS u
           if (Hfull(ii).lt.1.0.or.Hfull(iipn).lt.1.0.or.
     &         Hfull(iim1).lt.1.0.or.Hfull(iipnm1).lt.1.0.or.
     &         deppart(ipt).gt.((Hfull(ii)+Hfull(iipn)+Hfull(iim1)
     &         +Hfull(iipnm1))*0.25)) then
                s2=0.0
            goto 15
           endif
           zt(k)=zzv(ii,k)*0.25*(Hfull(ii)+Hfull(iipn)+Hfull(iim1)+
     &                           Hfull(iipnm1))
           s1(k)=u(ii,k)
          elseif (flag.eq.9) then
c Option for POLCOMS v
           if (Hfull(ii).lt.1.0.or.Hfull(iipn).lt.1.0.or.
     &         Hfull(iim1).lt.1.0.or.Hfull(iipnm1).lt.1.0.or.
     &         deppart(ipt).gt.((Hfull(ii)+Hfull(iipn)+Hfull(iim1)
     &         +Hfull(iipnm1))*0.25)) then
                s2=0.0
            goto 15
           endif
           zt(k)=zzv(ii,k)*0.25*(Hfull(ii)+Hfull(iipn)+Hfull(iim1)+
     &                           Hfull(iipnm1))
           s1(k)=v(ii,k)
          endif
   20    Continue
         Call Intpol(zt,s1,deppart(ipt),s2,kb-1)
   15    if (count.eq.1) vsp1=s2
         if (count.eq.2) vsp2=s2
         if (count.eq.3) vsp3=s2
         if (count.eq.4) vsp4=s2
         count=count+1
   10   Continue

 999    Return

        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Calcvel2 (flag,x1,y1,vsp1,vsp2,vsp3,vsp4,Hfull,
     &                       Landflag)

c This subroutine calculates depth-mean velocities/temperatures/
c diffusivities for 4 surrounding points
c Flag : 1=u; 2=v; 3=w; 4=Tem; 5=PKh; 6=PAh; 7=Sal; 8=POLCOMS u;
c 9=POLCOMS v

        use param

        implicit none

        Real :: zt(kb),s1(kb),vz(kb),vsum,Hav,dzt(kb)
        Real :: Hfull(Mx2)
        Real :: vsp1,vsp2,vsp3,vsp4

        Integer :: flag,count,x1,y1,Landflag,i,j,k,ii,iim1,iipn,iipnm1,
     &             ipt

C Vertically average variable at each surrounding point

        count=1
        Do 10 j=y1,y1+1
        Do 10 i=x1,x1+1
         ii=(j-1)*im + i
         iim1=(j-1)*im + i-1
         iipn=j*im + i
         iipnm1=j*im + i-1
         vsum=0.
         Hav=0.
         if (flag.eq.1) then
          Hav=0.5*(Hfull(ii)+Hfull(iim1))
         elseif (flag.eq.2) then
          Hav=0.5*(Hfull(ii)+Hfull(iipn))
         elseif (flag.ge.3.and.flag.le.7) then
          Hav=Hfull(ii)
         elseif (flag.eq.8.or.flag.eq.9) then
          Hav=0.25*(Hfull(ii)+Hfull(iipn)+Hfull(iim1)+Hfull(iipnm1))
         endif
         zt(kb)=-1.*Hav
         Do 20 k=1,kb-1
          if (flag.eq.1) then
           if (Hfull(ii).lt.1.0.or.Hfull(iim1).lt.1.0) then
                vsum=0.0
            goto 15
           endif
           zt(k)=zz(ii,k)*Hav
           s1(k)=u(ii,k)
          elseif (flag.eq.2) then
           if (Hfull(ii).lt.1.0.or.Hfull(iipn).lt.1.0) then
                vsum=0.0
            goto 15
           endif
           zt(k)=zz(ii,k)*Hav
           s1(k)=v(ii,k)
          elseif (flag.eq.3) then
           if (Hfull(ii).lt.1.0) then
                vsum=0.0
                Landflag=1
            goto 999
           endif
           zt(k)=z(ii,k)*Hav
           s1(k)=w(ii,k)
          elseif (flag.eq.4) then
           if (Hfull(ii).lt.1.0) then
                vsum=0.0
                Landflag=1
            goto 999
           endif
           zt(k)=zz(ii,k)*Hav
           s1(k)=tem(ii,k)
          elseif (flag.eq.5) then
           if (Hfull(ii).lt.1.0) then
                vsum=0.0
                Landflag=1
            goto 999
           endif
           zt(k)=z(ii,k)*Hav
           s1(k)=PKh(ii,k)
          elseif (flag.eq.6) then
          if (Hfull(ii).lt.1.0) then
                vsum=0.0
                Landflag=1
            goto 999
           endif
           zt(k)=z(ii,k)*Hav
           s1(k)=PAh(ii,k)
          elseif (flag.eq.7) then
           if (Hfull(ii).lt.1.0) then
                vsum=0.0
                Landflag=1
            goto 999
           endif
           zt(k)=zz(ii,k)*Hav
           s1(k)=sal(ii,k)
          elseif (flag.eq.8) then
           if (Hfull(ii).lt.1.0.or.Hfull(iim1).lt.1.0.or.
     &         Hfull(iipn).lt.1.0.or.Hfull(iipnm1).lt.1.0) then
                vsum=0.0
            goto 15
           endif
           zt(k)=zzv(ii,k)*Hav
           s1(k)=u(ii,k)
          elseif (flag.eq.9) then
           if (Hfull(ii).lt.1.0.or.Hfull(iim1).lt.1.0.or.
     &         Hfull(iipn).lt.1.0.or.Hfull(iipnm1).lt.1.0) then
                vsum=0.0
            goto 15
           endif
           zt(k)=zzv(ii,k)*Hav
           s1(k)=v(ii,k)
          endif
   20    Continue
         Do 30 k=1,kb-1
          dzt(k)=zt(k)-zt(k+1)
          vz(k)=s1(k)*dzt(k)
          vsum=vsum+vz(k)
   30    Continue
         vsum=vsum/Hav
   15    if (count.eq.1) vsp1=vsum
         if (count.eq.2) vsp2=vsum
         if (count.eq.3) vsp3=vsum
         if (count.eq.4) vsp4=vsum
         count=count+1
   10   Continue

 999    Return

        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Calcvel3 (ipt,flag,x1,y1,vsp1,vsp2,vsp3,vsp4,
     &                       Hfull,Landflag,ilev)

c This subroutine extracts velocities at particle-specific OCCAM levels
c at 4 surrounding points
c Flag : 1=u; 2=v;

        use param

        implicit none

        Real :: s2
        Real :: Hfull(Mx2)
        Real :: vsp1,vsp2,vsp3,vsp4

        Integer :: flag,count,x1,y1,Landflag,ilev,ii,iimn,iip1,iimnp1
        Integer :: i,j,ipt

C Extract velocity at each surrounding point at depth of particle

        count=1
        Do 10 j=y1,y1+1
        Do 10 i=x1,x1+1
         ii=(j-1)*im + i
         iimn=(j-2)*im + i
         iip1=(j-1)*im + i+1
         iimnp1=(j-2)*im + i+1
         s2=0.
         If (circumpolar.and.i.gt.im) then
          ii=(j-1)*im + i - im
          iimn=(j-2)*im + i - im
         Endif
         If (circumpolar.and.(i+1).gt.im) then
          iip1=(j-1)*im + i+1-im
          iimnp1=(j-2)*im + i+1-im
         Endif
         if (flag.eq.1) then
c Option for OCCAM u
           if (Hfull(ii).lt.1.0.or.Hfull(iimn).lt.1.0.or.
     &         Hfull(iip1).lt.1.0.or.Hfull(iimnp1).lt.1.0.or.
     &          deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           s2=u(ii,ilev)
         elseif (flag.eq.2) then
c Option for OCCAM v
           if (Hfull(ii).lt.1.0.or.Hfull(iimn).lt.1.0.or.
     &         Hfull(iip1).lt.1.0.or.Hfull(iimnp1).lt.1.0.or.
     &          deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           s2=v(ii,ilev)
         endif
   15    if (count.eq.1) vsp1=s2
         if (count.eq.2) vsp2=s2
         if (count.eq.3) vsp3=s2
         if (count.eq.4) vsp4=s2
         count=count+1
   10   Continue

 999    Return

        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Calcvel4 (ipt,flag,x1,y1,vsp1,vsp2,vsp3,vsp4,
     &                       Hfull,Landflag,ilev)

c This subroutine extracts velocities at particle-specific NEMO levels
c at 4 surrounding points
c Flag : 1=u; 2=v; 3=w; 4=ice; 5=uice; 6=vice

        use param

        implicit none

        Real :: s2
        Real :: Hfull(Mx2)
        Real :: vsp1,vsp2,vsp3,vsp4

        Integer :: flag,count,x1,y1,Landflag,ilev,ii,iimn,iip1,i,j,ipt

C Extract velocity at each surrounding point at depth of particle

        count=1
        Do 10 j=y1,y1+1
        Do 10 i=x1,x1+1
         ii=(j-1)*im + i
         iimn=(j-2)*im + i
         iip1=(j-1)*im + i+1
         s2=0.
         if (flag.eq.1) then
c Option for NEMO u
           if (Hfull(ii).lt.1.0.or.Hfull(iip1).lt.1.0.or.
     &          deppart(ipt).gt.((Hfull(ii)+Hfull(iip1))*0.5)) then
                s2=0.0
            goto 15
           endif
           s2=u(ii,ilev)
         elseif (flag.eq.2) then
c Option for NEMO v
           if (Hfull(ii).lt.1.0.or.Hfull(iimn).lt.1.0.or.
     &          deppart(ipt).gt.((Hfull(ii)+Hfull(iimn))*0.5)) then
                s2=0.0
            goto 15
           endif
           s2=v(ii,ilev)
         elseif (flag.eq.3) then
c Option for NEMO w
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           s2=w(ii,ilev)
         elseif (flag.eq.4) then
c Option for NEMO ice concentration
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           s2=Ice(ii)
         elseif (flag.eq.5) then
c Option for NEMO ice u
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           s2=Uice(ii)
         elseif (flag.eq.6) then
c Option for NEMO ice v
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           s2=Vice(ii)
         endif
   15    if (count.eq.1) vsp1=s2
         if (count.eq.2) vsp2=s2
         if (count.eq.3) vsp3=s2
         if (count.eq.4) vsp4=s2
         count=count+1
   10   Continue

 999    Return

        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Calcvel5 (ipt,flag,x1,y1,vsp1,vsp2,vsp3,vsp4,
     &                       Hfull,Landflag)

c This subroutine interpolates velocities at 4 surrounding points to a
c specified depth
c Also does T and S
c NEMO-specific
c Flag : 1=u; 2=v; 3=w; 4=T; 5=S

        use param

        implicit none

        Real :: zt(kb),s1(kb),s2
        Real :: Hfull(Mx2)
        Real :: vsp1,vsp2,vsp3,vsp4

        Integer :: flag,count,x1,y1,ii,iimn,iip1,i,j,k,ipt,nlev,Landflag

C Extract velocity at each surrounding point at depth of particle

        count=1
        Do 10 j=y1,y1+1
        Do 10 i=x1,x1+1
         ii=(j-1)*im + i
         iimn=(j-2)*im + i
         iip1=(j-1)*im + i+1
         s2=0.
         nlev=kb-1
         Do 20 k=1,kb-1
          if (flag.eq.1) then
c Option for NEMO u
           if (Hfull(ii).lt.1.0.or.Hfull(iip1).lt.1.0.or.
     &          deppart(ipt).gt.((Hfull(ii)+Hfull(iip1))*0.5)) then
                s2=0.0
            goto 15
           endif
           if (nemolev(ii,k).lt.0.) then
             nlev=k-1
             goto 30
           endif
           zt(k)=nemolev(ii,k)*(-1)
           s1(k)=u(ii,k)
          elseif (flag.eq.2) then
c Option for NEMO v
           if (Hfull(ii).lt.1.0.or.Hfull(iimn).lt.1.0.or.
     &          deppart(ipt).gt.((Hfull(ii)+Hfull(iimn))*0.5)) then
                s2=0.0
            goto 15
           endif
           if (nemolev(ii,k).lt.0.) then
             nlev=k-1
             goto 30
           endif
           zt(k)=nemolev(ii,k)*(-1)
           s1(k)=v(ii,k)
          elseif (flag.eq.3) then
c Option for NEMO w
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
                s2=0.0
            goto 15
           endif
           if (nemowlev(ii,k).lt.0.) then
             nlev=k-1
             goto 30
           endif
           zt(k)=nemowlev(ii,k)*(-1)
           s1(k)=w(ii,k)
          elseif (flag.eq.4) then
c Option for NEMO temperature
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
             Landflag=1
             goto 999
           endif
           if (nemolev(ii,k).lt.0.) then
             nlev=k-1
             goto 30
           endif
           zt(k)=nemolev(ii,k)*(-1)
           s1(k)=tem(ii,k)
          elseif (flag.eq.5) then
c Option for NEMO salinity
           if (Hfull(ii).lt.1.0.or.deppart(ipt).gt.Hfull(ii)) then
             Landflag=1
             goto 999
           endif
           if (nemolev(ii,k).lt.0.) then
             nlev=k-1
             goto 30
           endif
           zt(k)=nemolev(ii,k)*(-1)
           s1(k)=sal(ii,k)
          endif
   20    Continue
   30    Call Intpol(zt,s1,deppart(ipt),s2,nlev)
   15    if (count.eq.1) vsp1=s2
         if (count.eq.2) vsp2=s2
         if (count.eq.3) vsp3=s2
         if (count.eq.4) vsp4=s2
         count=count+1
   10   Continue

 999    Return

        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Findpts(x1,x2,x3,x4,y1,y2,y3,y4,i1,i2,i3,i4)

C Routine to identify 4 surrounding points

        use param

        implicit none

        Integer :: x1, x2, x3, x4, y1, y2, y3, y4, i1, i2, i3, i4

        x2=x1+1
        y2=y1
        x3=x1
        y3=y1+1
        x4=x1+1
        y4=y1+1
        i1=(y1-1)*im + x1
        i2=(y2-1)*im + x2
        i3=(y3-1)*im + x3
        i4=(y4-1)*im + x4
        If (circumpolar.and.x1.eq.0) then
         i1=(y1-1)*im + im
         i3=(y3-1)*im + im
        Endif
        If (circumpolar.and.x2.gt.im) then
         i2=(y2-1)*im + (x2-im)
         i4=(y4-1)*im + (x4-im)
        Endif

        Return
        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Intpol(X,A,Y,B,N)
C
C       THIS SUBROUTINE CALCULATES B AT CERTAIN DEPTHS FROM A
C        X(N) MUST BE DESCENDING
C        A(X) GIVEN FUNCTION
C        Y THE DESIRED DEPTH
C        B(Y) FOUND BY LINEAR INTERPOLATION AND EXTRAPOLATION
C                          NO EXTRAPOLATION NOW...
C        N NUMBER OF POINTS GIVEN IN X
C
C *******************************************************************

        implicit none

C
C Global variables.
C
        Real :: X(N),A(N),Y,B
        Integer :: N
C
C Local variables.
C
        Integer :: J
        Real :: NegY

        Negy=-1.*Y

C
C EXTRAPOLATION CASES
C
        If (NegY.GT.X(1)) then
          B=A(1)
          Return
        Else if (NegY.LT.X(N)) then
          B=A(N)
          Return
        Endif
C
C INTERPOLATION CASES
C
        Do 60 J=1,N
          If (NegY.GT.X(J).OR.NegY.LT.X(J+1)) goto 60
          B=A(J)-(A(J)-A(J+1))*(X(J)-NegY)/(X(J)-X(J+1))
          Return
 60     Continue

        Return

        End

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        Subroutine Tinterpvel(flag,i1,i2,i3,i4,velcount)

        use param

        implicit none

        Integer :: flag,velcount,i1,i2,i3,i4,ipt,ii(4),tcount,k
        Real :: Dtr, Dtr2, a, b

        Dtr2=float(Dt2)
        Dtr=float(Dt)
        tcount=velcount-Dt
        ii(1)=i1
        ii(2)=i2
        ii(3)=i3
        ii(4)=i4

        Do 10 ipt=1,4

c Calculate vals on half timestep

          if (mu(ii(ipt)).eq.0) goto 15
          if (flag.eq.1) then
            a=Dtr2-(Dtr/2.)-tcount
            b=tcount+(Dtr/2.)
            U(ii(ipt),1:kb)=(a*U1(ii(ipt),1:kb)+b*U2(ii(ipt),1:kb))/Dtr2
          endif

   15	  Continue

          if ((polcoms.or.occam).and.mu(ii(ipt)).eq.0) goto 25
          if (nemo.and.mv(ii(ipt)).eq.0) goto 25
          if (flag.eq.2) then
            a=Dtr2-(Dtr/2.)-tcount
            b=tcount+(Dtr/2.)
            V(ii(ipt),1:kb)=(a*V1(ii(ipt),1:kb)+b*V2(ii(ipt),1:kb))/Dtr2
          endif

   25	  Continue

          if (flag.lt.3.or.mz(ii(ipt)).eq.0) goto 35

          if (flag.eq.3) then
            a=Dtr2-(Dtr/2.)-tcount
            b=tcount+(Dtr/2.)
            W(ii(ipt),1:kb)=(a*W1(ii(ipt),1:kb)+b*W2(ii(ipt),1:kb))/Dtr2
          elseif (flag.eq.4) then
            a=Dtr2-tcount
            b=float(tcount)
            Tem(ii(ipt),1:kb)=(a*T1(ii(ipt),1:kb)+b*T2(ii(ipt),1:kb))/
     &                         Dtr2
          elseif (flag.eq.5) then
            a=Dtr2-tcount
            b=float(tcount)
            PKh(ii(ipt),1:kb)=(a*PKh1(ii(ipt),1:kb)+
     &                         b*PKh2(ii(ipt),1:kb))/Dtr2
          elseif (flag.eq.6) then
            a=Dtr2-tcount
            b=float(tcount)
            PAh(ii(ipt),1:kb)=(a*PAh1(ii(ipt),1:kb)+
     &                         b*PAh2(ii(ipt),1:kb))/Dtr2
          elseif (flag.eq.7) then
            a=Dtr2-tcount
            b=float(tcount)
            Sal(ii(ipt),1:kb)=(a*Sa1(ii(ipt),1:kb)+
     &                         b*Sa2(ii(ipt),1:kb))/Dtr2
          elseif (flag.eq.8) then
            a=Dtr2-tcount
            b=float(tcount)
            Ice(ii(ipt))=(a*Ice1(ii(ipt))+b*Ice2(ii(ipt)))/Dtr2
          elseif (flag.eq.9) then
            a=Dtr2-(Dtr/2.)-tcount
            b=tcount+(Dtr/2.)
            Uice(ii(ipt))=(a*Uice1(ii(ipt))+b*Uice2(ii(ipt)))/Dtr2
          elseif (flag.eq.10) then
            a=Dtr2-(Dtr/2.)-tcount
            b=tcount+(Dtr/2.)
            Vice(ii(ipt))=(a*Vice1(ii(ipt))+b*Vice2(ii(ipt)))/Dtr2
          endif

   35	  Continue

   10	Continue

        return
        end

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx





