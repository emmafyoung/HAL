        Subroutine Input_bio2

c Routine to read biological parameters from file

        use param

        implicit none

        Integer :: istage,i,uniti
        Real :: cumul
        Character :: dia_distr

        uniti=luinput(3)

c Initialise biological data

        stageprogressincrement=0.
        swimflaga=0
        swimflagr=0

        Do i=1,Maxp
         part_st(i)=1
         stageprogress(i)=Pstart_T(i)
        Enddo

        cumul=0.

c Read general biological information

        Read (uniti,101) true_no_stages
        if (true_no_stages.gt.10) then
          print *,'Stages greater than 10 not allowed!'
          stop
        endif
        Read (uniti,106) species
        Read (uniti,102) growthYN
        Read (uniti,102) mortalityYN
        Read (uniti,102) migrationYN
        Read (uniti,102) diapauseYN
        Read (uniti,102) swimming
        Read (uniti,109) swim_base
        Read (uniti,103) hatchlength
        Read (uniti,103) settlelength
        Read (uniti,108) mind
        Read (uniti,108) maxd
        Read (uniti,108) maxdi
        Read (uniti,108) mind_dia
        Read (uniti,108) maxd_dia
        Read (uniti,102) dia_distr
   70   Format (F8.3)

  101   Format (I2)
  102   Format (A1)
  103   Format (F4.1)
  104   Format (I5)
  105   Format (F7.5)
  106   Format (A7)
  107   Format (F8.2)
  108   Format (F7.2)
  109   Format (F5.2)
  110   Format (F7.4)

c Don't allow selection of diapause with any other behaviours

        if ((growthYN.eq.'y'.or.mortalityYN.eq.'y'.or.
     &       migrationYN.eq.'y').and.diapauseYN.eq.'y') then
          print *,'Diapause not allowed with any other behaviours'
          stop
        endif

c If absolute swimming speed, convert to m/s from cm/s
c Also set swim flags

        if (swimming.eq.'a') then
         swim_base=swim_base/100.
         swimflaga=1
        elseif (swimming.eq.'r') then
         swimflagr=1
        endif

c Loop round all stages

   10    Read (uniti,101,end=999) istage
         Read (uniti,102) stage_spec(istage)
         Read (uniti,107) stagelen(istage)
         Read (uniti,104) biocalcfreq(istage)
         Read (uniti,102) buoytyp(istage)
         Read (uniti,105) buoyancy(istage)
         Read (uniti,109) eggdiam(istage)
         Read (uniti,102) VMtyp(istage)
         Read (uniti,110) VMvel(istage)
         Read (uniti,102) diapause_cue(istage)
         Read (uniti,108) dia_down(istage)
         Read (uniti,108) dia_up(istage)
c Make sure dia_up Jday occurs after dia_down
         If (diapause_cue(istage).eq.'j'.and.dia_up(istage).lt.
     &       dia_down(istage)) dia_up(istage)=dia_up(istage)+365.

         goto 10

c In case the first stage is a larval stage
  999   if (stage_spec(1).eq.'l') then
         Do i=1,Maxp
          stageprogress(i)=hatchlength
         Enddo
        endif

c Change length of stages to cumulative seconds
c Buoyancy type v = buoyant velocity (already in m/s)
c Buoyancy type d = egg density (convert to kg/m^3)
c Change egg diameter from mm to m
c Convert vertical migration speeds to m/s

        Do i=1,true_no_stages
         cumul=cumul + stagelen(i)*24.*3600.
         stagelen(i)=cumul
         if (buoytyp(i).eq.'v'.or.buoytyp(i).eq.'V') then
          buoyancy(i)=buoyancy(i)
         elseif (buoytyp(i).eq.'d'.or.buoytyp(i).eq.'D') then
          buoyancy(i)=1000.*buoyancy(i)
         else
          buoyancy(i)=0.0
         endif
         eggdiam(i)=eggdiam(i)/1000.
         VMvel(i)=VMvel(i)/100.
c If species is passive set buotyp to v and buoyancy to zero
         if (species.eq.'passive') then
          buoytyp(i)='v'
          buoyancy(i)=0.
         endif
        Enddo

c Call subroutine to assign target depths to particles with diapause
c behaviour

        if (diapauseYN.eq.'y') then
         Call Diapause_depth (dia_distr)
        endif

        close(uniti)

        Return

        End subroutine

