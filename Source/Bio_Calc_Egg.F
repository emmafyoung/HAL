      Subroutine Bio_Calc_Egg(i,T)

c This subroutine calculates temperature dependant growth and mortality
c for different egg stages.
c This is specific to plaice eggs and new routines will need to be
c written for other species
c Calculations bypassed on first call to prevent implementation of full
c day's growth/mortality in one Dt

        use param

        implicit none

        Integer, INTENT(IN) :: i, T

c Declare local variables

        Real :: dailymortality, mortality
        Real :: calpha(6),cbeta(6)

        calpha = (/ 5.186, 8.002, 12.819, 25.398, 29.880, 43.853 /)
        cbeta = (/ -1.612, -2.540, -4.098, -8.078, -9.313, -14.427 /)

c Mortality taken from Pepin (1991)

        if (mortalityYN.eq.'n'.or.(T-Pstart_T(i)).eq.0) goto 3002
        dailymortality=1.0-(0.03*EXP(0.18*Tempav_bio(i)))
        mortality=EXP(LOG(dailymortality)*biocalcfreq(part_st(i))
     &                /(86400.))
        partwt(i)=partwt(i)*mortality

 3002   continue

c Stage dependant growth taken from Fox et al. (2006)

        if ((T-Pstart_T(i)).eq.0) goto 3004
        stageprogressincrement=biocalcfreq(part_st(i))
        if (growthYN.eq.'n') goto 3003
        stageprogressincrement=biocalcfreq(part_st(i))*(
     1   LOG(calpha(part_st(i))+cbeta(part_st(i))*LOG(10.))/
     2   LOG(calpha(part_st(i))+cbeta(part_st(i))*LOG(Tempav_bio(i))))
 3003   stageprogress(i)=stageprogress(i)+stageprogressincrement
        if (stageprogress(i).gt.(Pstart_T(i)+stagelen(part_st(i))))
     1   then
         part_st(i)=part_st(i)+1
         if (stage_spec(part_st(i)).eq.'l')
     &        stageprogress(i)=hatchlength
        endif

 3004   Tempav_bio(i)=0

c Buoyancy terms

        eggd(i)=eggdiam(part_st(i))
        bwegg(i)=buoyancy(part_st(i))

        Return

        End subroutine

