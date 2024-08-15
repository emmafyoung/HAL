       Subroutine Bio_Calc_Larvae(i,i1,i2,i3,i4,x1,y1,pth,T)

c This routine calculates temperature dependent growth and mortality
c for plaice larvae. Other species will require new routines to be
c written for mortality and growth.

c Routine includes code for vertical migration and directed swimming
c Swimming rates are currently fixed: future developments should
c consider dependence of swimming speed on larval length.

        use param

        implicit none

        Integer, INTENT(IN) :: i,i1,i2,i3,i4,x1,y1,pth,T

        Real :: dailymortality, mortality
        Real :: rise
        Real :: p, q, alpha, beta, gamma, delta

        if (mortalityYN.eq.'n'.or.(T-Pstart_T(i)).eq.0) goto 3002

c Mortality formula taken from Pepin (1991)
c       dailymortality=1.0-(0.044*EXP(0.077*Tempav_bio(i)))

c Mortality formula taken from Nash (1997)
        dailymortality=(1.0)-(0.048)

        mortality=1-((-1*LOG(dailymortality))*
     &               biocalcfreq(part_st(i))/(86400.))
        partwt(i)=partwt(i)*mortality

 3002   continue

c Larval growth section

        if (growthYN.eq.'n'.or.(T-Pstart_T(i)).eq.0) goto 3004

c Growth Rate for larvae taken from Hovenkamp (1991)
c       stageprogressincrement=biocalcfreq(part_st(i))*
c     1 (0.00464*EXP(0.4319*Tempav_bio(i)))/86400.

c Growth Rate for larvae taken from Nash (1997)
        stageprogressincrement=(biocalcfreq(part_st(i))/86400.)*
     1   4./(122.81*EXP(0.092*(-1.0)*Tempav_bio(i)))

        stageprogress(i)=stageprogress(i)+stageprogressincrement
        if (stageprogress(i).gt.settlelength) then
          part_st(i)=part_st(i)+1
          print *,'End of Larval Stage. Weight = ',partwt(i)
        endif

 3004   Tempav_bio(i)=0

c Length dependent vertical migration
c Based on a variety of sources including personal communication with
c Geoff Arnold and larval vertical distributions published in the
c literature.

c Separate larvae into behavioural size groups

        if (migrationYN.eq.'n') goto 3006

        if (stageprogress(i).ge.(7.0)) then
         if (stageprogress(i).lt.(10.5)) then
           if (stageprogress(i).ge.(9.0)) then

c Find weighted distances for sea level rise at particle location
            p=Xelev(i)-x1
            q=Yelev(i)-y1
            alpha=p*q
            beta=p*(1-q)
            gamma=(1-p)*(1-q)
            delta=(1-p)*q
            rise=alpha*sealevelrise(i4)+beta*sealevelrise(i2)
     &           +gamma*sealevelrise(i1)+delta*sealevelrise(i3)

c Larvae between 9 and 10.5mm
c Light and tide-dependent migration

c Original migration rates
            if (rise.gt.(0.0).and.light(i).eq.'y')
     &        bwlarvae(i)=buoyancy(part_st(i))+0.001
            if (rise.gt.(0.0).and.light(i).eq.'n')
     &        bwlarvae(i)=buoyancy(part_st(i))+0.005
            if (rise.lt.(0.0).and.light(i).eq.'y')
     &        bwlarvae(i)=buoyancy(part_st(i))-0.005
            if (rise.lt.(0.0).and.light(i).eq.'n')
     &        bwlarvae(i)=buoyancy(part_st(i))-0.005

c Tidal only migration (for testing)
c           if (rise.gt.(0.0))
c     &         bwlarvae(i)=buoyancy(part_st(i))+0.005
c           if (rise.lt.(0.0))
c     &         bwlarvae(i)=buoyancy(part_st(i))-0.005

           Else

c Larvae between 7 and 9mm
c Light dependent migration for smaller larvae
c Type of DVM set in biology_spec.dat
c : 'd' for larvae at the surface during the day
c : 'n' for larvae at the surface during the night

            If (VMtyp(part_st(i)).eq.'d') then
             If (light(i).eq.'y') bwlarvae(i)=VMvel(part_st(i))
             If (light(i).eq.'n') bwlarvae(i)=VMvel(part_st(i))*(-1.)
            Elseif (VMtyp(part_st(i)).eq.'n') then
             If (light(i).eq.'y') bwlarvae(i)=VMvel(part_st(i))*(-1.)
             If (light(i).eq.'n') bwlarvae(i)=VMvel(part_st(i))
            Endif
           Endif

         Else

c Larvae longer than 10.5mm (settle out at 11mm)

           bwlarvae(i)=buoyancy(part_st(i))-0.01

         Endif
        Endif

c Directed Swimming
c Currently only used to hold late stage larvae in place on the seabed
c Needs to be light so larvae have a visual reference to orientate by
c Need to be close to the surface to be able to see it.
c The distance realistically is cm's or 10's of cm's but because of the
c model timestep and bathymetry 2m is used.
c The relative effects are somewhat untested
c Stage dependent speed should be converted to a length dependent
c formula
c Currently swim_speed of 0=no swimming,1=swimming into the current at
c the current speed, i.e no motion, etc.

 3006   swim_speed(i)=0.0
        if (swimming.eq.'n') goto 3008
        if (light(i).eq.'n') goto 3008
        if (ABS(pth-deppart(i)).gt.2) goto 3008
        if (stageprogress(i).ge.(7.5).and.stageprogress(i).lt.(8.5))
     &   swim_speed(i)=0.5*swim_base
        if (stageprogress(i).ge.(8.5).and.stageprogress(i).lt.(10.5))
     &   swim_speed(i)=1.1*swim_base
        if (stageprogress(i).ge.(10.5).and.stageprogress(i).lt.(11.5))
     &   swim_speed(i)=1.0*swim_base

 3008   Return

        End subroutine

