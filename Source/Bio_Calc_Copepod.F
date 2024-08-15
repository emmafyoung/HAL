       Subroutine Bio_Calc_Copepod(i,mod_jday)

c This routine calculates seasonal diapause behaviour for copepods.

        use param

        implicit none

        Integer, INTENT(IN) :: i, mod_jday

c Declare local variables

        Integer :: short_day

c Vertical diapause movement

        If (diapauseYN.eq.'n') goto 4001

c First set dia_flag to determine vertical movement
c Assumes single diapause cycle

        If (diapause_cue(part_st(i)).eq.'j') then
         If (mod_jday.ge.dia_down(part_st(i)).and.
     &       mod_jday.lt.dia_up(part_st(i)).and.
     &       dia_flag(i).ne.1) then
          dia_flag(i)=3
         Elseif (mod_jday.ge.dia_up(part_st(i)).and.
     &           dia_flag(i).ne.0) then
          dia_flag(i)=4
         Endif
        Elseif (diapause_cue(part_st(i)).eq.'d') then
c Specify shortest day depending on northern/southern hemisphere
c for daylight hours diapause cue
         If (Ydeg(i).lt.0) then
          short_day=172
         else
          short_day=355
         endif
         If (nhours(i).le.dia_down(part_st(i)).and.
     &       mod_jday.lt.short_day.and.
     &       dia_flag(i).ne.1) then
          dia_flag(i)=3
         Elseif (nhours(i).ge.dia_up(part_st(i)).and.
     &       mod_jday.ge.short_day.and.
     &       dia_flag(i).ne.0) then
          dia_flag(i)=4
         Endif
        Endif

        If (dia_flag(i).eq.0.or.dia_flag(i).eq.1) then
         bwlarvae(i)=0.
        Elseif (dia_flag(i).eq.3) then
         bwlarvae(i)=VMvel(part_st(i))*(-1.)
        Elseif (dia_flag(i).eq.4) then
         bwlarvae(i)=VMvel(part_st(i))
        Endif

 4001   Return

        End subroutine

