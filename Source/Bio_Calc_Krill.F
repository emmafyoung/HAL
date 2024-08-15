        Subroutine Bio_Calc_Krill(i)

c This routine calculates light-dependent vertical migration for krill.
c Also includes an option for directed swimming.

        use param

        implicit none

        Integer, INTENT(IN) :: i

c Vertical migration

        If (migrationYN.eq.'n') goto 3001
        If (VMtyp(part_st(i)).eq.'d') then
         If (light(i).eq.'y') bwlarvae(i)=VMvel(part_st(i))
         If (light(i).eq.'n') bwlarvae(i)=VMvel(part_st(i))*(-1.)
        Elseif (VMtyp(part_st(i)).eq.'n') then
         If (light(i).eq.'y') bwlarvae(i)=VMvel(part_st(i))*(-1.)
         If (light(i).eq.'n') bwlarvae(i)=VMvel(part_st(i))
        Endif

c Directed Swimming

 3001   swim_speed(i)=0.0
        if (swimming.eq.'n') goto 4001
        swim_speed(i)=swim_base

 4001   Return

        End subroutine

