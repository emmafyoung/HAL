        Subroutine Diapause_depth (dia_distr)

c Routine to assign initial target diapause depth for each particle
c using weighted random depth assignment within a defined depth range
c Also initialise dia_flag to 0 (surface phase)

        use param

        implicit none

        Character, INTENT(IN) :: dia_distr

        Integer :: ndep, ii, idx, iseed2, ierror
        Integer, allocatable, dimension(:) :: interd
        Real :: ran1
        Real*8 :: ixr, psum, prob, xvv
        Real, allocatable, dimension(:) :: pdist, pn, pdistn
        Data iseed2 /-1/

        ndep=maxd_dia-mind_dia-1

        allocate(interd(Npart), stat=ierror)
        if (ierror /=0) STOP "***Error allocating interd***"
        allocate(pdistn(Npart),stat=ierror)
        if (ierror /=0) STOP "***Error allocating pdistn***"
        allocate(pdist(ndep),stat=ierror)
        if (ierror /=0) STOP "***Error allocating pdist***"
        allocate(pn(ndep),stat=ierror)
        if (ierror /=0) STOP "***Error allocating pdn***"

        select case (dia_distr)

        case ('n')
c Random distribution

         do ii=1,Npart
          xvv=Ran(iseed2)
          pdistn(ii)=(ndep+1)*xvv
         enddo

         forall(ii=1:Npart) partd_dia(ii)=mind_dia+pdistn(ii)

        case ('y')
c Weighted random distribution within sin^2

         Do ii=1,ndep
          ixr=(ii*Pi)/(ndep+1)
          pdist(ii)=sin(ixr)*sin(ixr)
         Enddo

         psum=sum(pdist)

         forall(ii=1:ndep) pn(ii)=pdist(ii)/psum

c Loop around particles randomly assigning depth within distribution

         do ii=1,Npart

          xvv=Ran(iseed2)

          prob = 0
          do idx = 1, ndep
           prob = prob + pn( idx )   !! 0 < prob < 1
           if ( xvv <= prob ) exit
          enddo

          interd(ii)=idx

         enddo

         forall(ii=1:Npart) partd_dia(ii)=mind_dia+interd(ii)

        end select

        do ii=1,Npart
          write(54,*) (ii,partd_dia(ii))
        enddo

c Initialise dia_flag

        forall(ii=1:Npart) dia_flag(ii)=0

        deallocate(interd)
        deallocate(pdist)
        deallocate(pn)

        Return

        End subroutine


        
