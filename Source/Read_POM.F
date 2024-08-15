        Subroutine Read_POM(velcount,timecount,endflag,switch)

c Read in velocity data for each output of the hydrodynamical model.
c Velocities are interpolated for the intervening timesteps.
c Note velocities are also staggered to be half a timestep out
c of phase with the position calculations.

        use param

        implicit none

        Integer, INTENT(IN) :: timecount,switch
        Integer, INTENT(INOUT) :: velcount,endflag

c Declare local variables

        Real*8 :: ztmp(kb)
        Real :: sig(maxterms),f(maxterms),uv(maxterms)
        Real :: Dtr2, Dtr
        Integer :: lo, lo1, ngap, nrtape, maxt
        Integer :: ierror, i, j, k, l, readcount
        Integer :: nsize, unitst
        Real, allocatable, dimension(:,:) :: Ua,Ub,Va,Vb,Wa,Wb
        Real, allocatable, dimension(:,:) :: Ta,Tb,Sa,Sb
        Real, allocatable, dimension(:,:) :: PKha,PKhb,PAha,PAhb
        Logical :: OK

c Allocate arrays

        allocate(Ua(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ua***"

        allocate(Ub(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ub***"

        allocate(Va(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Va***"

        allocate(Vb(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Vb***"

        allocate(Wa(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Wa***"

        allocate(Wb(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Wb***"

        allocate(Ta(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ta***"

        allocate(Tb(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Tb***"

        allocate(Sa(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Sa***"

        allocate(Sb(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Sb***"

        allocate(PKha(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PKha***"

        allocate(PKhb(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PKhb***"

        allocate(PAha(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PAha***"

        allocate(PAhb(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for PAhb***"

c Find a free unit to assign to the input file

        unitst=26
        Do
         Inquire(unit=unitst,opened=OK)
         If (.not.OK) exit
         unitst=unitst+1
        Enddo

        nrtape=unitst

        Dtr2=float(Dt2)
        Dtr=float(Dt)
        nsize=nrr*ncc

        If (timecount.eq.0) then
         Open(nrtape,file=velfiles(velnum),form='unformatted',
     &        action='read')
         do 1 i=1,mx2
          Hzeta_a(i) = 0.0
          Hzeta_b(i) = 0.0
          sealevelrise(i) = 0.0
          mz(i)= 0
          do 2 j=1,kb
            ua(i,j) = 0.0
            ub(i,j) = 0.0
            va(i,j) = 0.0
            vb(i,j) = 0.0
            ta(i,j) = 0.0
            tb(i,j) = 0.0
            sa(i,j) = 0.0
            sb(i,j) = 0.0
            PKha(i,j) = 0.0
            PKhb(i,j) = 0.0
            PAha(i,j) = 0.0
            PAhb(i,j) = 0.0
  2       continue
  1      continue
c Do preliminary read to get information on output gap in POM
         read(nrtape) Dt2,maxt,nvl
         if (nvl.ne.kb) then
          print *,'Levels in tape file different to included file!'
          print *,'Levels : ',nvl,kb
          endflag=-1
          goto 500
         endif
         read(nrtape)
         read(nrtape) (ztmp(l),l=1,nvl)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
c Convert ztmp to 2D array
c Calculate z levels from zz
         do i=1,Mx2
           do k=1,nvl
             zz(i,k)=ztmp(k)
           enddo
           if (nvl.gt.1) then
            z(i,1)=0.
            do k=2,nvl
              z(i,k)=2*zz(i,k-1) - z(i,k-1)
            enddo
           else
            z(i,1)=0.5
           endif
         enddo
         if (nvl.eq.1) nvl=2
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape) lo1
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape) lo
         ngap = lo - lo1
c Read POM tape file
         rewind(nrtape)
         read(nrtape) Dt2,Maxt,nvl
         dt2=dt2*ngap
         print *,'Tape file output time interval = ',dt2
         if (mod(int(Dt2),Dt).ne.0) then
          print *,'Timestep must be a divisor of hydrodynamic model 
     &                  output'
          print *,'Dt2, ngap: ',Dt2,ngap
          print *,'Dt: ',Dt
          endflag=-1
          goto 500
         endif
         read(nrtape) (sig(i),f(i),uv(i),i=1,Maxt)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape)
         read(nrtape) (cspz(i),i=1,itot)
        Endif

        If (velcount.eq.0) then

          read(nrtape,end=1000) ((Hzeta_a(l),l=iz1(k),iz2(k)),k=1,nz)
          read(nrtape) (((Ua(l,j),l=iu1(k),iu2(k)),k=1,nzu),j=1,nvl-1)
          read(nrtape) (((Va(l,j),l=iv1(k),iv2(k)),k=1,nzu),j=1,nvl-1)
          read(nrtape) (((Wa(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((Ta(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((Sa(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((PKha(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((PAha(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          do 160 k=1,nz
            do 160 l=iz1(k),iz2(k)
              Hzeta(l)=Hzeta_a(l)
              do 162 i=1,nvl
                Tem(l,i)=Ta(l,i)
                Sal(l,i)=Sa(l,i)
                PKh(l,i)=PKha(l,i)
                PAh(l,i)=PAha(l,i)
  162         continue
  160     continue
          read(nrtape) lo1

c Read in the information for the next input too.
          readcount=readcount+1
          read(nrtape,end=1000) ((Hzeta_b(l),l=iz1(k),iz2(k)),k=1,nz)
          read(nrtape) (((Ub(l,j),l=iu1(k),iu2(k)),k=1,nzu),j=1,nvl-1)
          read(nrtape) (((Vb(l,j),l=iv1(k),iv2(k)),k=1,nzu),j=1,nvl-1)
          read(nrtape) (((Wb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((Tb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((Sb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((PKhb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((PAhb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) lo1
          do 170 i=1,nsize
           if (mz(i).eq.0) goto 170
            do 175 j=1,nvl
             U(i,j)=(((Dtr2-(Dtr/2.))*Ua(i,j))+((Dtr/2.)*Ub(i,j)))/Dtr2
             V(i,j)=(((Dtr2-(Dtr/2.))*Va(i,j))+((Dtr/2.)*Vb(i,j)))/Dtr2
             W(i,j)=(((Dtr2-(Dtr/2.))*Wa(i,j))+((Dtr/2.)*Wb(i,j)))/Dtr2
  175       continue
  170     continue
          velcount=velcount+Dt

        Elseif (velcount.ge.Dt.and.velcount.lt.Dt2) then
         If (switch.eq.1) then
          Do 180 i=1,nsize
            if (mz(i).eq.0) goto 180
            Do 185 j=1,nvl
             U(i,j)=(((Dtr2-(Dtr/2.)-velcount)*Ua(i,j))+
     &               ((velcount+(Dtr/2.))*Ub(i,j)))/Dtr2
             V(i,j)=(((Dtr2-(Dtr/2.)-velcount)*Va(i,j))+
     &               ((velcount+(Dtr/2.))*Vb(i,j)))/Dtr2
             W(i,j)=(((Dtr2-(Dtr/2.)-velcount)*Wa(i,j))+
     &               ((velcount+(Dtr/2.))*Wb(i,j)))/Dtr2
             Tem(i,j)=(((Dtr2-velcount)*Ta(i,j))+(float(velcount)*
     &                   Tb(i,j)))/Dtr2
             Sal(i,j)=(((Dtr2-velcount)*Sa(i,j))+(float(velcount)*
     &                   Sb(i,j)))/Dtr2
             PKh(i,j)=(((Dtr2-velcount)*PKha(i,j))+(float(velcount)*
     &                   PKhb(i,j)))/Dtr2
             PAh(i,j)=(((Dtr2-velcount)*PAha(i,j))+(float(velcount)*
     &                   PAhb(i,j)))/Dtr2
  185       Continue
            Hzeta(I)=(((Dtr2-velcount)*Hzeta_a(I))+(float(velcount)*
     &                  Hzeta_b(I)))/Dtr2
c Tidal height change in m/hour
            sealevelrise(I)=(Hzeta(I)-Hzeta_a(I))*3600./velcount
  180     continue
         Endif
         velcount=velcount+Dt

        Elseif (velcount.eq.Dt2) then
          Do 190 i=1,nsize
            if (mz(i).eq.0) goto 190
            sealevelrise(I)=(Hzeta_b(i)-Hzeta_a(i))*3600./Dt2
            Hzeta_a(i)=Hzeta_b(i)
            Do 191 j=1,nvl
             Ua(i,j)=Ub(i,j)
             Va(i,j)=Vb(i,j)
             Wa(i,j)=Wb(i,j)
             Ta(i,j)=Tb(i,j)
             Sa(i,j)=Sb(i,j)
             PKha(i,j)=PKhb(i,j)
             PAha(i,j)=PAhb(i,j)
  191       continue
  190     continue
          readcount=readcount+1
  201     read(nrtape,end=1000) ((Hzeta_b(l),l=iz1(k),iz2(k)),k=1,nz)
          read(nrtape) (((Ub(l,j),l=iu1(k),iu2(k)),k=1,nzu),j=1,nvl-1)
          read(nrtape) (((Vb(l,j),l=iv1(k),iv2(k)),k=1,nzu),j=1,nvl-1)
          read(nrtape) (((Wb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((Tb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((Sb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((PKhb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) (((PAhb(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
          read(nrtape) lo1
          Do 195 i=1,nsize
            if (mz(i).eq.0) goto 195
            Hzeta(i)=Hzeta_a(i)
            Do 196 j=1,nvl
             U(i,j)=(((Dtr2-(Dtr/2.))*Ua(i,j))+((Dt/2.)*Ub(i,j)))/Dtr2
             V(i,j)=(((Dtr2-(Dtr/2.))*Va(i,j))+((Dt/2.)*Vb(i,j)))/Dtr2
             W(i,j)=(((Dtr2-(Dtr/2.))*Wa(i,j))+((Dt/2.)*Wb(i,j)))/Dtr2
             Tem(i,j)=Ta(i,j)
             Sal(i,j)=Sa(i,j)
             PKh(i,j)=PKha(i,j)
             PAh(i,j)=PAha(i,j)
  196       continue
  195     continue
          velcount=Dt
        Endif

        Goto 500

 1000   print *,'End of velocity data'
        print *,'Closing first velocity file and opening next...'
        close (nrtape)
        velnum=velnum+1
        Open(nrtape,file=velfiles(velnum),status='old',
     &       form='unformatted',action='read')

c get to start of data and dummy read first set of data as it is
c duplicated...

        Read (nrtape)
        Read (nrtape)
        Read (nrtape)
        Read (nrtape)
        Read (nrtape)
        Read (nrtape)
        Read (nrtape)
        read(nrtape,end=5000) ((Hzeta_a(l),l=iz1(k),iz2(k)),k=1,nz)
        read(nrtape) (((Ua(l,j),l=iu1(k),iu2(k)),k=1,nzu),j=1,nvl-1)
        read(nrtape) (((Va(l,j),l=iv1(k),iv2(k)),k=1,nzu),j=1,nvl-1)
        read(nrtape) (((Wa(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
        read(nrtape) (((Ta(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
        read(nrtape) (((Sa(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
        read(nrtape) (((PKha(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
        read(nrtape) (((PAha(l,j),l=iz1(k),iz2(k)),k=1,nz),j=1,nvl-1)
        read(nrtape) lo1

c go back to main calculation

        goto 201

 5000   write (6,*) 'Problem with next file read! File : ',
     &               velfiles(velnum)
        endflag=-1

  500   continue

c Deallocate arrays

        deallocate(Ua)
        deallocate(Ub)
        deallocate(Va)
        deallocate(Vb)
        deallocate(Wa)
        deallocate(Wb)
        deallocate(Ta)
        deallocate(Tb)
        deallocate(Sa)
        deallocate(Sb)
        deallocate(PKha)
        deallocate(PKhb)
        deallocate(PAha)
        deallocate(PAhb)

        return

        End subroutine

