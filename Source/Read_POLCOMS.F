        Subroutine Read_POLCOMS(velcount,timecount,endflag,switch)

c Read in velocity data for each output of the hydrodynamical model.
c Velocities are interpolated for the intervening timesteps.
c Note velocities are also staggered to be half a timestep out
c of phase with the position calculations.
c Setting values at level 1 (above sea surface) to those at level 2, 
c to avoid interpolation with zero values.

c *** Note, reading diffusivities currently commented out

        use param

        implicit none

        Integer, INTENT(IN) :: timecount,switch
        Integer, INTENT(INOUT) :: velcount,endflag

c Declare local variables

        Real :: Dtr, Dtr2
        Integer :: i, j, k, itmp, unitst
        Integer :: lo1, nrtape
        Integer, allocatable, dimension(:) :: isea,jsea,iusea,jusea
        Logical :: OK

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

c Zero Hzeta array before each calculation

        do i=1,mx2
         Hzeta(i) = 0.0
        enddo

        If (timecount.eq.0) then
         Open(nrtape,file=velfiles(velnum),form='unformatted',
     &        action='read')
         do 1 i=1,mx2
          Hzeta_a(i) = 0.0
          Hzeta_b(i) = 0.0
          sealevelrise(i) = 0.0
          mz(i)= 0
          mu(i)= 0
  1      continue

c Read header info.
         read (nrtape) ncc,nrr,nvl,nz
         if (ncc.ne.im.or.nrr.ne.jm.or.nvl.ne.kb) then
          print *,'Grid array sizes do not match those in included
     &          file!'
          print *,'Rows : ',nrr,jm
          print *,'Cols : ',ncc,im
          print *,'Vertical levels : ',nvl,kb
          endflag=-1
          goto 500
         endif

c Allocate arrays
         allocate(isea(nz))
         allocate(jsea(nz))
         read (nrtape) isea
         read (nrtape) jsea
         do i=1,nz
          iz(i)=(jm-jsea(i))*im + isea(i)
          mz(iz(i))=1
         enddo

c Read header info. for u and v arrays
         read (nrtape) ncc,nrr,nvl,nzu
         if (ncc.ne.im.or.nrr.ne.jm.or.nvl.ne.kb) then
          print *,'Grid array sizes do not match those in included
     &          file!'
          print *,'Rows : ',nrr,jm
          print *,'Cols : ',ncc,im
          print *,'Vertical levels : ',nvl,kb
          endflag=-1
          goto 500
         endif

c Allocate arrays
         allocate(iusea(nzu))
         allocate(jusea(nzu))
         read (nrtape) iusea
         read (nrtape) jusea
         do i=1,nzu
          iu(i)=(jm-jusea(i))*im + iusea(i)
          mu(iu(i))=1
         enddo

c Deallocate arrays
         deallocate(isea)
         deallocate(jsea)
         deallocate(iusea)
         deallocate(jusea)

        Endif

        If (velcount.eq.0) then

          read(nrtape,end=1000) lo1
          write (6,*) lo1
          read(nrtape) (Hzeta_a(iz(k)),k=1,nz)
          read(nrtape) ((U1(iu(k),j),j=kb-1,2,-1),k=1,nzu)
          read(nrtape) ((V1(iu(k),j),j=kb-1,2,-1),k=1,nzu)
          read(nrtape) ((W1(iz(k),j),j=kb-1,2,-1),k=1,nz)
          if (TSopt) then
           read(nrtape) ((T1(iz(k),j),j=kb-1,2,-1),k=1,nz)
           read(nrtape) ((Sa1(iz(k),j),j=kb-1,2,-1),k=1,nz)
          endif
c          read(nrtape) ((PKh1(iz(k),j),j=kb-1,2,-1),k=1,nz)
c          read(nrtape) ((PAh1(iz(k),j),j=kb-1,2,-1),k=1,nz)
          do 160 k=1,nz
            itmp=iz(k)
            Hzeta(itmp)=Hzeta_a(itmp)
            if (TSopt) then
             T1(itmp,1)=T1(itmp,2)
             Sa1(itmp,1)=Sa1(itmp,2)
            endif
c           PKh1(itmp,1)=PKh1(itmp,2)
c           PAh1(itmp,1)=PAh1(itmp,2)
  160     continue

c Read in the information for the next input too.

          read(nrtape,end=1000) lo1
          write (6,*) lo1
          read(nrtape) (Hzeta_b(iz(k)),k=1,nz)
          read(nrtape) ((U2(iu(k),j),j=kb-1,2,-1),k=1,nzu)
          read(nrtape) ((V2(iu(k),j),j=kb-1,2,-1),k=1,nzu)
          read(nrtape) ((W2(iz(k),j),j=kb-1,2,-1),k=1,nz)
          if (TSopt) then
           read(nrtape) ((T2(iz(k),j),j=kb-1,2,-1),k=1,nz)
           read(nrtape) ((Sa2(iz(k),j),j=kb-1,2,-1),k=1,nz)
          endif
c          read(nrtape) ((PKh2(iz(k),j),j=kb-1,2,-1),k=1,nz)
c          read(nrtape) ((PAh2(iz(k),j),j=kb-1,2,-1),k=1,nz)
          do 170 k=1,nz
            itmp=iz(k)
            W2(itmp,1)=W2(itmp,2)
            if (TSopt) then
             T2(itmp,1)=T2(itmp,2)
             Sa2(itmp,1)=Sa2(itmp,2)
            endif
c           PKh2(itmp,1)=PKh2(itmp,2)
c           PAh2(itmp,1)=PAh2(itmp,2)
  170     continue
          do 174 k=1,nzu
            itmp=iu(k)
            U2(itmp,1)=U2(itmp,2)
            V2(itmp,1)=V2(itmp,2)
  174     continue
          velcount=velcount+Dt

        Elseif (velcount.ge.Dt.and.velcount.lt.Dt2) then

         if (switch.eq.1) then
          Do 180 i=1,mx2
            if (mz(i).eq.0) goto 180
            Hzeta(i)=(((Dtr2-velcount)*Hzeta_a(i))+(float(velcount)*
     &          Hzeta_b(i)))/Dtr2
c Tidal height change in m/hour
            sealevelrise(i)=(Hzeta(i)-Hzeta_a(i))*3600./velcount
  180     continue
         endif
         velcount=velcount+Dt

        Elseif (velcount.eq.Dt2) then

          Do 190 i=1,mx2
            if (mz(i).eq.0) goto 190
            sealevelrise(i)=(Hzeta_b(i)-Hzeta_a(i))*3600./Dt2
            Hzeta_a(i)=Hzeta_b(i)
            Do 192 j=1,kb-1
             W1(i,j)=W2(i,j)
             if (TSopt) then
              T1(i,j)=T2(i,j)
              Sa1(i,j)=Sa2(i,j)
             endif
c             PKh1(i,j)=PKh2(i,j)
c             PAh1(i,j)=PAh2(i,j)
  192       continue
  190     continue
          Do 194 i=1,mx2
            if (mu(i).eq.0) goto 194
            Do 196 j=1,kb-1
             U1(i,j)=U2(i,j)
             V1(i,j)=V2(i,j)
  196       continue
  194     continue
  201     read(nrtape,end=1000) lo1
          write (6,*) lo1
          read(nrtape) (Hzeta_b(iz(k)),k=1,nz)
          read(nrtape) ((U2(iu(k),j),j=kb-1,2,-1),k=1,nzu)
          read(nrtape) ((V2(iu(k),j),j=kb-1,2,-1),k=1,nzu)
          read(nrtape) ((W2(iz(k),j),j=kb-1,2,-1),k=1,nz)
          if (TSopt) then
           read(nrtape) ((T2(iz(k),j),j=kb-1,2,-1),k=1,nz)
           read(nrtape) ((Sa2(iz(k),j),j=kb-1,2,-1),k=1,nz)
          endif
c          read(nrtape) ((PKh2(iz(k),j),j=kb-1,2,-1),k=1,nz)
c          read(nrtape) ((PAh2(iz(k),j),j=kb-1,2,-1),k=1,nz)

          Do 205 i=1,mx2
            if (mz(i).eq.0) goto 205
            Hzeta(i)=Hzeta_a(i)
            W2(i,1)=W2(i,2)
            if (TSopt) then
             T2(i,1)=T2(i,2)
             Sa2(i,1)=Sa2(i,2)
            endif
c           PKh2(i,1)=PKh2(i,2)
c           PAh2(i,1)=PAh2(i,2)
  205     continue
          Do 209 i=1,mx2
            if (mu(i).eq.0) goto 209
            U2(i,1)=U2(i,2)
            V2(i,1)=V2(i,2)
  209     continue
          velcount=Dt

        Endif

        Goto 500

 1000   print *,'End of velocity data'
        print *,'Closing first velocity file and opening next...'
        close (nrtape)
        velnum=velnum+1
        Open(nrtape,file=velfiles(velnum),form='unformatted',
     &       action='read')

c read header info

        read (nrtape,end=5000) ncc,nrr,nvl,nz
        allocate(isea(nz))
        allocate(jsea(nz))
        read (nrtape) isea
        read (nrtape) jsea
        read (nrtape) ncc,nrr,nvl,nzu
        allocate(iusea(nzu))
        allocate(jusea(nzu))
        read (nrtape) iusea
        read (nrtape) jusea

c Deallocate arrays
        deallocate(isea)
        deallocate(jsea)
        deallocate(iusea)
        deallocate(jusea)

c go back to main calculation

        goto 201

 5000   write (6,*) 'Problem with next file read! File : ',
     &               velfiles(velnum)
        endflag=-1

  500   return

        End subroutine

