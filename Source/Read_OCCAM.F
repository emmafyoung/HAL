        Subroutine Read_OCCAM(velcount,timecount,endflag,switch)

c Read in velocity data for each output of the hydrodynamical model.
c Velocities now interpolated in Tinterpvel for the intervening
c timesteps.
c OCCAM fields start at bottom left so flip latitude.

        use param

        implicit none

        Integer, INTENT(IN) :: timecount,switch
        Integer, INTENT(INOUT) :: velcount,endflag

c Declare local variables

        Real :: Dtr, Dtr2
        Integer :: lo1, nrtape, ii, i, j, k, ierror, unitst
        Real, allocatable, dimension(:,:,:) :: Ua,Ub,Va,Vb
        Logical :: OK

c Allocate arrays

        allocate(Ua(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ua***"

        allocate(Ub(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ub***"

        allocate(Va(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Va***"

        allocate(Vb(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Vb***"

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

        If (timecount.eq.0) then
         Open(nrtape,file=velfiles(velnum),form='unformatted',
     &        action='read')
         do 1 i=1,im
          do 2 j=1,jm
           do 3 k=1,kb
            ua(i,j,k) = 0.0
            va(i,j,k) = 0.0
  3        continue
  2       continue
  1      continue

c Read header info.
         read (nrtape) ncc,nrr,nvl
         if (ncc.ne.im.or.nrr.ne.jm.or.nvl.ne.kb) then
          print *,'Grid array sizes do not match those in included
     &          file!'
          print *,'Rows : ',nrr,jm
          print *,'Cols : ',ncc,im
          print *,'Vertical levels : ',nvl,kb
          endflag=-1
          goto 500
         endif
        Endif

        If (velcount.eq.0) then

          read(nrtape,end=1000) lo1
          write (6,*) lo1
          read(nrtape) (((Ua(i,j,k),i=1,im),j=1,jm),k=1,kb)
          read(nrtape) (((Va(i,j,k),i=1,im),j=1,jm),k=1,kb)

c Read in the information for the next input too.

          read(nrtape,end=1000) lo1
          write (6,*) lo1
          read(nrtape) (((Ub(i,j,k),i=1,im),j=1,jm),k=1,kb)
          read(nrtape) (((Vb(i,j,k),i=1,im),j=1,jm),k=1,kb)

c Convert to snake array
c Convert velocities from cm/s to m/s

          do 174 i=1,im
           do 175 j=1,jm
            ii=i+(jm-j)*im
            if (mu(ii).eq.0) goto 175
            do 176 k=1,kb
             U1(ii,k)=Ua(i,j,k)/100.
             U2(ii,k)=Ub(i,j,k)/100.
             V1(ii,k)=Va(i,j,k)/100.
             V2(ii,k)=Vb(i,j,k)/100.
  176       continue
  175      continue
  174     continue
          velcount=velcount+Dt

        Elseif (velcount.ge.Dt.and.velcount.lt.Dt2) then

          velcount=velcount+Dt

        Elseif (velcount.eq.Dt2) then

  201     read(nrtape,end=1000) lo1
          write (6,*) lo1
          read(nrtape) (((Ub(i,j,k),i=1,im),j=1,jm),k=1,kb)
          read(nrtape) (((Vb(i,j,k),i=1,im),j=1,jm),k=1,kb)

          Do 209 i=1,im
           Do 210 j=1,jm
            ii=i+(jm-j)*im
            if (mu(ii).eq.0) goto 210
            Do 211 k=1,kb
             U1(ii,k)=U2(ii,k)
             U2(ii,k)=Ub(i,j,k)/100.
             V1(ii,k)=V2(ii,k)
             V2(ii,k)=Vb(i,j,k)/100.
  211       continue
  210      continue
  209     continue
         velcount=Dt

        Endif

        Goto 500

 1000   print *,'End of velocity data'
        print *,'Closing first velocity file and opening next...'
        close (nrtape)
        velnum=velnum+1
        Open(nrtape,file=velfiles(velnum),status='old',
     &       form='unformatted',action='read')

c Read header info.
        read (nrtape) ncc,nrr,nvl

c get to start of data and dummy read first set of data as it is
c duplicated...

        read(nrtape,end=5000) lo1
        read(nrtape) (((Ub(i,j,k),i=1,im),j=1,jm),k=1,kb)
        read(nrtape) (((Vb(i,j,k),i=1,im),j=1,jm),k=1,kb)

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

        return

        End subroutine

