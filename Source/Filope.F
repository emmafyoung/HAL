        Subroutine Filope (warmst,filenm,numnamesi)

        use netcdf
        use param

        implicit none

c Routine to open input and output files

        Character, INTENT(INOUT) :: warmst
        Character*50, INTENT(IN) :: filenm
        Integer, INTENT(OUT) :: numnamesi

        Integer :: i,ii,lendat,lenin,lenout,iotype,lenid,lenin0,status
        Integer :: numnames,unit1,ncid
        Character*90 :: iname,iname0
        Character*90 :: idname,fname
        Character*256 :: indir
        Integer :: oldmode,endflag

c Specify unit number to start assigning from

        unit1=13

c Input files...

        open (1,file=filenm)
        read (1,'(a)') iname
        read (1,'(a)') indir

c Find the length of the directory string

        do i=1,len(indir)
         if (indir(i:i).eq.' ') then
          lenin=i-1
          goto 10
         endif
        enddo
        lenin=len(indir)
  10    continue

c Open files

        read (1,*) numnames
        write (6,*) 'Input file names: '
        do i=1,numnames
         read (1,'(i2,1x,a)') iotype,iname
         write (6,'(i2,1x,a)') iotype,iname(1:len_trim(iname))
         fname=indir(1:lenin)//'/'//iname
         if (iotype.eq.1) then
          open(unit1,file=fname)
         elseif (iotype.eq.0) then
          open(unit1,file=fname,form='unformatted')
         elseif (iotype.eq.2) then
          open(unit1,file=iname)
         elseif (iotype.eq.3) then
          open(unit1,file=iname,form='unformatted')
         elseif (iotype.eq.4) then
          status=nf90_open(fname,nf90_nowrite,ncid)
          call Check(status,endflag,'opening input netcdf file')
         elseif (iotype.eq.5) then
          status=nf90_open(iname,nf90_nowrite,ncid)
          call Check(status,endflag,'opening input netcdf file')
         endif
         if (iname(1:3).eq.'res') warmst='y'
         if (iotype.lt.4) then
          luinput(i)=unit1
          unit1=unit1+1
         else
          luinput(i)=ncid
          if (ncid.eq.unit1) unit1=unit1+1
         endif
        enddo

        numnamesi=numnames

c Output files...

        read (1,'(a)') iname
        read (1,*) numnames
        read (1,'(a)') Idname

c Find the length of the directory string

        do i=1,len(Idname)
         if (Idname(i:i).eq.' ') then
          lenid=i-1
          goto 20
         endif
        enddo
        lenid=len(Idname)
  20    continue

c Open files

        write (6,*) 'Output file names: '
        do ii=1,numnames
         read (1,'(i2,1x,a)') iotype,iname0
         do i=1,len(iname0)
          if (iname0(i:i).eq.' ') then
           lenin0=i-1
           goto 30
          endif
         enddo
         lenin0=len(iname0)
  30     continue
         iname=Idname(1:lenid)//'/'//iname0
         write(*,'(I2,1x,A)') iotype,iname0(1:len_trim(iname0))
c Only allow netcdf output format
         if (ii.eq.1.and.iotype.lt.4) then
          print *,'Only netcdf output format allowed'
          STOP 1
         endif
         if (iotype.eq.0) then
          open(unit1,file=iname,form='unformatted')
         elseif(iotype.eq.1) then
          open(unit1,file=iname)
         elseif (iotype.eq.2) then
          open(unit1,file=iname0)
         elseif (iotype.eq.3) then
          open(unit1,file=iname0,form='unformatted')
         elseif (iotype.eq.4) then
          status=nf90_create(iname,nf90_netcdf4,ncid)
          call Check(status,endflag,'opening output netcdf file')
c Set to no fill mode
          status=nf90_set_fill(ncid,nf90_nofill,oldmode)
          call Check(status,endflag,'setting nofill in output file')
         elseif (iotype.eq.5) then
          status=nf90_create(iname0,nf90_netcdf4,ncid)
          call Check(status,endflag,'opening output netcdf file')
c Set to no fill mode
          status=nf90_set_fill(ncid,nf90_nofill,oldmode)
          call Check(status,endflag,'setting nofill in output file')
         endif
         if (iotype.lt.4) then
          luoutput(ii)=unit1
          unit1=unit1+1
         else
          luoutput(ii)=ncid
         endif
        enddo

        print *,'Allocated input unit numbers : ',luinput
        print *,'Allocated output unit numbers : ',luoutput
        close(1)

        if (endflag.lt.0) STOP 1

        return
        end

