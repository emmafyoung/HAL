        Subroutine Read_NEMO_netcdf(velcount,timecount,endflag,
     &switch,fcount,numdt,ncidu,ncidv,ncidw,ncidts,ncidi,varidu,varidv,
     &varidw,varidt,varids,varidi,varidiu,varidiv,numlons,numlats,
     &numlevs)

c Read in velocity data for each output of the hydrodynamical model.
c Also options for reading temperature and salinity, and ice data.
c Variables now interpolated in Tinterpvel for the intervening
c timesteps.
c NEMO fields start at bottom left so flip latitude.

        use netcdf
        use param

        implicit none

        Integer, INTENT(IN) :: timecount,switch
        Integer, INTENT(INOUT) :: velcount,endflag,fcount,numdt,ncidu,
     &ncidv,ncidw,ncidts,ncidi,varidu,varidv,varidw,varidt,varids,
     &varidi,varidiu,varidiv,numlons,numlats,numlevs

c Declare local variables

        Real Dtr, Dtr2
        Integer lo1, nrtape,ii,status,chlons
        Integer i,j,k,ierror
        Character (len=100) :: fsub,fnameu,fnamev,fnamew,fnamets,fnamei
        Integer, dimension(4) :: dimIDs

        Real, allocatable, dimension(:,:,:) :: Ua,Ub,Va,Vb,Wa,Wb
        Real, allocatable, dimension(:,:,:) :: Ta,Tb,Sa,Sb
        Real, allocatable, dimension(:,:) :: Ia,Ib,Iua,Iub,Iva,Ivb
        Real, allocatable, dimension(:,:,:,:) :: vtmp
        Real, allocatable, dimension(:,:,:) :: itmp

c Allocate arrays

        allocate(Ua(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ua***"

        allocate(Ub(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ub***"

        allocate(Va(im,jm,kb),stat=ierror)
       if (ierror /= 0) STOP "***Not enough memory for Va***"

        allocate(Vb(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Vb***"

        allocate(Wa(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Wa***"

        allocate(Wb(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Wb***"

        allocate(Ta(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ta***"

        allocate(Tb(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Tb***"

        allocate(Sa(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Sa***"

        allocate(Sb(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Sb***"

        allocate(Ia(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ia***"

        allocate(Ib(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ib***"

        allocate(Iua(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Iua***"

        allocate(Iub(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Iub***"

        allocate(Iva(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Iva***"

        allocate(Ivb(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for Ivb***"

        allocate(vtmp(im+2,jm,kb,1),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for vtmp***"

        allocate(itmp(im+2,jm,1),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for itmp***"

        Dtr2=float(Dt2)
        Dtr=float(Dt)

        fsub=velfiles(velnum)
        fnameu=fsub(1:len_trim(fsub))//'U.nc'
        fnamev=fsub(1:len_trim(fsub))//'V.nc'
        fnamew=fsub(1:len_trim(fsub))//'W.nc'
        fnamets=fsub(1:len_trim(fsub))//'T.nc'
        fnamei=fsub(1:len_trim(fsub))//'I.nc'

        If (velcount.eq.0) then

         status=nf90_open(fnameu,nf90_nowrite,ncidu)
         call Check(status,endflag,'opening input u file')
         if (endflag.lt.0) goto 500
         print *,'Opened velfile: ',fnameu
         if (newname) then
          status=nf90_inq_varid(ncidu,"uo",varidu)
         else
          status=nf90_inq_varid(ncidu,"vozocrtx",varidu)
         endif
         status=nf90_inquire_variable(ncidu,varidu,dimids=dimIDs)
         status=nf90_inquire_dimension(ncidu,dimIDs(1),len=numlons)
         status=nf90_inquire_dimension(ncidu,dimIDs(2),len=numlats)
         status=nf90_inquire_dimension(ncidu,dimIDs(3),len=numlevs)
         status=nf90_inquire_dimension(ncidu,dimIDs(4),len=numdt)
         print *,'Dimensions of velocity fields: ',numlons,numlats,
     &numlevs,numdt

         if (nemo.and.circumpolar) then
          chlons=numlons-2
         else
          chlons=numlons
         endif

         if (chlons.ne.im.or.numlats.ne.jm.or.numlevs.ne.kb) then
          print *,'Velocity array sizes do not match those in included
     &          file!'
          print *,'Rows : ',numlats,jm
          print *,'Cols : ',chlons,im
          print *,'Vertical levels : ',numlevs,kb
          endflag=-1
          goto 500
         endif

         do 1 i=1,im
          do 2 j=1,jm
           do 3 k=1,kb
            ua(i,j,k) = 0.0
            va(i,j,k) = 0.0
            wa(i,j,k) = 0.0
            ub(i,j,k) = 0.0
            vb(i,j,k) = 0.0
            wb(i,j,k) = 0.0
  3        continue
  2       continue
  1      continue

         if (TSopt) then
          do i=1,im
           do j=1,jm
            do k=1,kb
             ta(i,j,k) = 0.0
             sa(i,j,k) = 0.0
             tb(i,j,k) = 0.0
             sb(i,j,k) = 0.0
            enddo
           enddo
          enddo
         endif

         if (Iopt) then
          do i=1,im
           do j=1,jm
            ia(i,j) = 0.0
            ib(i,j) = 0.0
            Iua(i,j) = 0.0
            Iub(i,j) = 0.0
            Iva(i,j) = 0.0
            Ivb(i,j) = 0.0
           enddo
          enddo
         endif

         status=nf90_open(fnamev,nf90_nowrite,ncidv)
         call Check(status,endflag,'opening input v file')
         if (endflag.lt.0) goto 500
         if (newname) then
          status=nf90_inq_varid(ncidv,"vo",varidv)
         else
          status=nf90_inq_varid(ncidv,"vomecrty",varidv)
         endif
         status=nf90_open(fnamew,nf90_nowrite,ncidw)
         call Check(status,endflag,'opening input w file')
         if (endflag.lt.0) goto 500
         if (newname) then
          status=nf90_inq_varid(ncidw,"wo",varidw)
         else
          status=nf90_inq_varid(ncidw,"vovecrtz",varidw)
         endif

         if (TSopt) then
          status=nf90_open(fnamets,nf90_nowrite,ncidts)
          call Check(status,endflag,'opening input t file')
          if (endflag.lt.0) goto 500
          if (newname) then
           status=nf90_inq_varid(ncidts,"potemp",varidt)
           status=nf90_inq_varid(ncidts,"salin",varids)
          else
           status=nf90_inq_varid(ncidts,"votemper",varidt)
           status=nf90_inq_varid(ncidts,"vosaline",varids)
          endif
         endif

         if (Iopt) then
          status=nf90_open(fnamei,nf90_nowrite,ncidi)
          call Check(status,endflag,'opening input i file')
          if (endflag.lt.0) goto 500
          if (newname) then
           status=nf90_inq_varid(ncidi,"ice_pres",varidi)
           status=nf90_inq_varid(ncidi,"uice_ipa",varidiu)
           status=nf90_inq_varid(ncidi,"vice_ipa",varidiv)
          else
           status=nf90_inq_varid(ncidi,"siconc",varidi)
           status=nf90_inq_varid(ncidi,"uice",varidiu)
           status=nf90_inq_varid(ncidi,"vice",varidiv)
          endif
         endif

c Reshape arrays after reading in if using circumpolar nemo

         if (nemo.and.circumpolar) then

          status=nf90_get_var(ncidu,varidu,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) ua(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) ua(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(ua)) ua=0.
          status=nf90_get_var(ncidv,varidv,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) va(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) va(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(va)) va=0.
          status=nf90_get_var(ncidw,varidw,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) wa(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) wa(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(wa)) wa=0.

          if (TSopt) then
           status=nf90_get_var(ncidts,varidt,vtmp,
     &start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           forall(j=1:jm,k=1:kb) ta(1:1012,j,k)=vtmp(431:1442,j,k,1)
           forall(j=1:jm,k=1:kb) ta(1013:1440,j,k)=vtmp(3:430,j,k,1)
           where (isnan(ta)) ta=0.
           status=nf90_get_var(ncidts,varids,vtmp,
     &start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           forall(j=1:jm,k=1:kb) sa(1:1012,j,k)=vtmp(431:1442,j,k,1)
           forall(j=1:jm,k=1:kb) sa(1013:1440,j,k)=vtmp(3:430,j,k,1)
           where (isnan(sa)) sa=0.
          endif

          if (Iopt) then
           status=nf90_get_var(ncidi,varidi,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) ia(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) ia(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(ia)) ia=0.
           status=nf90_get_var(ncidi,varidiu,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) iua(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) iua(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(iua)) iua=0.
           status=nf90_get_var(ncidi,varidiv,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) iva(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) iva(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(iva)) iva=0.
          endif

          print *,'Successfully read first input files'

         else

          status=nf90_get_var(ncidu,varidu,ua,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(ua)) ua=0.
          status=nf90_get_var(ncidv,varidv,va,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(va)) va=0.
          status=nf90_get_var(ncidw,varidw,wa,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(wa)) wa=0.

          if (TSopt) then
           status=nf90_get_var(ncidts,varidt,ta,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           where (isnan(ta)) ta=0.
           status=nf90_get_var(ncidts,varids,sa,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           where (isnan(sa)) sa=0.
          endif

          if (Iopt) then
           status=nf90_get_var(ncidi,varidi,ia,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(ia)) ia=0.
           status=nf90_get_var(ncidi,varidiu,iua,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(iua)) iua=0.
           status=nf90_get_var(ncidi,varidiv,iva,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(iva)) iva=0.
          endif

         endif

         fcount=fcount+1

c If single 5-day mean files, close current files and open next ones

         if (numdt.eq.1) then
           status=nf90_close(ncidu)
           status=nf90_close(ncidv)
           status=nf90_close(ncidw)
           if (TSopt) then
            status=nf90_close(ncidts)
           endif
           if (Iopt) then
            status=nf90_close(ncidi)
           endif
           velnum=velnum+1
           fsub=velfiles(velnum)
           fnameu=fsub(1:len_trim(fsub))//'U.nc'
           fnamev=fsub(1:len_trim(fsub))//'V.nc'
           fnamew=fsub(1:len_trim(fsub))//'W.nc'
           fnamets=fsub(1:len_trim(fsub))//'T.nc'
           fnamei=fsub(1:len_trim(fsub))//'I.nc'
           status=nf90_open(fnameu,nf90_nowrite,ncidu)
           call Check(status,endflag,'opening input u file')
           if (endflag.lt.0) goto 500
           print *,'Opened velfile: ',fnameu
           if (newname) then
            status=nf90_inq_varid(ncidu,"uo",varidu)
           else
            status=nf90_inq_varid(ncidu,"vozocrtx",varidu)
           endif
           status=nf90_open(fnamev,nf90_nowrite,ncidv)
           call Check(status,endflag,'opening input v file')
           if (endflag.lt.0) goto 500
           if (newname) then
            status=nf90_inq_varid(ncidv,"vo",varidv)
           else
            status=nf90_inq_varid(ncidv,"vomecrty",varidv)
           endif
           status=nf90_open(fnamew,nf90_nowrite,ncidw)
           call Check(status,endflag,'opening input w file')
           if (endflag.lt.0) goto 500
           if (newname) then
            status=nf90_inq_varid(ncidw,"wo",varidw)
           else
            status=nf90_inq_varid(ncidw,"vovecrtz",varidw)
           endif
           if (TSopt) then
            status=nf90_open(fnamets,nf90_nowrite,ncidts)
            call Check(status,endflag,'opening input t file')
            if (endflag.lt.0) goto 500
            if (newname) then
             status=nf90_inq_varid(ncidts,"potemp",varidt)
             status=nf90_inq_varid(ncidts,"salin",varids)
            else
             status=nf90_inq_varid(ncidts,"votemper",varidt)
             status=nf90_inq_varid(ncidts,"vosaline",varids)
            endif
           endif
           if (Iopt) then
            status=nf90_open(fnamei,nf90_nowrite,ncidi)
            call Check(status,endflag,'opening input i file')
            if (endflag.lt.0) goto 500
            if (newname) then
             status=nf90_inq_varid(ncidi,"ice_pres",varidi)
             status=nf90_inq_varid(ncidi,"uice_ipa",varidiu)
             status=nf90_inq_varid(ncidi,"vice_ipa",varidiv)
            else
             status=nf90_inq_varid(ncidi,"siconc",varidi)
             status=nf90_inq_varid(ncidi,"uice",varidiu)
             status=nf90_inq_varid(ncidi,"vice",varidiv)
            endif
           endif
           fcount=1
         endif

c Read in the information for the next input too.

         if (nemo.and.circumpolar) then

          status=nf90_get_var(ncidu,varidu,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) ub(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) ub(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(ub)) ub=0.
          status=nf90_get_var(ncidv,varidv,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) vb(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) vb(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(vb)) vb=0.
          status=nf90_get_var(ncidw,varidw,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) wb(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) wb(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(wb)) wb=0.

          if (TSopt) then
           status=nf90_get_var(ncidts,varidt,vtmp,
     &start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           forall(j=1:jm,k=1:kb) tb(1:1012,j,k)=vtmp(431:1442,j,k,1)
           forall(j=1:jm,k=1:kb) tb(1013:1440,j,k)=vtmp(3:430,j,k,1)
           where (isnan(tb)) tb=0.
           status=nf90_get_var(ncidts,varids,vtmp,
     &start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           forall(j=1:jm,k=1:kb) sb(1:1012,j,k)=vtmp(431:1442,j,k,1)
           forall(j=1:jm,k=1:kb) sb(1013:1440,j,k)=vtmp(3:430,j,k,1)
           where (isnan(sb)) sb=0.
          endif

          if (Iopt) then
           status=nf90_get_var(ncidi,varidi,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) ib(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) ib(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(ib)) ib=0.
           status=nf90_get_var(ncidi,varidiu,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) iub(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) iub(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(iub)) iub=0.
           status=nf90_get_var(ncidi,varidiv,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) ivb(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) ivb(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(ivb)) ivb=0.
          endif

          print *,'Successfully read second input files'

         else

          status=nf90_get_var(ncidu,varidu,ub,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(ub)) ub=0.
          status=nf90_get_var(ncidv,varidv,vb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(vb)) vb=0.
          status=nf90_get_var(ncidw,varidw,wb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(wb)) wb=0.

          if (TSopt) then
           status=nf90_get_var(ncidts,varidt,tb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           where (isnan(tb)) tb=0.
           status=nf90_get_var(ncidts,varids,sb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           where (isnan(sb)) sb=0.
          endif

          if (Iopt) then
           status=nf90_get_var(ncidi,varidi,ib,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(ib)) ib=0.
           status=nf90_get_var(ncidi,varidiu,iub,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(iub)) iub=0.
           status=nf90_get_var(ncidi,varidiv,ivb,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(ivb)) ivb=0.
          endif

         endif
         fcount=fcount+1

c Convert to snake array and flip latitude

         do 174 i=1,im
           do 175 j=1,jm
            ii=i+(jm-j)*im
            do 176 k=1,kb
             U1(ii,k)=Ua(i,j,k)
             U2(ii,k)=Ub(i,j,k)
             V1(ii,k)=Va(i,j,k)
             V2(ii,k)=Vb(i,j,k)
             W1(ii,k)=Wa(i,j,k)
             W2(ii,k)=Wb(i,j,k)
             if (TSopt) then
              T1(ii,k)=Ta(i,j,k)
              T2(ii,k)=Tb(i,j,k)
              Sa1(ii,k)=Sa(i,j,k)
              Sa2(ii,k)=Sb(i,j,k)
             endif
  176       continue
            if (Iopt) then
             Ice1(ii)=Ia(i,j)
             Ice2(ii)=Ib(i,j)
             Uice1(ii)=Iua(i,j)
             Uice2(ii)=Iub(i,j)
             Vice1(ii)=Iva(i,j)
             Vice2(ii)=Ivb(i,j)
            endif
  175      continue
  174    continue
         velcount=velcount+Dt

c If single 5-day mean files, close current files and open next ones

         if (numdt.eq.1) then
           status=nf90_close(ncidu)
           status=nf90_close(ncidv)
           status=nf90_close(ncidw)
           if (TSopt) then
            status=nf90_close(ncidts)
           endif
           if (Iopt) then
            status=nf90_close(ncidi)
           endif
           velnum=velnum+1
           fsub=velfiles(velnum)
           fnameu=fsub(1:len_trim(fsub))//'U.nc'
           fnamev=fsub(1:len_trim(fsub))//'V.nc'
           fnamew=fsub(1:len_trim(fsub))//'W.nc'
           fnamets=fsub(1:len_trim(fsub))//'T.nc'
           fnamei=fsub(1:len_trim(fsub))//'I.nc'
           status=nf90_open(fnameu,nf90_nowrite,ncidu)
           call Check(status,endflag,'opening input u file')
           if (endflag.lt.0) goto 500
           print *,'Opened velfile: ',fnameu
           if (newname) then
            status=nf90_inq_varid(ncidu,"uo",varidu)
           else
            status=nf90_inq_varid(ncidu,"vozocrtx",varidu)
           endif
           status=nf90_open(fnamev,nf90_nowrite,ncidv)
           call Check(status,endflag,'opening input v file')
           if (endflag.lt.0) goto 500
           if (newname) then
            status=nf90_inq_varid(ncidv,"vo",varidv)
           else
            status=nf90_inq_varid(ncidv,"vomecrty",varidv)
           endif
           status=nf90_open(fnamew,nf90_nowrite,ncidw)
           call Check(status,endflag,'opening input w file')
           if (endflag.lt.0) goto 500
           if (newname) then
            status=nf90_inq_varid(ncidw,"wo",varidw)
           else
            status=nf90_inq_varid(ncidw,"vovecrtz",varidw)
           endif
           if (TSopt) then
            status=nf90_open(fnamets,nf90_nowrite,ncidts)
            call Check(status,endflag,'opening input t file')
            if (endflag.lt.0) goto 500
            if (newname) then
             status=nf90_inq_varid(ncidts,"potemp",varidt)
             status=nf90_inq_varid(ncidts,"salin",varids)
            else
             status=nf90_inq_varid(ncidts,"votemper",varidt)
             status=nf90_inq_varid(ncidts,"vosaline",varids)
            endif
           endif
           if (Iopt) then
            status=nf90_open(fnamei,nf90_nowrite,ncidi)
            call Check(status,endflag,'opening input i file')
            if (endflag.lt.0) goto 500
            if (newname) then
             status=nf90_inq_varid(ncidi,"ice_pres",varidi)
             status=nf90_inq_varid(ncidi,"uice_ipa",varidiu)
             status=nf90_inq_varid(ncidi,"vice_ipa",varidiv)
            else
             status=nf90_inq_varid(ncidi,"siconc",varidi)
             status=nf90_inq_varid(ncidi,"uice",varidiu)
             status=nf90_inq_varid(ncidi,"vice",varidiv)
            endif
           endif
           fcount=1
         endif

        Elseif (velcount.ge.Dt.and.velcount.lt.Dt2) then

          velcount=velcount+Dt

        Elseif (velcount.eq.Dt2) then

         print *,'In Dt2 Velcount: ',velcount

         if (nemo.and.circumpolar) then

          status=nf90_get_var(ncidu,varidu,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          call Check(status,endflag,'reading input u file')
          if (endflag.lt.0) goto 500
          forall(j=1:jm,k=1:kb) ub(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) ub(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(ub)) ub=0.
          status=nf90_get_var(ncidv,varidv,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) vb(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) vb(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(vb)) vb=0.
          status=nf90_get_var(ncidw,varidw,vtmp,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          forall(j=1:jm,k=1:kb) wb(1:1012,j,k)=vtmp(431:1442,j,k,1)
          forall(j=1:jm,k=1:kb) wb(1013:1440,j,k)=vtmp(3:430,j,k,1)
          where (isnan(wb)) wb=0.

          if (TSopt) then
           status=nf90_get_var(ncidts,varidt,vtmp,
     &start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           forall(j=1:jm,k=1:kb) tb(1:1012,j,k)=vtmp(431:1442,j,k,1)
           forall(j=1:jm,k=1:kb) tb(1013:1440,j,k)=vtmp(3:430,j,k,1)
           where (isnan(tb)) tb=0.
           status=nf90_get_var(ncidts,varids,vtmp,
     &start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           forall(j=1:jm,k=1:kb) sb(1:1012,j,k)=vtmp(431:1442,j,k,1)
           forall(j=1:jm,k=1:kb) sb(1013:1440,j,k)=vtmp(3:430,j,k,1)
           where (isnan(sb)) sb=0.
          endif

          if (Iopt) then
           status=nf90_get_var(ncidi,varidi,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) ib(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) ib(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(ib)) ib=0.
           status=nf90_get_var(ncidi,varidiu,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) iub(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) iub(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(iub)) iub=0.
           status=nf90_get_var(ncidi,varidiv,itmp,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           forall(j=1:jm) ivb(1:1012,j)=itmp(431:1442,j,1)
           forall(j=1:jm) ivb(1013:1440,j)=itmp(3:430,j,1)
           where (isnan(ivb)) ivb=0.
          endif

          print *,'Successfully read next input files'

         else

          status=nf90_get_var(ncidu,varidu,ub,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(ub)) ub=0.
          status=nf90_get_var(ncidv,varidv,vb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(vb)) vb=0.
          status=nf90_get_var(ncidw,varidw,wb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
          where (isnan(wb)) wb=0.

          if (TSopt) then
           status=nf90_get_var(ncidts,varidt,tb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           where (isnan(tb)) tb=0.
           status=nf90_get_var(ncidts,varids,sb,start=(/1,1,1,fcount/),
     &count=(/numlons,numlats,numlevs,1/))
           where (isnan(sb)) sb=0.
          endif

          if (Iopt) then
           status=nf90_get_var(ncidi,varidi,ib,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(ib)) ib=0.
           status=nf90_get_var(ncidi,varidiu,iub,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(iub)) iub=0.
           status=nf90_get_var(ncidi,varidiv,ivb,start=(/1,1,fcount/),
     &count=(/numlons,numlats,1/))
           where (isnan(ivb)) ivb=0.
          endif

         endif

         fcount=fcount+1

c Convert to snake array and flip latitude

         do 209 i=1,im
           do 210 j=1,jm
            ii=i+(jm-j)*im
            do 211 k=1,kb
             U1(ii,k)=U2(ii,k)
             U2(ii,k)=Ub(i,j,k)
             V1(ii,k)=V2(ii,k)
             V2(ii,k)=Vb(i,j,k)
             W1(ii,k)=W2(ii,k)
             W2(ii,k)=Wb(i,j,k)
             if (TSopt) then
              T1(ii,k)=T2(ii,k)
              T2(ii,k)=Tb(i,j,k)
              Sa1(ii,k)=Sa2(ii,k)
              Sa2(ii,k)=Sb(i,j,k)
             endif
  211       continue
            if (Iopt) then
             Ice1(ii)=Ice2(ii)
             Ice2(ii)=Ib(i,j)
             Uice1(ii)=Uice2(ii)
             Uice2(ii)=Iub(i,j)
             Vice1(ii)=Vice2(ii)
             Vice2(ii)=Ivb(i,j)
            endif
  210      continue
  209    continue
         velcount=Dt
         print *,'In Dt2 reset Velcount: ',velcount

c If single 5-day mean files, or reached the end of the file,
c close current files and open next ones

         if (numdt.eq.1.or.(fcount-1).eq.numdt) then
           status=nf90_close(ncidu)
           status=nf90_close(ncidv)
           status=nf90_close(ncidw)
           if (TSopt) then
            status=nf90_close(ncidts)
           endif
           if (Iopt) then
            status=nf90_close(ncidi)
           endif
           velnum=velnum+1
           fsub=velfiles(velnum)
           fnameu=fsub(1:len_trim(fsub))//'U.nc'
           fnamev=fsub(1:len_trim(fsub))//'V.nc'
           fnamew=fsub(1:len_trim(fsub))//'W.nc'
           fnamets=fsub(1:len_trim(fsub))//'T.nc'
           fnamei=fsub(1:len_trim(fsub))//'I.nc'
           status=nf90_open(fnameu,nf90_nowrite,ncidu)
           call Check(status,endflag,'opening input u file')
           if (endflag.lt.0) goto 500
           print *,'Opened velfile: ',fnameu
           if (newname) then
            status=nf90_inq_varid(ncidu,"uo",varidu)
           else
            status=nf90_inq_varid(ncidu,"vozocrtx",varidu)
           endif
           status=nf90_open(fnamev,nf90_nowrite,ncidv)
           call Check(status,endflag,'opening input v file')
           if (endflag.lt.0) goto 500
           if (newname) then
            status=nf90_inq_varid(ncidv,"vo",varidv)
           else
            status=nf90_inq_varid(ncidv,"vomecrty",varidv)
           endif
           status=nf90_open(fnamew,nf90_nowrite,ncidw)
           call Check(status,endflag,'opening input w file')
           if (endflag.lt.0) goto 500
           if (newname) then
            status=nf90_inq_varid(ncidw,"wo",varidw)
           else
            status=nf90_inq_varid(ncidw,"vovecrtz",varidw)
           endif
           if (TSopt) then
            status=nf90_open(fnamets,nf90_nowrite,ncidts)
            call Check(status,endflag,'opening input t file')
            if (endflag.lt.0) goto 500
            if (newname) then
             status=nf90_inq_varid(ncidts,"potemp",varidt)
             status=nf90_inq_varid(ncidts,"salin",varids)
            else
             status=nf90_inq_varid(ncidts,"votemper",varidt)
             status=nf90_inq_varid(ncidts,"vosaline",varids)
            endif
           endif
           if (Iopt) then
            status=nf90_open(fnamei,nf90_nowrite,ncidi)
            call Check(status,endflag,'opening input i file')
            if (endflag.lt.0) goto 500
            if (newname) then
             status=nf90_inq_varid(ncidi,"ice_pres",varidi)
             status=nf90_inq_varid(ncidi,"uice_ipa",varidiu)
             status=nf90_inq_varid(ncidi,"vice_ipa",varidiv)
            else
             status=nf90_inq_varid(ncidi,"siconc",varidi)
             status=nf90_inq_varid(ncidi,"uice",varidiu)
             status=nf90_inq_varid(ncidi,"vice",varidiv)
            endif
           endif
           fcount=1
         endif
        Endif

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
        deallocate(Ia)
        deallocate(Ib)
        deallocate(Iua)
        deallocate(Iub)
        deallocate(Iva)
        deallocate(Ivb)
        deallocate(vtmp)
        deallocate(itmp)

        return

        End subroutine

