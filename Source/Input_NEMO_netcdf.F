        Subroutine Input_NEMO_netcdf(endflag)

c Subroutine to input physical and model parameters - NEMO netcdf
c version

        use netcdf
        use param

        implicit none

        Integer, INTENT(INOUT) :: endflag

        Integer :: icn,idum,i,j,k,ierror,varid,status
        Integer, allocatable, dimension(:,:) :: mbathy,kmt,kmu,kmv,
     &          tmaskutil,umaskutil,vmaskutil
        Real :: ndep,edep,sdep,wdep,cdep,shallowpt
        Real, allocatable, dimension(:) :: vlat
        Real, allocatable, dimension(:,:) :: glamt,gphit,glamu,gphiv,hin
        Real, allocatable, dimension(:,:) :: e1t,e2t
        Real, allocatable, dimension(:,:) :: gradu,gradv
        Real, allocatable, dimension(:,:,:) :: gdept,gdepw,e3t
        Real, allocatable, dimension(:,:,:) :: temp,temp2

        write (6,*) 'In Input_NEMO_netcdf'

c Allocate temporary arrays

        allocate(vlat(jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for vlat***"
        allocate(glamt(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for glamt***"
        allocate(gphit(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for gphit***"
        allocate(glamu(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for glamu***"
        allocate(gphiv(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for gphiv***"
        allocate(e1t(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for e1t***"
        allocate(e2t(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for e2t***"
        allocate(mbathy(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for mbathy***"
        allocate(tmaskutil(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for tmaskutil***"
        allocate(umaskutil(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for umaskutil***"
        allocate(vmaskutil(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for vmaskutil***"
        allocate(hin(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for hin***"
        allocate(gradu(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for gradu***"
        allocate(gradv(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for gradv***"
        allocate(kmt(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for kmt***"
        allocate(kmu(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for kmu***"
        allocate(kmv(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for kmv***"
        allocate(gdept(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for gdept***"
        allocate(gdepw(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for gdepw***"
        allocate(temp(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for temp***"
        allocate(temp2(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for temp2***"
        allocate(e3t(im,jm,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for e3t***"

        write(6,*) 'Allocated arrays, size of temp: ',SIZE(temp)
        write(6,*) 'Global allocated arrays, size of nemolev: ',
     &SIZE(nemolev)

c Input the grid information

c Read required variables from mesh_hgr
        status=nf90_inq_varid(luinput(5),"glamt",varid)
        status=nf90_get_var(luinput(5),varid,glamt)
        call Check(status,endflag,'glamt read')
        status=nf90_inq_varid(luinput(5),"gphit",varid)
        status=nf90_get_var(luinput(5),varid,gphit)
        call Check(status,endflag,'gphit read')
        status=nf90_inq_varid(luinput(5),"glamu",varid)
        status=nf90_get_var(luinput(5),varid,glamu)
        call Check(status,endflag,'glamu read')
        status=nf90_inq_varid(luinput(5),"gphiv",varid)
        status=nf90_get_var(luinput(5),varid,gphiv)
        call Check(status,endflag,'gphiv read')
        if (curvi.eq.1) then
         status=nf90_inq_varid(luinput(5),"e1t",varid)
         status=nf90_get_var(luinput(5),varid,e1t)
         call Check(status,endflag,'e1t read')
         status=nf90_inq_varid(luinput(5),"e2t",varid)
         status=nf90_get_var(luinput(5),varid,e2t)
         call Check(status,endflag,'e2t read')
        endif

c Read required variables from mesh_zgr
        if (nemo.and.circumpolar) then
c Different variable names in circumpolar nemo files
         status=nf90_inq_varid(luinput(6),"gdept",varid)
         status=nf90_get_var(luinput(6),varid,gdept)
         call Check(status,endflag,'gdept read')
         status=nf90_inq_varid(luinput(6),"gdepw",varid)
         status=nf90_get_var(luinput(6),varid,gdepw)
         call Check(status,endflag,'gdepw read')
         status=nf90_inq_varid(luinput(6),"e3t",varid)
         status=nf90_get_var(luinput(6),varid,e3t)
         call Check(status,endflag,'e3t read')
        else
         status=nf90_inq_varid(luinput(6),"gdept_0",varid)
         status=nf90_get_var(luinput(6),varid,gdept)
         call Check(status,endflag,'gdept_0 read')
         status=nf90_inq_varid(luinput(6),"gdepw_0",varid)
         status=nf90_get_var(luinput(6),varid,gdepw)
         call Check(status,endflag,'gdepw_0 read')
         status=nf90_inq_varid(luinput(6),"e3t_0",varid)
         status=nf90_get_var(luinput(6),varid,e3t)
         call Check(status,endflag,'e3t_0 read')
        endif
        status=nf90_inq_varid(luinput(6),"mbathy",varid)
        status=nf90_get_var(luinput(6),varid,mbathy)
        call Check(status,endflag,'mbathy read')

c Read required variables from mask
        status=nf90_inq_varid(luinput(7),"tmaskutil",varid)
        status=nf90_get_var(luinput(7),varid,tmaskutil)
        call Check(status,endflag,'tmaskutil read')
        status=nf90_inq_varid(luinput(7),"umaskutil",varid)
        status=nf90_get_var(luinput(7),varid,umaskutil)
        call Check(status,endflag,'umaskutil read')
        status=nf90_inq_varid(luinput(7),"vmaskutil",varid)
        status=nf90_get_var(luinput(7),varid,vmaskutil)
        call Check(status,endflag,'vmaskutil read')

c Get ncc, nrr, nvl from size(gdept)

        ncc=size(gdept,1)
        nrr=size(gdept,2)
        nvl=size(gdept,3)

        write(6,*) 'ncc',ncc
        write(6,*) 'nrr',nrr
        write(6,*) 'nvl',nvl
        if (ncc.ne.im.or.nrr.ne.jm.or.nvl.ne.kb) then
          print *,'Grid array sizes do not match those in included
     &          file!'
          print *,'Rows : ',nrr,jm
          print *,'Cols : ',ncc,im
          print *,'Vertical levels : ',nvl,kb
          endflag=-1
          goto 999
        endif

c If circumpolar set and the number of columns doesn't match the
c 1/4 deg model, stop!

        if (circumpolar.and.(ncc.ne.1440)) then
         print *,'Circumpolar option currently only available for the
     &         1/4 deg N205 model!'
         endflag=-1
         goto 999
        endif

c Extract top left t-point coordinates from glamt, gphit
        chst=minval(glamt)
        write(6,*) 'chst',chst
        fist=maxval(gphit)
        write(6,*) 'fist',fist

c Following 4 sections of code for curvilinear in y, not x

c Calculate longitude resolution (grid cells per degree)
        ndch=nint(1/(glamt(2,1)-glamt(1,1)))
        write(6,*) 'ndch',ndch

c Calculate dch
        dch=1./ndch
        write(6,*) 'dch',dch

c Extract 1D latitude array and change origin to NW corner
        forall(j=1:jm) lat(j)=gphit(1,jm-j+1)
        write(6,*) 'lat(1),lat(jm)',lat(1),lat(jm)

c Extract 1D dlat array (centred on t-points, origin NW corner)
        forall(j=1:jm) vlat(j)=gphiv(1,jm-j+1)
        forall(j=1:jm-1) dlat(j)=vlat(j)-vlat(j+1)
        dlat(jm)=dlat(jm-1)
        write(6,*) 'dlat(1),dlat(jm)',dlat(1),dlat(jm)

c Following 3 sections of code for curvilinear in x and y

c Create 2D lon, lat arrays referenced to top left
        forall (i=1:im, j=1:jm)
         lons(i,jm-j+1)=glamt(i,j)
         lats(i,jm-j+1)=gphit(i,j)
        end forall

        if (curvi.eq.1) then

c Create 2D dchi array, centred on t-points, origin NW corner

         forall (i=1:im,j=1:jm)
          dchcu(i,jm-j+1)=(180.*e1t(i,j))/(Pi*Radius_E*
     &                     cos(Deg2Rad*gphit(i,j)))
         end forall

c Create 2D dfi array, centred on t-points, origin NW corner

         forall (i=1:im,j=1:jm)
          dfic(i,jm-j+1)=(180.*e2t(i,j))/(Pi*Radius_E)
         end forall

        endif

c Create cspz array
        icn=0
        do j=1,jm
         do i=1,im
          icn=icn+1
          cspz(icn)=cos(lats(i,j)*Deg2Rad)
         enddo
        enddo

c Create 2D and 3D depth arrays (origin top left)
        do i=1,im
         do j=1,jm
          do k=1,kb
           gdepw(i,j,k)=sum(e3t(i,j,1:k));
          enddo
          if (mbathy(i,j).eq.0) then
           hin(i,jm-j+1)=0
          else
           hin(i,jm-j+1)=gdepw(i,j,(mbathy(i,j)))
          endif
         enddo
        enddo
        write(6,*) 'Created hin array'

c Create 2d gradient flag arrays for directed swimming option
        forall(i=1:im, j=1:jm)
         gradu(i,j)=0.
         gradv(i,j)=0.
        end forall
        do j=2,jm-1
         do i=2,im-1
          ndep=hin(i,j-1)
          edep=hin(i+1,j)
          sdep=hin(i,j+1)
          wdep=hin(i-1,j)
          cdep=hin(i,j)
          shallowpt=min(ndep,edep,sdep,wdep,cdep)
          if (ndep.eq.shallowpt) gradv(i,j)=1.
          if (edep.eq.shallowpt) gradu(i,j)=1.
          if (sdep.eq.shallowpt) gradv(i,j)=-1.
          if (wdep.eq.shallowpt) gradu(i,j)=-1.
         enddo
        enddo
        write(6,*) 'Created ugrad and vgrad array'

c Create 3D arrays of t and w-pt depth (ensure origin is top left)
c Mask points below land
        forall(i=1:im, j=1:jm, k=1:kb) temp(i,jm-j+1,k)=gdept(i,j,k)
        forall(i=1:im, j=1:jm) temp2(i,jm-j+1,1)=0
        forall(i=1:im, j=1:jm, k=2:kb) temp2(i,jm-j+1,k)=gdepw(i,j,k-1)
        do j=1,jm
         do i=1,im
          if (mbathy(i,j).le.kb-1) then
           temp(i,jm-j+1,mbathy(i,j)+1:kb)=-1
          endif
          if (mbathy(i,j).le.kb-2) then
           temp2(i,jm-j+1,mbathy(i,j)+2:kb)=-1
          endif
         enddo
        enddo
        write(6,*) 'Created temp and temp2 arrays'

c Extract mask arrays (ensure origin is top left)
        forall(i=1:im, j=1:jm)
         kmt(i,j)=tmaskutil(i,jm-j+1)
         kmu(i,j)=umaskutil(i,jm-j+1)
         kmv(i,j)=vmaskutil(i,jm-j+1)
        end forall
        write(6,*) 'Extracted mask arrays'

c Convert 2D arrays to 1D snake arrays
        icn=0
        do j=1,jm
          do i=1,im
           icn=icn + 1
           H(icn)=hin(i,j)
           swimu(icn)=gradu(i,j)
           swimv(icn)=gradv(i,j)
           if (kmt(i,j).eq.0) then
            mz(icn)=0
           else
            mz(icn)=1
           endif
           if (kmu(i,j).eq.0) then
            mu(icn)=0
           else
            mu(icn)=1
           endif
           if (kmv(i,j).eq.0) then
            mv(icn)=0
           else
            mv(icn)=1
           endif
           do k=1,kb
            nemolev(icn,k)=temp(i,j,k)
            nemowlev(icn,k)=temp2(i,j,k)
           enddo
          enddo
        enddo
        write(6,*) 'Created snake arrays'

c Deallocate arrays

        deallocate(glamt)
        deallocate(gphit)
        deallocate(glamu)
        deallocate(gphiv)
        deallocate(mbathy)
        deallocate(tmaskutil)
        deallocate(umaskutil)
        deallocate(vmaskutil)
        deallocate(hin)
        deallocate(kmt)
        deallocate(kmu)
        deallocate(kmv)
        deallocate(gdept)
        deallocate(gdepw)
        deallocate(temp)
        deallocate(temp2)
        deallocate(e3t)

c Close source netcdf files

        status=nf90_close(luinput(5))
        status=nf90_close(luinput(6))
        status=nf90_close(luinput(7))

 999    return

        End subroutine

