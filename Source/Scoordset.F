        Subroutine Scoordset

c Subroutine to set up sigma levels in s coordinates consistent with
c POLCOMS
c N.B. z=sigo, zz=sig
c *** Location and name of file containing scoord information is
c *** hardwired and will need changing for future simulations
c *** STOP command included for now to alert the user

        use param

        implicit none

        Real*8 :: dsc,sds
        Real*8 :: sval(kb)
        Real*8 :: hc,cc,theta,bb,ffh,cs
        Integer :: icn, ierror, i, j, k, ii, iim1, iipn, iipnm1,unitst
        Integer :: nr
        character*100 :: filsco
        Real, allocatable, dimension(:,:) :: dsu,ds,ztmp,zztmp,zzvtmp
        Real, allocatable, dimension(:,:) :: zvtmp,dsv,dsuv
        Logical :: OK

c Allocate arrays

        allocate(dsu(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for dsu***"

        allocate(ds(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for ds***"

        allocate(ztmp(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for ztmp***"

        allocate(zztmp(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for zztmp***"

        allocate(zzvtmp(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for zzvtmp***"

        allocate(zvtmp(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for zvtmp***"

        allocate(dsv(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for dsv***"

        allocate(dsuv(Mx2,kb),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for dsuv***"

c Stop programme and print warning statement

        print *,'WARNING! Using scoord option - check filename in code'
        stop

c Evenly spaced S values

        sval(1) = -1.0d0
        sval(kb-1) = 0.0d0
        dsc = 1.0d0/(kb-2)
        do k=2,kb-2
          sval(k) = sval(k-1)+dsc
        enddo

c Read S coordinate parameters from 'scoord_params.dat'

        filsco = '~/scoord_params.dat'

        unitst=26
        Do
         Inquire(unit=unitst,opened=OK)
         If (.not.OK) exit
         unitst=unitst+1
        Enddo

        nr=unitst

        open(nr,file=filsco,form="formatted")
        read (nr,*) hc
        read (nr,*) cc
        read (nr,*) theta
        read (nr,*) bb
        close (nr)

        icn=0
        do j=1,jm
          do i=1,im
            icn=icn+1
            ztmp(icn,1) = -1.0d0
            do k=2,kb-2
              if(H(icn).gt.hc) then
                ffh=(H(icn)-hc)/H(icn)
                cs = (1.0d0-bb)*(sinh(theta*sval(k)))/sinh(theta)+
     &                bb*(tanh(theta*(sval(k)+0.5d0))
     &                -tanh(0.5d0*theta))/(2.*tanh(0.5d0*theta))
                ztmp(icn,k) = sval(k)+ffh*(cs-sval(k))
              else
                ztmp(icn,k) = sval(k)
              endif
            enddo
            ztmp(icn,kb-1) = 0.0d0
          enddo
        enddo

c Define coordinates on B points

        ds = 0.0
        icn=0
        do j=1,jm
          do i=1,im
            icn=icn+1
            do k=1,kb-2
              ds(icn,k) = ztmp(icn,k+1)-ztmp(icn,k)
            enddo
c Set surface level to ensure correct sum for ds
            sds = 0.0d0
            do k=1,kb-3
              sds = sds+ds(icn,k)
            enddo
            ds(icn,kb-2) = 1.0d0-sds
            dsu(icn,1) = ds(icn,1)
            do k=2,kb-2
              dsu(icn,k) = 0.5d0*(ds(icn,k)+ds(icn,k-1))
            enddo
            dsu(icn,kb-1) = ds(icn,kb-2)
            zztmp(icn,1) = -1.0d0+0.5d0*ds(icn,1)
            do k=2,kb-2
              zztmp(icn,k) = zztmp(icn,k-1)+dsu(icn,k)
            enddo
            zztmp(icn,kb-1) = 0.0d0
          enddo
        enddo

c Average coordinates onto U points

        do j=1,jm
          do i=1,im
            ii=(j-1)*im + i
            iim1=(j-1)*im + i-1
            iipn=j*im + i
            iipnm1=j*im + i-1
            if (i.eq.1) then
              iim1=ii
              iipnm1=iipn
            endif
            if (j.eq.jm) then
              iipnm1=iim1
              iipn=ii
            endif
            if (i.eq.1.and.j.eq.jm) iipnm1=ii
            zvtmp(ii,1) = -1.0d0
            if ((H(ii)+H(iim1)+H(iipnm1)+H(iipn)).ne.0.0) then
              do k=2,kb-2
                zvtmp(ii,k) = (
     &             H(ii)    *ztmp(ii,k)
     &            +H(iim1)  *ztmp(iim1,k)
     &            +H(iipnm1)*ztmp(iipnm1,k)
     &            +H(iipn)  *ztmp(iipn,k))/
     &            (H(ii)+H(iim1)+H(iipnm1)+H(iipn))
              enddo
            else
              do k=2,kb-2
                zvtmp(ii,k)=ztmp(ii,k)
              enddo
            endif
            zvtmp(ii,kb-1) = 0.0d0
            do k=1,kb-2
              dsv(ii,k) = zvtmp(ii,k+1)-zvtmp(ii,k)
            enddo
c Set surface level to ensure correct sum for DS
            sds=0.
            do k=1,kb-3
              sds=sds+dsv(ii,k)
            enddo
            dsv(ii,kb-2) = 1.0d0-sds
            dsuv(ii,1) = dsv(ii,1)
            dsuv(ii,kb-1) = dsv(ii,kb-2)
            do k=2,kb-2
              dsuv(ii,k) = 0.5d0*(dsv(ii,k)+dsv(ii,k-1))
            enddo
            zzvtmp(ii,1) = -1.0d0+0.5d0*dsv(ii,1)
            do k=2,kb-2
              zzvtmp(ii,k) = zzvtmp(ii,k-1)+dsuv(ii,k)
            enddo
            zzvtmp(ii,kb-1) = 0.0d0
            zzvtmp(ii,kb) = 0.0d0
          enddo
        enddo

c Now reverse z, zz and zzv so compatible with POM-based code

        do i=1,Mx2
          do k=1,kb-1
            z(i,k)=ztmp(i,kb-k)
            zz(i,k)=zztmp(i,kb-k)
            zzv(i,k)=zzvtmp(i,kb-k)
          enddo
        enddo

c Deallocate arrays

        deallocate(dsu)
        deallocate(ds)
        deallocate(ztmp)
        deallocate(zztmp)
        deallocate(zzvtmp)
        deallocate(zvtmp)
        deallocate(dsv)
        deallocate(dsuv)

        return

        End subroutine

