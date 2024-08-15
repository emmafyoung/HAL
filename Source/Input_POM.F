        Subroutine Input_POM(endflag)

c Subroutine to input physical and model parameters - POM version

        use param

        implicit none

        Integer, INTENT(INOUT) :: endflag

        Integer :: nrbath,nrd,iinz,nzob,nlat,nu,nv,izob(nlim)
        Integer :: i,j,l,ierror,nsize,i2,j2
        Real :: chend,fiend,ch0,fi0,dx,dy
        Real, allocatable, dimension(:,:) :: zchiphi,chix,chiy,phiy

c Allocate arrays

        allocate(zchiphi(2,Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for zchiphi***"

        allocate(chix(3,Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for chix***"

        allocate(chiy(3,Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for chiy***"

        allocate(phiy(3,Mx2),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for phiy***"

        nrbath=luinput(5)
        nrd=luinput(6)

        read(nrbath,*) nrr,ncc,itot,iinz,nzob,nlat
        nsize=nrr*ncc
        if (ncc.ne.im.or.nrr.ne.jm) then
          print *,'Grid array sizes do not match those in included
     &          file!'
          print *,'Rows : ',nrr,jm
          print *,'Cols : ',ncc,im
          endflag=-1
          goto 999
        endif

        read(nrbath,*) ndch,ndfi
        read(nrbath,*) chst,chend,fist,fiend
        read(nrbath,*) ch0,fi0
        read(nrbath,*) nz,nu,nv
        read(nrbath,*) (iz1(i),i=1,nz)
        read(nrbath,*) (iz2(i),i=1,nz)
        read(nrbath,*) (iu1(i),i=1,nu)
        read(nrbath,*) (iu2(i),i=1,nu)
        read(nrbath,*) (iv1(i),i=1,nv)
        read(nrbath,*) (iv2(i),i=1,nv)
        read(nrbath,*) dx,dy
        read(nrbath,*) (h(i),i=1,itot)
        read(nrbath,*) (izob(i),i=1,nzob)
        read(nrbath,*)
        read(nrbath,*)
        read(nrbath,*)

c Calculate dch,dfi
c Extra read for curvilinear models

        if (curvi.eq.0) then
          dch=dx/deg2rad
          dfi=dy/deg2rad
        else
          read(nrbath,*) ((zchiphi(i,j),i=1,2),j=1,nsize)
          read(nrbath,*)
          read(nrbath,*)
          read(nrbath,*) ((chix(i,j),i=1,3),j=1,nsize)
          read(nrbath,*) ((chiy(i,j),i=1,3),j=1,nsize)
          read(nrbath,*) ((phiy(i,j),i=1,3),j=1,nsize)
          do j=1,nrr
            do i=1,ncc
              l=(j-1)*ncc + i
              lons(i,j)=zchiphi(1,l)
              lats(i,j)=zchiphi(2,l)
              dchcu(i,j)=chix(2,l)*dx
              dfic(i,j)=phiy(1,l)*dy
            enddo
          enddo
        endif

        close (nrbath)

        do 50 j=1,nz
          do 50 i=iz1(j),iz2(j)
            mz(i) = 1
  50    continue

c Read depth file and overwrite array read from gbath file

        do j=1,nrr
         do i=1,ncc
          l=(nrr-j)*ncc + i
          read (nrd,*) i2,j2,h(l)
         enddo
        enddo
        close (nrd)

c Deallocate arrays

        deallocate(zchiphi)
        deallocate(chix)
        deallocate(chiy)
        deallocate(phiy)

 999    return

        End subroutine

