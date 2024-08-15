        Subroutine Input_phys(endflag)

c Subroutine to input physical and model parameters

        use param

        implicit none

        Integer, INTENT(INOUT) :: endflag

        Real :: hmin,lattmp
        Integer :: ind,idum,ierror,i,j,k,icn
        Real, allocatable, dimension(:,:) :: hin,kmt,kmu
        Integer :: uniti,uniti2,uniti3,uniti4

c Allocate arrays

        allocate(hin(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for hin***"

        allocate(kmt(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for kmt***"

        allocate(kmu(im,jm),stat=ierror)
        if (ierror /= 0) STOP "***Not enough memory for kmu***"

        uniti=luinput(6)
        uniti2=luinput(7)
        uniti3=luinput(5)
        uniti4=luinput(8)

c Input the grid information

        read (uniti3,10) ncc
        read (uniti3,10) nrr
        read (uniti3,10) nvl
  10    format(I5)
        if (ncc.ne.im.or.nrr.ne.jm.or.nvl.ne.kb) then
          print *,'Grid array sizes do not match those in included
     &          file!'
          print *,'Rows : ',nrr,jm
          print *,'Cols : ',ncc,im
          print *,'Vertical levels : ',nvl,kb
          endflag=-1
          goto 999
        endif
        read (uniti3,20) hmin
  20    format(F5.1)
        read (uniti3,10) ndch
        read (uniti3,10) ndfi
        read (uniti3,30) chst
        read (uniti3,30) fist
c If POLCOMS, change fist to top left box corner
        if (polcoms) then
         chst=chst-(0.5/ndch)
         fist=fist+(nrr-1)*(1./ndfi)+(0.5/ndfi)
        endif
  30    format(F8.4)
        read (uniti3,10) Dt2
c Convert Dt2 from hours to seconds
        Dt2=Dt2*3600

c If POLCOMS or OCCAM, create cspz array (from north to south)
        if (polcoms.or.occam) then
         ind=1
         do j=1,jm
          lattmp=fist-(j-1)*(1./ndfi)-(0.5/ndfi)
          do i=1,im
            cspz(ind)=cos(lattmp*Deg2Rad)
            ind=ind+1
          enddo
         enddo
        endif

        if (polcoms) then

c Read depths

         read(uniti2,*) ((hin(i,j),i=1,im),j=1,jm)

c Set depths to greater than minimum depth and reverse array

         icn=0
         nz=0
         do j=jm,1,-1
          do i=1,im
            icn=icn+1
            if (hin(i,j).gt.0.) then
              H(icn)=max(hmin,hin(i,j))
            else
              H(icn)=0.
            endif
          enddo
         enddo

        endif

        if (occam) then

c Read in list of grid box base depths

         do k=1,kb
          read (uniti,*) idum,hlev(k)
         enddo

c Read in KMT mask array and KMU array

         do j=1,jm
          do i=1,im
           read(uniti2,*) kmt(i,j)
           read(uniti4,*) kmu(i,j)
          enddo
         enddo

c Calculate depth array (and set index to start at top left)

         icn=0
         do j=jm,1,-1
          do i=1,im
           icn=icn + 1
           if (kmt(i,j).eq.0) then
            H(icn)=0
           else
            H(icn)=hlev(kmt(i,j))
           endif
          enddo
         enddo

c Calculate mask arrays (mz and mu)
         icn=0
         do j=jm,1,-1
          do i=1,im
           icn=icn + 1
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
          enddo
         enddo

        endif

c Calculate dch,dfi

        dch=1./ndch
        dfi=1./ndfi

        close (uniti)
        close (uniti2)
        if (occam) then
          close (uniti)
          close (uniti4)
        endif

c Deallocate arrays

        deallocate(hin)
        deallocate(kmt)
        deallocate(kmu)

 999    return

        End subroutine

