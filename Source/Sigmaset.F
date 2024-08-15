        Subroutine Sigmaset

c Subroutine to set up sigma levels consistent with POLCOMS

        use param

        implicit none

        Integer i, k
        Real*8 dsc,sds,dsu(kb),ds(kb),ztmp(kb),zztmp(kb)

c Evenly spaced sigma values

        ztmp(1) = -1.0d0
        ztmp(kb-1) = 0.0d0
        dsc = 1.0d0/(kb-2)
        do k=2,kb-2
          ztmp(k) = ztmp(k-1)+dsc
        enddo

        do k=1,kb-2
          ds(k) = ztmp(k+1)-ztmp(k)
        enddo

c Set surface level to ensure correct sum for ds

        sds = 0.0d0
        do k=1,kb-3
          sds = sds+ds(k)
        enddo
        ds(kb-2) = 1.0d0-sds

c Calculate zz(ki),value of sigma at u,v,b level given by ki.

        dsu(1)=ds(1)
        dsu(kb-1)=ds(kb-2)
        do k=2,kb-2
          dsu(k)=0.5d0*(ds(k)+ds(k-1))
        enddo
        zztmp(1)=-1.0d0+0.5d0*ds(1)
        do k=2,kb-2
          zztmp(k)=zztmp(k-1)+dsu(k)
        enddo
        zztmp(kb-1) = 0.0d0

c Now reverse z and zz so compatible with POM-based code
c Create dummy zzv array and convert all to 2D

        do i=1,Mx2
          do k=1,kb-1
            z(i,k)=ztmp(kb-k)
            zz(i,k)=zztmp(kb-k)
            zzv(i,k)=zztmp(kb-k)
          enddo
        enddo

        return

        End subroutine

