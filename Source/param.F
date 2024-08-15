      module param
c
c     Common declarations for the IBM model
c
#include "comblk_ibm.h"
c
      real, allocatable, dimension(:,:) :: nemolev,U,V,W,Tem,Sal
      real, allocatable, dimension(:) :: Ice,Uice,Vice
      real, allocatable, dimension(:,:) :: nemowlev
      real, allocatable, dimension(:,:) :: PKh,PAh
      real, allocatable, dimension(:,:) :: U1,U2,V1,V2,W1,W2
      real, allocatable, dimension(:,:) :: T1,T2,Sa1,Sa2
      real, allocatable, dimension(:) :: Ice1,Ice2,Uice1,Uice2
      real, allocatable, dimension(:) :: Vice1,Vice2
      real, allocatable, dimension(:,:) :: PAh1,PAh2,PKh1,PKh2
      real*8, allocatable, dimension(:,:) :: z,zz,zzv
      real, allocatable, dimension(:,:,:) :: Drift1,Drift2

      end

