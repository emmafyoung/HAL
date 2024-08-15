      Subroutine Check(status,endflag,operation)

      use netcdf
      implicit none

      integer, intent(in) :: status
      integer, intent(inout) :: endflag
      character(len=*), intent(in) :: operation

      if (status == NF90_NOERR) return
      print *, "Error encountered during ", operation
      print *, nf90_strerror(status)
      endflag=-1

      end subroutine Check
