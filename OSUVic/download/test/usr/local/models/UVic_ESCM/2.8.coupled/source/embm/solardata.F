       subroutine solardata
!=======================================================================
!     routine to read and interpolate solar forcing data

!     based on code by: M. Eby
!=======================================================================

      implicit none
# if defined uvic_embm_solardata || defined uvic_embm_solardata_transient

      integer iou, n, ln, ib(10), ic(10)

      logical exists

      real d(3), t(3)
      real c1e3, wt1, wt3

      real, allocatable :: data(:), time(:)

      save d, data, ln, t, time

# include "param.h"
# include "cembm.h"
# include "tmngr.h"

      character(120) :: fname, name, new_file_name

      c1e3 = 1.e3
      if (.not. allocated (time)) then
        fname = "solar.nc"
        name = new_file_name (fname)
        inquire (file=trim(name), exist=exists)
        if (.not. exists) then
          print*, "==> Warning: ", trim(name), " does not exist."
          ln = 3
          allocate ( time(ln) )
          allocate ( data(ln) )
          time(:) = year0
          data(:) = solarconst
        else
          call openfile (name, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          call getvara ('solar_constant', iou, ln, ib, ic, data
     &,     c1e3, c0)
          call closefile (iou)
        endif
        t(:) = time(1)
        d(:) = data(1)
      endif

#  if defined uvic_embm_solardata_transient
      t(2) = min(time(ln), max(time(1), year0 + relyr))
#  else
      t(2) = min(time(ln), max(time(1), solar_yr))
#  endif

      if (t(2) .le. time(1)) then
        d(2) = data(1)
      elseif (t(2) .ge. time(ln)) then
        d(2) = data(ln)
      else
        if (t(2) .gt. t(3)) then
          do n=2,ln
            if (time(n-1) .le. t(2) .and. time(n) .ge. t(2)) then
              t(1) = time(n-1)
              d(1) = data(n-1)
              t(3) = time(n)
              d(3) = data(n)
            endif
          enddo
        endif
        wt1 = 1.
        if (t(3) .ne. t(1)) wt1 = (t(3)-t(2))/(t(3)-t(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        d(2) = d(1)*wt1 + d(3)*wt3
      endif

      solarconst = d(2)
# endif

      return
      end
