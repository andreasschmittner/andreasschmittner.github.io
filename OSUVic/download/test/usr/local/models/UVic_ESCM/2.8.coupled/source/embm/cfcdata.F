      subroutine cfcdata
!=======================================================================
!     routine to read and interpolate one dimensional CFC forcing data

!     based on code by: M. Eby
!=======================================================================

      implicit none
#if defined uvic_cfc11 || defined uvic_cfc12

      integer iou, n, ln, ib(10), ic(10)

      logical exists

      real d(3,4), t(3), wt1, wt3

      real, allocatable :: data(:,:), time(:)

      save d, data, ln, t, time

# include "param.h"
# include "calendar.h"
# include "cembm.h"
# include "switch.h"
# include "tmngr.h"

      character(120) :: fname, name, new_file_name

      if (.not. allocated (time)) then
        name = "cfcccn.nc"
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "==> Warning: ", trim(fname), " does not exist."
          ln = 3
          allocate ( time(ln) )
          allocate ( data(ln,4) )
          time(:) = year0
          data(:,1) = cfc11ccnn
          data(:,2) = cfc11ccns
          data(:,3) = cfc12ccnn
          data(:,4) = cfc12ccns
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln,4) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          call getvara ('cfc11ccnn', iou, ln, ib, ic, data(1,1), c1, c0)
          call getvara ('cfc11ccns', iou, ln, ib, ic, data(1,2), c1, c0)
          call getvara ('cfc12ccnn', iou, ln, ib, ic, data(1,3), c1, c0)
          call getvara ('cfc12ccns', iou, ln, ib, ic, data(1,4), c1, c0)
          call closefile (iou)
        endif
        t(:) = time(1)
        d(:,1) = data(1,1)
        d(:,2) = data(1,2)
        d(:,3) = data(1,3)
        d(:,4) = data(1,4)
      endif

      t(2) = min(time(ln), max(time(1), year0 + relyr))

      if (t(2) .le. time(1)) then
        d(2,:) = data(1,:)
      elseif (t(2) .ge. time(ln)) then
        d(2,:) = data(ln,:)
      else
        if (t(2) .gt. t(3)) then
          do n=2,ln
            if (time(n-1) .le. t(2) .and. time(n) .ge. t(2)) then
              t(1) = time(n-1)
              d(1,:) = data(n-1,:)
              t(3) = time(n)
              d(3,:) = data(n,:)
            endif
          enddo
        endif
        wt1 = 1.
        if (t(3) .ne. t(1)) wt1 = (t(3)-t(2))/(t(3)-t(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        d(2,:) = d(1,:)*wt1 + d(3,:)*wt3
      endif

      cfc11ccnn = d(2,1)
      cfc11ccns = d(2,2)
      cfc12ccnn = d(2,3)
      cfc12ccns = d(2,4)
#endif

      return
      end
