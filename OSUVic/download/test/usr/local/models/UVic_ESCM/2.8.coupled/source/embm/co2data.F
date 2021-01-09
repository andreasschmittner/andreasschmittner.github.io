       subroutine co2data
!=======================================================================
!     routine to read and interpolate one dimensional forcing data
!     data is in ppmv CO2 or if uvic_carbon_coupled, Gt carbon/year

!     based on code by: M. Eby
!=======================================================================

      implicit none
#if defined uvic_embm_co2data || defined uvic_embm_co2data_transient

      integer iou, n, ln, ib(10), ic(10)

      logical exists

      real d(3), t(3), wt1, wt3

      real, allocatable :: data(:), time(:)

      save d, data, ln, t, time

#include "param.h"
#include "calendar.h"
#include "cembm.h"
#include "switch.h"
#include "tmngr.h"
#if defined uvic_carbon_co2_2d
# include "atm.h"
      real tmp, dmsk(imt,jmt)
#endif

      character(120) :: fname, name, new_file_name

      if (.not. allocated (time)) then
# if defined uvic_carbon_coupled
        name = "co2emit.nc"
# else
        name = "co2ccn.nc"
# endif
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "==> Warning: ", trim(fname), " does not exist."
          ln = 3
          allocate ( time(ln) )
          allocate ( data(ln) )
          time(:) = year0
# if defined uvic_carbon_coupled
          data(:) = 0.
# else
          data(:) = co2ccn
# endif
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
# if defined uvic_carbon_coupled
          call getvara ('co2emit', iou, ln, ib, ic, data, c1, c0)
# else
          call getvara ('co2ccn', iou, ln, ib, ic, data, c1, c0)
# endif
          call closefile (iou)
        endif
        t(:) = time(1)
        d(:) = data(1)
      endif

# if defined uvic_embm_co2data_transient
      t(2) = min(time(ln), max(time(1), year0 + relyr))
# else
      t(2) = min(time(ln), max(time(1), co2_yr))
# endif

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

# if defined uvic_carbon_coupled
!     convert flux from Gt year-1 to g cm-2 s-1 (surface area = 5.1e18 cm2)
      d(2) = d(2)*1.e15/(5.1e18*yrlen*daylen)
#  if defined uvic_carbon_co2_2d
!     convert flux from g cm-2 s-1 to ppmv s-1
!     4.138e-7 => 12e-6 g (umol carbon)-1 / 29 g (mol air)-1
      flux(:,:,ico2) = flux(:,:,ico2) - d(2)/(4.138e-7*rhoatm*shc)
#  else
!     convert flux from g cm-2 s-1 to umol cm-2 s-1
      co2flx = co2flx - d(2)/12.e-6
#  endif
# else
#  if defined uvic_carbon_co2_2d
      dmsk(:,:) = 1.
      call areaavg (at(1,1,1,ico2), dmsk, tmp)
      at(:,:,1,ico2) = at(:,:,1,ico2) + d(2) - tmp
      dmsk(:,:) = 1.
      call areaavg (at(1,1,2,ico2), dmsk, tmp)
      at(:,:,2,ico2) = at(:,:,2,ico2) + d(2) - tmp
#  else
      co2ccn = d(2)
#  endif
# endif
#endif

      return
      end
