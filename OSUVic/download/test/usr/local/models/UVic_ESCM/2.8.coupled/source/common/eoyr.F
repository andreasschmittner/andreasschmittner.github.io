
      subroutine eoyr (is, ie, js, je)
!=======================================================================
!     things to do at the end of a year

!     input:
!       is, ie, js, je = starting and ending indicies for i and j

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "param.h"
# include "calendar.h"
# include "csbc.h"
# include "levind.h"
# include "atm.h"
# include "coord.h"
# include "grdvar.h"
# include "cembm.h"
# include "iounit.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
# if defined uvic_ice
#  include "ice.h"
# endif
# if defined uvic_mtlm
#  include "mtlm.h"
# endif
# if defined uvic_embm_orbit_transient
#  include "insolation.h"
# endif

      character(120) :: fname, new_file_name
      character(3) :: a3

      integer i, ie, is, iou, j, je, js, L, ib(10), ic(10), n

      logical exists, inqvardef

      real wt, f, dmsk(imt,jmt)

# if defined uvic_embm_orbit_transient

      orbit_yr = year0 + relyr
      call orbit (orbit_yr, eccen, obliq, mvelp, lambm0)
# endif
# if defined uvic_embm_icedata_transient

      call icedata

      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (tmsk(i,j), tmsk(i+1,j), tmsk(i,j+1)
     &,                    tmsk(i+1,j+1))
        enddo
      enddo
      call embmbc (umsk)
!     remove isolated bays
      do j=2,jmtm1
        do i=2,imtm1
          tmsk(i,j) = max (umsk(i,j), umsk(i-1,j), umsk(i,j-1)
     &,                    umsk(i-1,j-1))
        enddo
      enddo
      call embmbc (tmsk)
      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (tmsk(i,j), tmsk(i+1,j), tmsk(i,j+1)
     &,                    tmsk(i+1,j+1))
        enddo
      enddo
      call embmbc (umsk)

#  if uvic_ice && !defined uvic_ice_cpts
      do j=1,jmt
        do i=1,imt
          if (tmsk(i,j) .ge. 0.5) then
            if (hice(i,j,1) .le. 0.) aice(i,j,1) = 0.
            if (hice(i,j,2) .le. 0.) aice(i,j,2) = 0.
          endif
        enddo
      enddo
      call embmbc (aice(1,1,1))
      call embmbc (aice(1,1,2))

#  endif
!     recalculate rivers to be consistent with land masks
      call  rivinit
#  if defined uvic_mtlm

!     set flag to recalculate land index
      NEW_LAND = .true.
#  endif
# endif
# if defined uvic_embm_running_average || defined uvic_embm_awind

!     make new annual and running average
      call embmbc (atbar)
      if (totaltime .ge. yrlen*daylen - dtatm*0.5) then
        wt = 0.05
        do j=js,je
          do i=is,ie
            rtbar(i,j) = (1.-wt)*rtbar(i,j) + wt*atbar(i,j)/totaltime
            atbar(i,j) = 0.0
          enddo
        enddo
      endif
      totaltime = 0.0
#  if defined uvic_embm_awind

      call calc_awind (is, ie, js, je)
#  endif
# endif
#endif

      return
      end
