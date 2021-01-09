      subroutine def_tavg
!=======================================================================
!     defines tavg files for UVic_ESCM

!     based on code by: M. Eby
!=======================================================================

      character(120) :: fname

#if defined uvic_embm
      call def_tavg_embm (fname)
#endif
#if defined uvic_mtlm
      call def_tavg_mtlm (fname)
#endif
#if defined uvic_ism
      call def_tavg_ism (fname)
#endif
#if defined uvic_mom
      call def_tavg_mom (fname)
#endif

      return
      end

#if defined uvic_embm
      subroutine def_tavg_embm (fname)
!=======================================================================
!     defines tavg files for the embm

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "size.h"
# include "coord.h"
# include "atm.h"
# include "ice.h"
# include "iounit.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

# if defined uvic_multiple_files
      name = ' '
# endif
      if (name .eq. ' ') then
        name = file_stamp ('tavg_embm',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
        call embm_snap_def (name, imt, jmt, nat, ncat, xt, yt, timunit
     &,                     expnam, runstamp, mapat)
      endif
      fname = name

      return
      end
#endif

#if defined uvic_mtlm
      subroutine def_tavg_mtlm (fname)
!=======================================================================
!     defines tavg files for the mtlm

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "size.h"
# include "coord.h"
# include "iounit.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

# if defined uvic_multiple_files
      name = ' '
# endif
      if (name .eq. ' ') then
        name = file_stamp ('tavg_mtlm',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
        call mtlm_snap_def (name, imt, jmt, NPFT, NTYPE, xt, yt
     &,                     timunit, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined uvic_ism
      subroutine def_tavg_ism (fname)
!=======================================================================
!     defines tavg files for the ism

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "size.h"
# include "coord.h"
# include "iounit.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

# if defined uvic_multiple_files
      name = ' '
# endif
      if (name .eq. ' ') then
        name = file_stamp ('tavg_ism',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
!        call ism_snap_def (name, imt, jmt, xt, yt
!     &,                    timunit, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined uvic_mom
      subroutine def_tavg_mom (fname)
!=======================================================================
!     defines tavg files for the mom

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

!      implicit none

# include "param.h"
# include "coord.h"
# include "iounit.h"
# include "mw.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

# if defined uvic_multiple_files
      name = ' '
# endif
      if (name .eq. ' ') then
        name = file_stamp ('tavg_mom',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
        call mom_snap_def (name, imt, jmt, km, nt, kpzd, xt, yt, timunit
     &,                    expnam, runstamp, mapt)
      endif
      fname = name

      return
      end
#endif
