      subroutine def_tsi
!=======================================================================
!     defines tsi files for UVic_ESCM

!     based on code by: M. Eby
!=======================================================================

      character(120) :: fname

#if defined uvic_embm
      call def_tsi_embm (fname)
#endif
#if defined uvic_mtlm
      call def_tsi_mtlm (fname)
#endif
#if defined uvic_ism
      call def_tsi_ism (fname)
#endif
#if defined uvic_mom
      call def_tsi_mom (fname)
#endif

      return
      end

#if defined uvic_embm
      subroutine def_tsi_embm (fname)
!=======================================================================
!     defines tsi files for the embm

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "iounit.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

      if (name .eq. ' ') then
        name = file_stamp ('tsi_embm',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
        call embm_tsi_def (name, timunit, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined uvic_mtlm
      subroutine def_tsi_mtlm (fname)
!=======================================================================
!     defines tsi files for the mtlm

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "iounit.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

      if (name .eq. ' ') then
        name = file_stamp ('tsi_mtlm',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
        call mtlm_tsi_def (name, timunit, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined uvic_ism
      subroutine def_tsi_ism (fname)
!=======================================================================
!     defines tsi files for the ism

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "iounit.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

      if (name .eq. ' ') then
        name = file_stamp ('tsi_ism',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
!        call ism_tsi_def (name, timunit, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined uvic_mom
      subroutine def_tsi_mom (fname)
!=======================================================================
!     defines tsi files for the mom

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

!      implicit none

# include "iounit.h"
# include "tmngr.h"

      character(120) :: fname, file_stamp, name, new_file_name
      integer iou, ntrec
      save name
      data name /' '/

      if (name .eq. ' ') then
        name = file_stamp ('tsi_mom',stamp,'.nc')
        name = new_file_name (name)
        call opennext (name, relyr, ntrec, iou)
        call closefile (iou)
        if (ntrec .gt. 1) then
          call opennew (name, iou)
          call closefile (iou)
        endif
        call mom_tsi_def (name, timunit, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif
