      subroutine def_rest (last)
!=======================================================================
!     defines rest files for UVic_ESCM

!     input:
!       last = last rest flag (1 = last)

!     based on code by: M. Eby
!=======================================================================

      character(120) :: fname
      integer last

#if defined uvic_embm
      call def_rest_embm (last, fname)
#endif
#if defined uvic_mtlm
      call def_rest_mtlm (last, fname)
#endif
#if defined uvic_ism
      call def_rest_ism (last, fname)
#endif
#if defined uvic_mom
      call def_rest_mom (last, fname)
#endif

      return
      end

#if defined uvic_embm
      subroutine def_rest_embm (last, fname)
!=======================================================================
!     defines rest file for the embm

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "tmngr.h"

      character(120) :: fname, name, file_stamp, new_file_name
      integer last

      if (last .eq. 1) then
        name = 'restart_embm.nc'
      else
        name = file_stamp ('rest_embm',stamp,'.nc')
      endif
      name = new_file_name (name)
      call embm_rest_def (name)
      fname = name

      return
      end
#endif

#if defined uvic_mtlm
      subroutine def_rest_mtlm (last, fname)
!=======================================================================
!     defines rest file for the mtlm

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "tmngr.h"

      character(120) :: fname, name, file_stamp, new_file_name
      integer last

      if (last .eq. 1) then
        name = 'restart_mtlm.nc'
      else
        name = file_stamp ('rest_mtlm',stamp,'.nc')
      endif
      name = new_file_name (name)
      call mtlm_rest_def (name)
      fname = name

      return
      end
#endif

#if defined uvic_ism
      subroutine def_rest_ism (last, fname)
!=======================================================================
!     defines rest file for the ism

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "tmngr.h"

      character(120) :: fname, name, file_stamp, new_file_name
      integer last

      if (last .eq. 1) then
        name = 'restart_ism.nc'
      else
        name = file_stamp ('rest_ism',stamp,'.nc')
      endif
      name = new_file_name (name)
      call ism_rest_def (name)
      fname = name

      return
      end
#endif

#if defined uvic_mom
      subroutine def_rest_mom (last, fname)
!=======================================================================
!     defines rest file for the mom

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "tmngr.h"

      character(120) :: fname, name, file_stamp, new_file_name
      integer last

      if (last .eq. 1) then
        name = 'restart_mom.nc'
      else
        name = file_stamp ('rest_mom',stamp,'.nc')
      endif
      name = new_file_name (name)
      call mom_rest_def (name)
      fname = name

      return
      end
#endif
