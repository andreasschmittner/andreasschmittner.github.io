      subroutine embmbc (data)

!=======================================================================
!     set boundary values

!     input:
!       data = array to be set

!     based on code by: A. Fanning and M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "param.h"

      integer i, j

      real data(imt,jmt)

!     set cyclic east west

      do j=1,jmt
        data(1,j) = data(imtm1,j)
        data(imt,j) = data(2,j)
      enddo

!     set closed North South

      do i=1,imt
        data(i,1) = data(i,2)
        data(i,jmt) = data(i,jmtm1)
      enddo
#endif

      return
      end
