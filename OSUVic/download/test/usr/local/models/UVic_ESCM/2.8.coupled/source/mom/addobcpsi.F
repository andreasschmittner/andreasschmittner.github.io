      subroutine addobcpsi

!-----------------------------------------------------------------------
!     add open boundary values to the solution (active open boundary)
!-----------------------------------------------------------------------

#  include "param.h"
#  include "emode.h"
#  include "iounit.h"
#  include "obc_data.h"
#  include "tmngr.h"

#if defined obc_south || defined obc_north

!-----------------------------------------------------------------------
!       decide whether to read obc data or not
!-----------------------------------------------------------------------

        begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
        methodobc = 3
        call timeinterp (begtim, indxob1p, tobc1p, obc1pdpm, 12, .true.
     &,           methodobc, inextdobc1p, iprevdobc1p, wprevobc1p,
     &            readob1p, inextobc1p, iprevobc1p)

!       read in the next data record from disk when needed

        if (readob1p) then
          call getunit (ionew8, 'obcpsi1', opt_obcpsi1)
          read (ionew8, rec=inextdobc1p) obctnext, spdpmn, im, kk, jm
# if defined obc_south
     &,                 (psiwall_south(i,inextobc1p),i=1,imt)
# endif
# if defined obc_north
     &,                 (psiwall_north(i,inextobc1p),i=1,imt)
# endif
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obcpsi1 record =',inextdobc1p,' into buffer ='
     &,   inextobc1p
     &,   ' method #',methodobc,' at day =', begtim, stamp

        call relunit (ionew8)
        endif

      do i=1,imt
          tnext = c1-wprevobc1p
# if defined obc_south
          data           = tnext     *psiwall_south(i,inextobc1p)
     &                   + wprevobc1p*psiwall_south(i,iprevobc1p)
          psi(i,1,1)     = data
          psi(i,2,1)     = data
# endif
# if defined obc_north
          data           = tnext     *psiwall_north(i,inextobc1p)
     &                   + wprevobc1p*psiwall_north(i,iprevobc1p)
          psi(i,jmtm1,1) = data
          psi(i,jmt  ,1) = data
# endif
      enddo
#endif

#if defined obc_west || defined obc_east

!-----------------------------------------------------------------------
!       decide whether to read obc data or not
!-----------------------------------------------------------------------

        begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
        methodobc = 3
        call timeinterp (begtim, indxob2p, tobc2p, obc2pdpm, 12, .true.
     &,           methodobc, inextdobc2p, iprevdobc2p, wprevobc2p,
     &            readob2p, inextobc2p, iprevobc2p)

!       read in the next data record from disk when needed

        if (readob2p) then
          call getunit (ionew10, 'obcpsi2', opt_obcpsi2)
          read (ionew10, rec=inextdobc2p) obctnext, spdpmn, im, kk, jm
# if defined obc_west
     &,                 (psiwall_west(j,inextobc2p),j=1,jmt)
# endif
# if defined obc_east
     &,                 (psiwall_east(j,inextobc2p),j=1,jmt)
# endif
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obcpsi2 record =',inextdobc2p,' into buffer ='
     &,   inextobc2p
     &,   ' method #',methodobc,' at day =', begtim, stamp

        call relunit (ionew10)
        endif

      do j=1,jmt
          tnext = c1-wprevobc2p
# if defined obc_west
          data           = tnext     *psiwall_west(j,inextobc2p)
     &                   + wprevobc2p*psiwall_west(j,iprevobc2p)
          psi(1,j,1)     = data
          psi(2,j,1)     = data
# endif
# if defined obc_east
          data           = tnext     *psiwall_east(j,inextobc2p)
     &                   + wprevobc2p*psiwall_east(j,iprevobc2p)
          psi(imtm1,j,1) = data
          psi(imt,  j,1) = data
# endif
      enddo
#endif
#if defined agape_1x1 && defined obc_west && defined obc_east

!-----------------------------------------------------------------------
!     set southern land masses and perimeter points to southern values
!     of the meridional open boundaries
!-----------------------------------------------------------------------

      do jrow=1,jpsimax
        do i=1,imt
          if (map(i,jrow) .ne. 0) then
            psi(i,jrow,1) = psimax
            psi(i,jrow,2) = psimax
          endif
        enddo
      enddo
#endif
      return
      end
