#if defined held_larichev
      subroutine hlmixi

# include "size.h"
# include "accel.h"
# include "coord.h"
# include "hmixc.h"
# include "iounit.h"
# include "scalar.h"
# include "stdunits.h"

!=======================================================================
!     Initialization for ...
!     A Scaling Theory for Horizontally Homogeneous, Baroclinically
!     Unstable Flow on a Beta Plane (Submitted to JAS, March 1995).

!     output:
!       hl_depth = integration depth
!       hl_back  = minimum diffusion coefficient
!       hl_max   = maximum diffusion coefficient

!     based on code by: R. C. Pacanowski
!=======================================================================

      namelist /hlmix/ hl_depth, hl_back, hl_max

      write (stdout,'(/,20x,a,/)')
     & 'H L M I X    I N I T I A L I Z A T I O N'

!-----------------------------------------------------------------------
!     initialize variables (all mixing units are cm**2/sec.)
!-----------------------------------------------------------------------

      hl_depth = 500.0e2
      hl_back  = 1.e4
      hl_max   = 1.e9
      c0       = 0.0

!-----------------------------------------------------------------------
!     provide for namelist over-ride of above settings + documentation
!-----------------------------------------------------------------------

      call getunit (io, 'control.in', 'formatted sequential rewind')
      read (io,hlmix,end=100)
100   continue

!-----------------------------------------------------------------------
!     set no-flux condition on density difference across bottom level
!     initialize inverse richardson number = zero at all levels
!-----------------------------------------------------------------------

      do j=1,jmw
        do i=1,imt
          droz(i,km,j) = c0
        enddo
      enddo
      do j=1,jemw
        do k=1,km
          do i=1,imt
            rich_inv(i,k,j) = c0
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     check for problems
!-----------------------------------------------------------------------

# if !defined isopycmix && !defined isoneutralmix
        write (stdout,'(/,(1x,a))')
     & '==> Error:"isopycmix" must be enabled since "held_larichev"'
     &,'            predicts mixing coeffs used by "isopycmix"         '
        stop '=> hlmix'
# endif

!     write out namelist values

      write (stdout,hlmix)
      call relunit (io)
      return
      end

      subroutine hlmix (joff, js, je, is, ie)

!=======================================================================
!     Compute mixing coefficients based on...
!     A Scaling Theory for Horizontally Homogeneous, Baroclinically
!     Unstable Flow on a Beta Plane (Submitted to JAS, March 1995).

!     inputs:

!      joff   = offset between rows in the MW and latitude rows
!      js     = starting row for loading variables to calculate
!               coefficients. calculations start at jstrt=max(js-1,jsmw)
!      je     = ending row for loading variables to calculate
!               coefficients. calculations end at je-1
!      is     = starting index for calculating coefficients in the
!               longitude direction
!      ie     = ending index for calculating coefficients in the
!               longitude direction
!      grav   = gravity (cm/sec**2)
!      umask  = land/sea mask on "u" grid (land=0.0, sea=1.0)
!      tmask  = land/sea mask on "t" grid (land=0.0, sea=1.0)

!     outputs:

!      hl_coeff = held_larichev diffusion coefficient

!     based on code by: R. C. Pacanowski
!=======================================================================

      integer tlev
# include "size.h"
# include "coord.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
      dimension ro(imt,km,1:jmw)
      dimension sumz(imt,1:jemw)
!#define debug_hlmix
# if defined debug_hlmix
      common /debughl/ hlcoeff(imt,jmt)
      if (joff + js .eq. 1) then
        do j=1,jmt
          do i=1,imt
            hlcoeff(i,j) = 0.0
          enddo
        enddo
      endif
# endif

!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------

      p25   = 0.25
      c0    = 0.0
      p5    = 0.5
      epsln = 1.e-20
      fx    = -p25/grav
      istrt = max(2,is)
      iend  = min(imt-1,ie)
      tlev  = taum1
      kstart= 1
      kend  = max (1,min(indp(hl_depth, zt, km),km-1))

!-----------------------------------------------------------------------
!     set "ro" (density) at j=1 for 1st memory window otherwise ... move
!     variables from top two rows to bottom two rows to eliminate
!     redundant calculation
!-----------------------------------------------------------------------

      if (joff + js .eq. 1) then
        do k=1,km
          do i=istrt-1,iend+1
            ro(i,k,1) = c0
          enddo
        enddo
      else
        call movehlmix (istrt-1, iend+1)
      endif

!-----------------------------------------------------------------------
!     compute density difference across bottom of "t" cells at tau-1
!     for rows js through je in the MW. Set density difference = zero
!     across bottom and in land areas
!-----------------------------------------------------------------------

      do ks=1,2
        call statec (t(1,1,1,1,tlev), t(1,1,1,2,tlev), ro(1,1,jsmw)
     &,              max(js,jsmw), je, istrt-1, iend+1, ks)
        do j=js,je
          do k=ks,km-1,2
            do i=istrt-1,iend+1
              droz(i,k,j) = (ro(i,k,j) - ro(i,k+1,j))*tmask(i,k+1,j)
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     compute inverse richardson numbers on bottom of "u" cells
!-----------------------------------------------------------------------

      jsriu = max(js,jsmw)-1
      do j=jsriu,je-1
        do k=kstart,kend
          t1 = fx/dzw(k)
          do i=istrt,iend
            rich_inv(i,k,j) = t1*umask(i,k+1,j)*(
     &                    (u(i,k,j,1,tlev) - u(i,k+1,j,1,tlev))**2 +
     &                    (u(i,k,j,2,tlev) - u(i,k+1,j,2,tlev))**2) /
     &                    (droz(i,k,j+1) + droz(i+1,k,j+1) +
     &                     droz(i,k,j)   + droz(i+1,k,j) + epsln)
          enddo
        enddo
        call setbcx (rich_inv(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     integrate inverse Richardson number vertically over "u" cells
!     assume Ri at ocean surface is the same as at bottom of first cell
!-----------------------------------------------------------------------

      do j=jsriu,je-1
        k = kstart
        do i=istrt,iend
          hl_u(i,j) = rich_inv(i,k,j)*dzw(k-1)
          sumz(i,j) = dzw(k-1) + epsln
        enddo
        do k=kstart,kend
          do i=istrt,iend
            hl_u(i,j) = hl_u(i,j) + rich_inv(i,k,j)*dzw(k)
            sumz(i,j) = sumz(i,j) + dzw(k)*umask(i,k+1,j)
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     constrain integrated inverse Richardson numbers to be >= epsln
!-----------------------------------------------------------------------

      do j=jsriu,je-1
        do i=istrt,iend
          hl_u(i,j) = max(hl_u(i,j),epsln)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct Diffusion coefficient = 1/(beta**2*T**3) on "u" cells
!-----------------------------------------------------------------------

      do j=jsriu,je-1
        jrow = joff + j
        fsq = cori(jrow,1)**2
        beta_sq_r = 1.0/(2.0*omega*csu(jrow)/radius)**2
        do i=istrt,iend
          tm2 = fsq*hl_u(i,j)/sumz(i,j)
          hl_u(i,j) = beta_sq_r*tm2**1.5
        enddo
      enddo

!-----------------------------------------------------------------------
!     Limit Diffusion coefficient to be between hl_max and hl_back
!-----------------------------------------------------------------------

      do j=jsriu,je-1
        do i=istrt,iend
          hl_u(i,j) = max(min(hl_u(i,j),hl_max),hl_back)
        enddo
# if defined cyclic
        hl_u(1,j)   = hl_u(imt-1,j)
        hl_u(imt,j) = hl_u(2,j)
# else
        hl_u(1,j)   = c0
        hl_u(imt,j) = c0
# endif
      enddo

!-----------------------------------------------------------------------
!     construct Diffusion coefficient on eastern face of "T" cells
!-----------------------------------------------------------------------

      do j=jsmw,je-1
        do i=istrt,iend
          hl_e(i,j) = p5*(hl_u(i,j) + hl_u(i,j-1))
        enddo
# if defined cyclic
        hl_e(1,j)   = hl_e(imt-1,j)
        hl_e(imt,j) = hl_e(2,j)
# else
        hl_e(1,j)   = c0
        hl_e(imt,j) = c0
# endif
      enddo

!-----------------------------------------------------------------------
!     construct Diffusion coefficient centered in "T" cells
!-----------------------------------------------------------------------

      do j=jsmw,je-1
        do i=istrt,iend
          hl_b(i,j) = p25*(hl_u(i,j) + hl_u(i,j-1) + hl_u(i-1,j)
     &                   + hl_u(i-1,j-1))
        enddo
# if defined cyclic
        hl_b(1,j)   = hl_b(imt-1,j)
        hl_b(imt,j) = hl_b(2,j)
# else
        hl_b(1,j)   = c0
        hl_b(imt,j) = c0
# endif
      enddo

!-----------------------------------------------------------------------
!     construct Diffusion coefficient on northern face of "T" cells
!-----------------------------------------------------------------------

      do j=1,je-1
        do i=istrt,iend
          hl_n(i,j) = p5*(hl_u(i,j) + hl_u(i-1,j))
        enddo
# if defined cyclic
        hl_n(1,j)   = hl_n(imt-1,j)
        hl_n(imt,j) = hl_n(2,j)
# else
        hl_n(1,j)   = c0
        hl_n(imt,j) = c0
# endif
      enddo

#if defined debug_hlmix
      do j=jsriu,je-1
        do i=1,imt
          hlcoeff(i,jrow) = hl_u(i,j)
        enddo
        jrow = j + joff
        if (jrow .ge. 12 .and. jrow .le. 20 .and. prxzts .and. eots)then
          print *,' tp: jrow=',jrow
          call scope (rich_inv(1,1,j), imt, imt, km, 'rich_inv')
        endif
        if (jrow .eq. jmt-1 .and. prxzts .and. eots) then
          print *,'hlcoeff='
          call matrix (hlcoeff, imt, 1, imt, -1, -jmt, 0.0)
        endif
      enddo
#endif

      return
      end

      subroutine movehlmix (is, ie)

!=======================================================================
!     as the MW moves northward, move data from the last two rows
!     into the first two rows. (last 3 rows if using biharmonic option)

!     based on code by: R. C. Pacanowski
!=======================================================================

# include "param.h"
# include "hmixc.h"

      nrows = jmw - ncrows
      do move=1,nrows
        jfrom = jmw - (nrows - move)
        jto   = move

!-----------------------------------------------------------------------
!       move quantities with rows dimensioned (1:jmw)
!-----------------------------------------------------------------------

        do k=1,km
          do i=is,ie
            droz(i,k,jto)  = droz(i,k,jfrom)
          enddo
        enddo

!-----------------------------------------------------------------------
!       move quantities with rows dimensioned (1:jemw)
!-----------------------------------------------------------------------

        if (jfrom .le. jemw) then
          do k=1,km
            do i=is,ie
              rich_inv(i,k,jto)  = rich_inv(i,k,jfrom)
            enddo
          enddo
          do i=is,ie
            hl_u(i,jto) = hl_u(i,jfrom)
          enddo
        endif

!-----------------------------------------------------------------------
!       move quantities with rows dimensioned (jsmw:jmw)
!-----------------------------------------------------------------------

        if (jto .ge. jsmw) then

!         nothing to move

        endif

      enddo
      return
      end
#else
      subroutine hlmix
      return
      end
#endif
