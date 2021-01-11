      subroutine glsbc
!-----------------------------------------------------------------------
!     Get the boundary conditions for the land model.

!     based on code by: K. Meissner, M. Eby
!-----------------------------------------------------------------------
#if defined uvic_mtlm

       implicit none

# include "param.h"
# include "coord.h"
# include "cembm.h"
# include "csbc.h"
# include "calendar.h"
# include "tmngr.h"
# include "switch.h"
# include "levind.h"
# include "mtlm.h"
# include "insolation.h"
# include "atm.h"
# if defined uvic_embm_cropdata_transient
#  include "veg.h"
# endif
# if defined uvic_ice && defined uvic_embm_icedata
#  include "ice.h"
# endif

      integer i, j, L, n

      real calday, fc, fd, fe, ff, cosz(POINTS), dt, t, pi, degrad

      fc = 1.0/atatm
      fd = 1.e-2/atatm
      fe = 1.e-3/atatm
      ff = 10./atatm
      atlnd = 0.
      pi = 4.*atan(1.)
      degrad = pi/180.

      if (NEW_LAND) then
!----------------------------------------------------------------------
!       Create the VEG_INDEX array of land points with each type
!----------------------------------------------------------------------
        L = 0
        land_map(:,:) = 0
        LAND_PTS = 0
        LAND_INDEX(:) = 0
        do j=2,jmt-1
          do i=2,imt-1
            if (kmt(i,j) .le. klmax) then
              L = L + 1
# if defined uvic_ice && defined uvic_embm_icedata
              if (aicel(i,j,2) .lt. 0.5 .and. tmsk(i,j) .lt. 0.5) then
                land_map(i,j) = L
                LAND_PTS = LAND_PTS + 1
                LAND_INDEX(LAND_PTS) = L
              endif
# else
              if (tmsk(i,j) .lt. 0.5) then
                land_map(i,j) = L
                LAND_PTS = LAND_PTS + 1
                LAND_INDEX(LAND_PTS) = L
              endif
# endif
            endif
          enddo
        enddo
        do N=1,NPFT
          VEG_PTS(N) = 0
          do J=1,LAND_PTS
            L = LAND_INDEX(J)
            if (FRAC(L,N) .gt. FRAC_MIN + epsln) then
              VEG_PTS(N) = VEG_PTS(N) + 1
              VEG_INDEX(VEG_PTS(N),N) = L
            endif
          enddo
        enddo
        NEW_LAND = .false.
      endif

!----------------------------------------------------------------------
!     Calculate the diurnal cycle in the SW radiation
!----------------------------------------------------------------------
      calday = dayoyr*365.25/yrlen
      call decl (calday, eccen, obliq, mvelp, lambm0, sindec, eccf)
      dt = SEC_DAY/STEP_DAY
      do n=1,STEP_DAY
        t = real(n-1)*dt
        call zenith (POINTS, t, dt, SEC_DAY, LAT, LONG, sindec, cosz)
        SUN(:,n) = solarconst*eccf*cosz(:)
      enddo
      do j=2,jmt-1
        do i=2,imt-1
          L = land_map(i,j)
          if (L .ne. 0) then
            t = 0.
            do n=1,STEP_DAY
              t = t + SUN(L,n)
            enddo
!           make sure the daily insolations agree
            SUN(L,:) = SUN(L,:)*solins(i,j)*STEP_DAY/(t + epsln)
!           calculate downward shortwave at the surface
            SUN(L,:) = SUN(L,:)*1.e-3*sbc(i,j,iaca)*pass*sbc(i,j,isca)
          endif
        enddo
      enddo

!----------------------------------------------------------------------
!     Calculate the time of maximum temperature. Assume at local noon.
!----------------------------------------------------------------------
      do L=1,POINTS
        TIME_MAX(L) = SEC_DAY*(0.5 + LONG(L)/360.)
        if (TIME_MAX(L) .lt. 0) then
          TIME_MAX(L) = TIME_MAX(L)+SEC_DAY*int(1.+TIME_MAX(L)/SEC_DAY)
        elseif (TIME_MAX(L) .gt. SEC_DAY) then
          TIME_MAX(L) = TIME_MAX(L)-SEC_DAY*int(TIME_MAX(L)/SEC_DAY)
        endif
      enddo

!----------------------------------------------------------------------
!     Set other boundary conditions and zero accumulators
!----------------------------------------------------------------------
      do j=2,jmtm1
        do i=2,imtm1
          L = land_map(i,j)
          if (L .ne. 0) then
!           set boundary conditions for land
            RAIN(L) = ff*sbc(i,j,ipr)
            SNOW(L) = ff*sbc(i,j,ips)
            SW_C(L) = fe*sbc(i,j,iswr)
            T_C(L) = fc*sbc(i,j,iat) + 273.15
            WIND(L) = fd*sbc(i,j,iaws)
            RH_C(L) = fc*sbc(i,j,irh)
            DTEMP_DAY(L) = sbc(i,j,idtr)
            TSTAR(L,:) = T_C(L)
            TSOIL(L) = T_C(L)
# if defined uvic_carbon_co2_2d
            CO2(L) = 1.0E-6*at(i,j,2,ico2)*EPCO2
# else
            CO2(L) = 1.0E-6*co2ccn*EPCO2
# endif
# if defined uvic_embm_cropdata_transient
            FRACA(L) = crops(i,j,2)
# endif
!           zero atmosphere accumulation variables
            sbc(i,j,iro) = 0.
            sbc(i,j,ievap) = 0.
            sbc(i,j,ilwr) = 0.
            sbc(i,j,isens) = 0.
            sbc(i,j,isca) = 0.
# if defined uvic_carbon
            sbc(i,j,inpp) = 0.
            sbc(i,j,isr) = 0.
# endif
          endif
        enddo
      enddo
# if defined time_averages

!-----------------------------------------------------------------------
!     zero time averages if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. timavgperts) call ta_mtlm_snap (0)
# endif
# if defined time_step_monitor

!-----------------------------------------------------------------------
!     zero time step integrals if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. tsiperts) call ta_mtlm_tsi (0)
# endif
#endif

      return
      end
