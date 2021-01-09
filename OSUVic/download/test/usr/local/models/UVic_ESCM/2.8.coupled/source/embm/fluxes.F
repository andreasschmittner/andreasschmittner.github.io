      subroutine fluxes (is, ie, js, je)

!=======================================================================
!     calculate energy and moisture fluxes

!     Note: evaporation and precipitation are in g cm-2 s-1
!           and humidities are in g g-1

!     for Thompson and Warren outgoing radiation (see: Thompson S.J.,
!     and S.G. Warren 'parameterization of outgoing ...'J. Atmos. Sci.,
!     39, 2667-2680, 1982.

!     based on code by: A. Fanning, D. Matthews and M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "param.h"
# include "cembm.h"
# include "atm.h"
# include "csbc.h"
# if defined uvic_ice
#  include "ice.h"
# endif
# include "veg.h"

      integer i, ie, iem1, imax, is, isp1, iter
      integer j, je, jem1, jmax, js, jsp1, maxit

      real b00, b10, b20, b01, b11, b21, b02, b12, b22, b03, b13, b23
      real delta, df, dt, dultnt, dulwr, dusens, emax, f, fb, ff, fg
      real fh, fl, fm, qair, qlnd, rhrh, sr, scrit, ssh, tair, tlnd
      real tlold, tol, tol2, ultnt, ulwr, usens, wspd

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

!-----------------------------------------------------------------------
!     set appropriate constants
!-----------------------------------------------------------------------
      fb = 0.94*rhoatm*cpatm
      maxit = 10
      tol = 0.01
      emax = 0.0
      imax = 0
      jmax = 0
      scrit = 0.75*soilmax
      ff = rhoatm*vlocn

!     Thomson and Warren constants

      b00 = 2.375e2
      b10 = -3.47968e1
      b20 = 1.02790e1
      b01 = 2.60065
      b11 = -1.62064
      b21 = 6.34856e-1
      b02 = 4.40272e-3
      b12 = -2.26092e-2
      b22 = 1.12265e-2
      b03 = -2.05237e-5
      b13 = -9.67e-5
      b23 = 5.62925e-5

      tol2 = tol*2.0

      do j=jsp1,jem1
        do i=isp1,iem1

!-----------------------------------------------------------------------
!         set the incoming short wave
!-----------------------------------------------------------------------
          dnswr(i,j) = solins(i,j)*sbc(i,j,iaca)*pass*sbc(i,j,isca)

!-----------------------------------------------------------------------
!         set wind speed and elevated air temperature
!-----------------------------------------------------------------------
          wspd = sbc(i,j,iws)
          tair = at(i,j,2,isat) - elev(i,j)*rlapse*rfactor
# if defined uvic_embm_icedata
     &         - hicel(i,j,2)*rlapse*rfactor
# endif

!-----------------------------------------------------------------------
!         calculate outgoing longwave radiation
!-----------------------------------------------------------------------
          rhrh = rh(i,j)*rh(i,j)
          outlwr(i,j) = 1.0e3*(b00 + b10*rh(i,j) + b20*rhrh
     &                + (b01 + b11*rh(i,j) + b21*rhrh)*tair
     &                + (b02 + b12*rh(i,j) + b22*rhrh)*tair**2
     &                + (b03 + b13*rh(i,j) + b23*rhrh)*tair**3)
# if defined uvic_carbon && defined uvic_carbon_co2_2d
     &                - co2for*alog(at(i,j,2,ico2)/350.0)
# else
     &                - anthro
# endif

          tair = at(i,j,2,isat) - elev(i,j)*rlapse
# if defined uvic_embm_icedata
     &         - hicel(i,j,2)*rlapse
# endif

          if (tmsk(i,j) .ge. 0.5) then

!-----------------------------------------------------------------------
!           calculations only for ocean points
!-----------------------------------------------------------------------
            dt = sbc(i,j,isst) - at(i,j,2,isat)
            fg = dalt_o*wspd

!-----------------------------------------------------------------------
!           calculate evaporation or sublimation (ensure it is positive)
!-----------------------------------------------------------------------
            ssh = cssh*exp(17.67*sbc(i,j,isst)/(sbc(i,j,isst) + 243.5))
            evap(i,j) = max(c0, rhoatm*fg*(ssh - at(i,j,2,ishum)))
            upltnt(i,j) = vlocn*evap(i,j)

!-----------------------------------------------------------------------
!           calculate upward sensible heat flux
!-----------------------------------------------------------------------
            upsens(i,j) = fb*fg*(dt)

!-----------------------------------------------------------------------
!           calculate upward longwave re-radiation
!-----------------------------------------------------------------------
            uplwr(i,j) = esocn*(sbc(i,j,isst) + 273.15)**4
     &                 - esatm*(at(i,j,2,isat) + 273.15)**4
# if defined uvic_mtlm

          elseif (land_map(i,j) .ne. 0) then

!----------------------------------------------------------------------
!           set fluxes over land from the land model
!---------------------------------------------------------------------
            upltnt(i,j) = 0.0
            evap(i,j) = sbc(i,j,ievap)
            upsens(i,j) = sbc(i,j,isens)
            uplwr(i,j) = sbc(i,j,ilwr)
# endif

          else

!-----------------------------------------------------------------------
!            calculations only for land points

!           find land temperature by balancing the surface heat budget
!             dwsr = ultnt + usens + ulwr
!           using Newton's method:
!             t(i+1) = t(i) - f(t(i))/df(t(i))
!           where:
!             f(t(i)) = dwsr - ultnt - usens - ulwr
!             -df(t(i)) = dultnt - dusens - dulwr
!-----------------------------------------------------------------------
            tlnd = surf(i,j)
            tlold = tlnd
            fm = esatm*(tair + 273.15)**4
            fg = rhoatm
!           calculate stomatal resistance
# if defined uvic_embm_cropdata
            sr = (1.-crops(i,j,2))*veg_rs(iveg(i,j))
     &         + crops(i,j,2)*veg_rs(icrops)
            dalt_v = veg_dalt(i,j)
# else
            sr = veg_rs(iveg(i,j))
            dalt_v = veg_dalt(iveg(i,j))
# endif
!           add in aerodynamic resistance
            sr = sr + 1.0/(dalt_v*wspd + epsln)
!           set beta parameter for calculating actual evaporation
            fh = min(max(c0+epsln, (soilm(i,j,lf)/soilmax)**(0.25)),c1)
!           set coefficients for latent heat (fl) and evaporation (fg)
            fl = fh*ff/(sr)
            fg = fh*fg/(sr)
            dusens = fb*dalt_v*wspd

!-----------------------------------------------------------------------
!           start loop for all land grid points
!-----------------------------------------------------------------------
            qair = rh(i,j)*cssh*exp(17.67*tair/(tair + 243.5))
            iter = 0
            delta = tol2
            do while (abs(delta) .gt. tol .and. iter .le. maxit)
              iter = iter + 1
              qlnd = cssh*exp(17.67*tlnd/(tlnd + 243.5))
              if (qlnd .gt. qair) then
                ultnt = fl*(qlnd - qair)
                dultnt = fl*qlnd*17.67*243.5/(tlnd + 243.5)**2
              else
                ultnt = 0.0
                dultnt = 0.0
              endif
              usens = dusens*(tlnd - tair)
              ulwr = eslnd*(tlnd + 273.15)**4 - fm
              dulwr = 4.0*eslnd*(tlnd + 273.15)**3
              f = dnswr(i,j) - ultnt - usens - ulwr
              df = dultnt + dusens + dulwr
              delta = f/df
              tlnd = tlnd + delta
            enddo
            if (iter .gt. maxit) then
!             if not converged, set to last converged temperature
              if (abs(delta) .gt. emax) then
                emax = abs(delta)
                imax = i
                jmax = j
                tlnd = tlold
              endif
            endif

!-----------------------------------------------------------------------
!           calculate fluxes on land
!-----------------------------------------------------------------------
            surf(i,j) = tlnd
            qlnd = cssh*exp(17.67*tlnd/(tlnd + 243.5))
            evap(i,j) = max(c0, fg*(qlnd - qair))
            evap(i,j) = max(c0, min(soilm(i,j,lf)/dts, evap(i,j)))
# if defined uvic_ice
            flux(i,j,ishum) = evap(i,j)*(1.0 - aice(i,j,2))
# else
            flux(i,j,ishum) = evap(i,j)
# endif
            upltnt(i,j) = vlocn*evap(i,j)
            upsens(i,j) = dusens*(tlnd - tair)
            uplwr(i,j) = eslnd*(tlnd + 273.15)**4 - fm

!           ensure fluxes are balanced since land can't absorb error
            upsens(i,j) = dnswr(i,j) - upltnt(i,j) - uplwr(i,j)

          endif

        enddo
      enddo

      if (emax .gt. 0.0) write (stdout,*)
     &  '==> Warning: land surface temperature not converging: '
     &, 'emax, i, j, soilm:', emax, imax, jmax, soilm(imax,jmax,2)

      return
      end

      subroutine precipitate (is, ie, js, je)

!=======================================================================
!     calculate precipitation explicitly and update humidity

!     based on code by: A. Fanning, D. Matthews and M. Eby
!=======================================================================

      implicit none

# include "param.h"
# include "scalar.h"
# include "cembm.h"
# include "atm.h"
# include "switch.h"
# if defined uvic_ice
#  include "ice.h"
# endif
# if defined uvic_mtlm
#  include "csbc.h"
#  include "mtlm.h"
# endif

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1, k, n, negq

      real fb, fc, qmax, rate, tair, soiltemp, pson, psot, ssh
      real tmp
# if defined uvic_mtlm
      real hs(imt,jmt)
# endif

      data negq /0/
      save negq

      if (eoyear) negq = 0

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

!-----------------------------------------------------------------------
!     set appropriate constants
!-----------------------------------------------------------------------
      fb = rhoatm*shq/dts
      fc = dts/rhosno
!     maximum relative humidity after rain
# if defined uvic_mtlm
      call unloadland (POINTS, LYING_SNOW, imt, jmt, land_map, hs)
!     convert from kg/m2 to cm
      hs(:,:) = hs(:,:)*0.1/rhosno
# endif

      precip(:,:) = 0.
      rh(:,:) = 0.
      k = 0
      do j=jsp1,jem1
        do i=isp1,iem1

!-----------------------------------------------------------------------
!         add rain due to orography
!-----------------------------------------------------------------------
          if (tmsk(i,j) .lt. 0.5) then
            tair = at(i,j,2,isat) - elev(i,j)*rlapse
# if defined uvic_embm_icedata
     &           - hicel(i,j,2)*rlapse
# endif
            ssh = cssh*exp(17.67*tair/(tair + 243.5))
            qmax = rhmax*ssh
            if (at(i,j,2,ishum) .gt. qmax) then
              tmp = fb*(at(i,j,2,ishum) - qmax)*rfactor
              precip(i,j) = precip(i,j) + tmp
              at(i,j,2,ishum) = at(i,j,2,ishum) - tmp/fb
              rh(i,j) = rhmax
            else
              rh(i,j) = max(c0, min(c1, at(i,j,2,ishum)/(ssh + epsln)))
            endif
          endif

!-----------------------------------------------------------------------
!         check if specific humidity is greater than rhmax of saturation
!-----------------------------------------------------------------------
          tair = at(i,j,2,isat)
          ssh = cssh*exp(17.67*tair/(tair + 243.5))
          qmax = rhmax*ssh
          if (at(i,j,2,ishum) .gt. qmax) then
            tmp = fb*(at(i,j,2,ishum) - qmax)
            precip(i,j) = precip(i,j) + tmp
            at(i,j,2,ishum) = at(i,j,2,ishum) - tmp/fb
            rh(i,j) = rhmax
          else
            rh(i,j) = max(rh(i,j), at(i,j,2,ishum)/(ssh + epsln))
            rh(i,j) = max(c0, min(c1, rh(i,j)))
          endif

!-----------------------------------------------------------------------
!         calculate snowfall (hsno at tau was set in the ice model)
!-----------------------------------------------------------------------
!         tair may be adjusted by a snowfall offset temperature tsno
          tair = at(i,j,2,isat) - tsno - elev(i,j)*rlapse
# if defined uvic_embm_icedata
     &         - hicel(i,j,2)*rlapse
# endif
          psno(i,j) = 0.0
# if defined uvic_ice_cpts
#  if defined uvic_mtlm
          hs(i,j) = hs(i,j) + hseff(i,j,idx,1)
          if (tair .le. c0 .and. hs(i,j) .lt. hsno_max)
     &      psno(i,j) = min(precip(i,j), hsno_max - hs(i,j))
#  else
          if (tair .le. c0 .and. hseff(i,j,idx,1) .lt. hsno_max)
     &      psno(i,j) = min(precip(i,j), hsno_max - hseff(i,j,idx,1))
#  endif
          if (tmsk(i,j) .ge. 0.5) then
!           only allow snow where there is sea ice
            psot = 0.
            do n=1,ncat
              pson = psno(i,j)*A(i,j,idx,n)
              psot = psot + pson
              hseff(i,j,idx,n) = hseff(i,j,idx,n) + fc*pson
            enddo
            psno(i,j) = psot
            if (addflux) flux(i,j,ishum) = flux(i,j,ishum)
     &                                   - dts*psno(i,j)
#  if defined uvic_mtlm
          elseif (land_map(i,j) .eq. 0) then
#  else
          else
#  endif
            hseff(i,j,idx,1) = hseff(i,j,idx,1) + fc*psno(i,j)
          endif
# elif uvic_ice
#  if defined uvic_mtlm
          hs(i,j) = hs(i,j) + hsno(i,j,2)
          if (tair .le. c0 .and. hs(i,j) .lt. hsno_max)
     &      psno(i,j) = min(precip(i,j), hsno_max - hs(i,j))
#  else
          if (tair .le. c0 .and. hsno(i,j,2) .lt. hsno_max)
     &      psno(i,j) = min(precip(i,j), hsno_max - hsno(i,j,2))
#  endif
          if (tmsk(i,j) .ge. 0.5) then
!           only allow snow where there is sea ice
            psno(i,j) = psno(i,j)*aice(i,j,2)
            hsno(i,j,2) = hsno(i,j,2) + fc*psno(i,j)
            if (addflux) flux(i,j,ishum) = flux(i,j,ishum)
     &                                   - dts*psno(i,j)
#  if defined uvic_mtlm
          elseif (land_map(i,j) .eq. 0) then
#  else
          else
#  endif
            hsno(i,j,2) = hsno(i,j,2) + fc*psno(i,j)
          endif
# endif

!-----------------------------------------------------------------------
!         update soilm and allocate surplus soil moisture to runoff
!-----------------------------------------------------------------------
# if defined uvic_mtlm
          if (tmsk(i,j) .lt. 0.5 .and. land_map(i,j) .eq. 0) then
# else
          if (tmsk(i,j) .lt. 0.5) then
# endif
            flux(i,j,ishum) = flux(i,j,ishum) - precip(i,j) + psno(i,j)
            soiltemp = soilm(i,j,2)
            soilm(i,j,2) = soilm(i,j,lf) - dts*flux(i,j,ishum)
            soilm(i,j,2) = max(c0, soilm(i,j,2))
            soilm(i,j,1) = soiltemp
            if (soilm(i,j,2) .gt. soilmax) then
              runoff(i,j) = (soilm(i,j,2) - soilmax)/dts
              soilm(i,j,2) = soilmax
            else
              runoff(i,j) = c0
            endif
          endif
        enddo
      enddo

      call embmbc (psno)
# if defined uvic_ice_cpts
      call embmbc (hseff(1,1,idx,1))
# else
      call embmbc (hsno(1,1,2))
# endif

      return
      end

      subroutine co2forc
!=======================================================================
!     calculate global average CO2 forcing

!     based on code by: A. Fanning and M. Eby
!=======================================================================
# if !defined uvic_carbon_co2_2d

      implicit none

#  include "cembm.h"
#  if defined uvic_embm_co2_lin || defined uvic_embm_co2_exp
#   include "tmngr.h"
#   include "switch.h"

      real yr

!-----------------------------------------------------------------------
!     set concentration
!-----------------------------------------------------------------------
      yr = year0 + relyr
      if (yr .gt. co2yrf) yr = co2yrf
      if (yr .lt. co2yri) yr = co2yri
#   if defined uvic_embm_co2_exp
      co2ccn = co2ccni*exp(co2rate*(yr - co2yri))
#   else
      co2ccn = co2ccni + co2rate*(yr - co2yri)
#   endif
#  endif

!-----------------------------------------------------------------------
!     current CO2 concentration is assumed to be 350 ppm
!     relative forcing (anthro) is added to current forcing
!-----------------------------------------------------------------------
      anthro = co2for*alog(co2ccn/350.0)
# endif
#endif

      return
      end
