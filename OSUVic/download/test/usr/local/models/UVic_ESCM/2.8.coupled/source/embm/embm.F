      subroutine embm (is, ie, js, je)

!=======================================================================

!     The atmospheric energy moisture balance model (EMBM) is based
!     primarily on the work of Augustus Fanning.

!     Fanning, A.F. and A.J. Weaver, An atmospheric energy-moisture
!       balance model: climatology, interpentadal climate change,
!       and coupling to an ocean general circulation model,
!       J. Geophys. Res., 101, 15,111-15,128, 1996

!     based on code by: A. Fanning and M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "size.h"
# include "csbc.h"
# include "cembm.h"
# if !defined uvic_mom
#  include "tmngr.h"
# endif

      integer ie, is, je, js, n

# if !defined uvic_mom
!-----------------------------------------------------------------------
!     update timestep counter, set corresponding model time, and set
!     time dependent logical switches which determine program flow.
!-----------------------------------------------------------------------

        itt = itt + 1
        call increment_time (dtatm)
        call set_time_switches

# endif
!-----------------------------------------------------------------------
!     increment counter and set the time step type
!-----------------------------------------------------------------------

      nats = nats + 1
      if (nats .gt. namix) then
        lf = 2
        dts = dtatm
        nats = 1
      else
        lf = 1
        dts = 2.0*dtatm
      endif

# if defined uvic_ice_cpts

!-----------------------------------------------------------------------
!     initialize a few things before each call to embm and ice
!-----------------------------------------------------------------------

      call init_ice_cpts (is, ie, js, je)
# endif

      addflux = .true.
# if defined uvic_embm_even_fluxes && defined uvic_mom
      if (mod(nats,2) .ne. 0) addflux = .false.
# endif

!-----------------------------------------------------------------------
!     calculate fluxes at tau
!-----------------------------------------------------------------------

      call fluxes (is, ie, js, je)
# if defined uvic_ice

!-----------------------------------------------------------------------
!     compute ice fluxes at tau and ice thickness and area at tau+1
!-----------------------------------------------------------------------

      call ice (is, ie, js, je)
# endif

!-----------------------------------------------------------------------
!     compute atmospheric tracers at tau+1. start with humidity so that
!     the precipitation flux can be calculated for latent heat
!-----------------------------------------------------------------------

      call solve (2)
      call precipitate (is, ie, js, je)
      call solve (1)
      do n=3,nat
        call solve (n)
      enddo

!-----------------------------------------------------------------------
!     calculate the total atmospheric fluxes for coupling
!-----------------------------------------------------------------------

      if (addflux) call sum_flux (is, ie, js, je)
# if defined time_averages

!-----------------------------------------------------------------------
!     accumulate time averages
!-----------------------------------------------------------------------
      call ta_embm_snap (is, ie, js, je, 1)
# endif
# if defined time_step_monitor

!-----------------------------------------------------------------------
!     accumulate time step integrals
!-----------------------------------------------------------------------
      call ta_embm_tsi (1)
# endif
# if !defined uvic_mom

!-----------------------------------------------------------------------
!     write output
!-----------------------------------------------------------------------

      call embmout (is, ie, js, je)
# endif
#endif

      return
      end

      subroutine sum_flux (is, ie, js, je)

!=======================================================================
!     sum fluxes over atmospheric time steps

!     based on code by:  A. Fanning and M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "param.h"
# include "cembm.h"
# include "atm.h"
# include "ice.h"
# include "csbc.h"
# if defined uvic_global_sums
#  include "coord.h"
#  include "grdvar.h"
# endif

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1
# if defined uvic_global_sums

      real fa, fb
      fb = scatter*(1. + pass)
# endif

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      do j=jsp1,jem1
# if defined uvic_global_sums
        fa = dts*cst(j)*dyt(j)*1.0e-7
# endif
        do i=isp1,iem1
          if (tmsk(i,j) .ge. 0.5) then
            flux(i,j,isat) = flux(i,j,isat) + dts*(dnswr(i,j)
     &                     - uplwr(i,j) - upltnt(i,j) - upsens(i,j))
            flux(i,j,ishum) = flux(i,j,ishum) + dts*(precip(i,j)
     &                      - evap(i,j))
# if defined uvic_mtlm
          elseif (land_map(i,j) .ne. 0) then
            sbc(i,j,iat) = sbc(i,j,iat) + dts*(at(i,j,2,1)
#  if defined uvic_embm_icedata
     &                   - hicel(i,j,2)*rlapse
#  endif
     &                   - elev(i,j)*rlapse)
            sbc(i,j,irh) = sbc(i,j,irh) + dts*rh(i,j)
            sbc(i,j,iaws) = sbc(i,j,iaws) + dts*sbc(i,j,iws)
            sbc(i,j,iswr) = sbc(i,j,iswr) + dts*dnswr(i,j)
            sbc(i,j,ipr) = sbc(i,j,ipr) + dts*precip(i,j)
            if (psno(i,j) .ge. 0.) then
              sbc(i,j,ips) = sbc(i,j,ips) + dts*psno(i,j)
              sbc(i,j,ipr) = sbc(i,j,ipr) - dts*psno(i,j)
            endif
# endif
          else
            sbc(i,j,iro) = sbc(i,j,iro) + dts*runoff(i,j)
          endif
# if defined  uvic_ice_evp || defined uvic_embm_awind
          if (umsk(i,j) .ge. 0.5) then
            flux(i,j,nat+1) = flux(i,j,nat+1) + dts*sbc(i,j,itaux)
            flux(i,j,nat+2) = flux(i,j,nat+2) + dts*sbc(i,j,itauy)
#  if defined uvic_ice_evp
            flux(i,j,nat+1) = flux(i,j,nat+1) + dts*xint(i,j)
            flux(i,j,nat+2) = flux(i,j,nat+2) + dts*yint(i,j)
#  endif
          endif
# endif
# if defined uvic_embm_running_average || defined uvic_embm_awind
          atbar(i,j) = atbar(i,j) + dts*at(i,j,2,isat)
# endif
# if defined uvic_global_sums
          dtoih = dtoih + fa*dxt(i)*(outlwr(i,j) - dnswr(i,j)*pass
     &          - solins(i,j)*sbc(i,j,iaca)*fb)
# endif
        enddo
      enddo

      totaltime = totaltime + dts
      atatm = atatm + dts
#endif

      return
      end
