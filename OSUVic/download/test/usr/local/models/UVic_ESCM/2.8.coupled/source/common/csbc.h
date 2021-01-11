!======================= include file "csbc.h" =========================
#include "derived_options.h"
!     surface boundary conditions (S.B.C.)

!     numsbc    = total number of surface boundary conditions
!     indices and names set in UVic_ESCM.F:
!       itaux     = x component of wind stress (dynes cm-2)
!       itauy     = y component of wind stress (dynes cm-2)
!       iws       = surface wind speed (cm s-1)
!       iaca      = atmospheric coalbedo
!       isca      = surface coalbedo (%/100)
!       ihflx     = heat flux
!       isflx     = salt flux
!       isst      = ocean model SST (C)
!       isss      = ocean model SSS (PSU*0.001-0.035)
!       iwa       = surface wind angle (degrees)
!       iro       = surface runoff (g cm-2 s-1)
!       iwxq      = x component of wind for moisture advection (cm s-1)
!       iwyq      = y component of wind for moisture advection (cm s-1)
!       iwxt      = x component of wind for temperature advection (cm s-1)
!       iwyt      = y component of wind for temperature advection (cm s-1)
!       ipsw      = penetrating shortwave (cal cm-2 s-1)
!       isu       = x component of ocean first layer velocity (cm s-1)
!       isv       = y component of ocean first layer velocity (cm s-1)
!       igu       = x component of ocean second layer velocity (cm s-1)
!       igv       = y component of ocean second layer velocity (cm s-1)
!       issdic    = sea surface concentration of dic (umol cm-3)
!       idicflx   = sea surface flux of carbon (umol cm-2 s-1)
!       issalk    = sea surface concentration of alkalinity (umol cm-3)
!       ialkflx   = sea surface flux of alkalinity (umol cm-2 s-1)
!       isso2     = sea surface concentration of  (umol cm-3)
!       io2flx    = sea surface flux of  (umol cm-2 s-1)
!       isspo4    = sea surface concentration of phosphate (nmol cm-3)
!       ipo4flx   = sea surface flux of phosphate (nmol cm-2 s-1)
!       issphyt   = sea surface concentration of phytoplankton (nmol P cm-3)
!       iphytflx  = sea surface flux of phytoplankton (nmol P cm-2 s-1)
!       isszoop   = sea surface concentration of zooplankton (nmol P cm-3)
!       izoopflx  = sea surface flux of zooplankton (nmol P cm-2 s-1)
!       issdetr   = sea surface concentration of detritus (nmol P cm-3)
!       idetrflx  = sea surface flux of detritus (nmol P cm-2 s-1)
!       issno3    = sea surface concentration of nitrate (nmol P cm-3)
!       ino3flx   = sea surface flux of nitrate (nmol P cm-2 s-1)
!       issdiaz   = sea surface concentration of diazotraphs (nmol P cm-3)
!       idiazflx  = sea surface flux of diazotraphs (nmol P cm-2 s-1)
!       issc14    = sea surface concentration of carbon 14 (umol cm-3)
!       ic14flx   = sea surface flux of carbon 14 (umol cm-2 s-1)
!       isscfc11  = sea surface concentration of cfc11 (umol cm-3)
!       icfc11flx = sea surface flux of cfc11 (umol cm-2 s-1)
!       isscfc12  = sea surface concentration of cfc12 (umol cm-3)
!       icfc12flx = sea surface flux of cfc12 (umol cm-2 s-1)
!       iat       = surface air temperature (C)
!       irh       = surface relative humidity (%/100)
!       ipr       = precipitation as rain (g cm-2 s-1)
!       ips       = precipitation as snow (g cm-2 s-1)
!       iaws      = averaged surface wind speed (cm s-1)
!       iswr      = surface shortwave radiation (g s-3)
!       ilwr      = surface longwave radiation (g s-3)
!       isens     = surface sensible heat flux (g s-3)
!       ievap     = surface evaporation (g cm-2 s-1)
!       idtr      = diurnal temperature range (C)
!       inpp      = net primary production
!       isr       = soil respiration
!     mapsbc      = surface boundary conditions names
!     sbc         = surface boundary condition data.

!     socn       = average ocean sea surface salinity
!     dicocn     = average ocean sea surface dic
!     alkocn     = average ocean sea surface alkalinity
!     o2ocn      = average ocean sea surface oxygen
!     po4ocn     = average ocean sea surface phosphate
!     phytocn    = average ocean sea surface phytoplankton
!     zoopocn    = average ocean sea surface zooplankton
!     detrocn    = average ocean sea surface detritus
!     no3ocn     = average ocean sea surface nitrogen
!     diazocn    = average ocean sea surface diazotrophs
!     c14ocn     = average ocean sea surface c14
!     cfc11ocn   = average ocean sea surface cfc11
!     cfc12ocn   = average ocean sea surface cfc12

!     ntspas      = the number of ocean time steps per atmosphere segment
!     ntspls      = the number of ocean time steps per land segment
!     ntspos      = the number of ocean time steps per ocean segment
!     dtatm       = atmosphere time step (s)
!     dtlnd       = land time step (s)
!     dtocn       = ocean time step (s)
!     dtatm       = atmosphere boundary condition averaging time (s)
!     dtlnd       = land boundary condition averaging time (s)
!     dtocn       = ocean boundary condition averaging time (s)
!     dampts      = time scale for damping surface tracers (days) to data
!     dampdz      = ocean level thickness for converting Newtonian damping
!                   to a surface flux
!     land_map    = map with indices for coupling to land arrays

      integer numsbc
      parameter (numsbc = 10
#if defined uvic_embm_awind
     &                  + 1
#endif
#if defined uvic_embm_adv_q
     &                  + 2
#endif
#if defined uvic_embm_adv_t
     &                  + 2
#endif
#if defined shortwave
     &                  + 1
#endif
#if defined uvic_ice_evp
     &                  + 4
#endif
#if defined uvic_carbon
     &                  + 2
# if defined uvic_carbon_14
     &                  + 2
# endif
#endif
#if defined uvic_alk
     &                  + 2
#endif
#if defined uvic_o2
     &                  + 2
#endif
#if defined uvic_cfc11
     &                  + 2
#endif
#if defined uvic_cfc12
     &                  + 2
#endif
#if defined uvic_npzd
# if defined uvic_npzd_vflux
     &                  + 8
# else
     &                  + 2
# endif
#endif
#if defined uvic_nitrogen
# if defined uvic_npzd_vflux
     &                  + 4
# else
     &                  + 2
# endif
#endif
#if defined uvic_mtlm
     &                  + 10
#endif
#if defined uvic_mtlm && defined uvic_carbon
     &                  + 2
#endif
     &                     )

      integer itaux, itauy, iws, iaca, isca, ihflx, isflx, isst, isss
      integer iwa, iro, iwxq, iwyq, iwxt, iwyt, ipsw, isu, isv, igu
      integer igv, issdic, idicflx, issalk, ialkflx, isso2, io2flx
      integer isspo4, ipo4flx, issphyt, iphytflx, isszoop, izoopflx
      integer issdetr, idetrflx, issno3, ino3flx, issdiaz, idiazflx
      integer issc14, ic14flx, isscfc11, icfc11flx, isscfc12, icfc12flx
      integer iat, irh, ipr, ips, iaws, iswr, ilwr, isens, ievap, idtr
      integer inpp, isr

      common /csbc_i/ itaux, itauy, iws, iaca, isca, ihflx, isflx, isst
      common /csbc_i/ isss, iwa, iro, iwxq, iwyq, iwxt, iwyt, ipsw, isu
      common /csbc_i/ isv, igu, igv, issdic, idicflx, issalk, ialkflx
      common /csbc_i/ isso2, io2flx, isspo4, ipo4flx, issphyt, iphytflx
      common /csbc_i/ isszoop, izoopflx, issdetr, idetrflx, issno3
      common /csbc_i/ ino3flx, issdiaz, idiazflx, issc14, ic14flx
      common /csbc_i/ isscfc11, icfc11flx, isscfc12, icfc12flx, iat
      common /csbc_i/ irh, ipr, ips, iaws, iswr, ilwr, isens, ievap
      common /csbc_i/ idtr, inpp, isr

      character(20) :: mapsbc
      common /csbc_c/ mapsbc(numsbc)

      real sbc
      common /csbc_r/ sbc(imt,jmt,numsbc)

      real socn, dicocn, alkocn, o2ocn, po4ocn, phytocn, zoopocn
      real detrocn, no3ocn, diazocn, c14ocn, cfc11ocn, cfc12ocn

      common /csbc_r/ socn, dicocn, alkocn, o2ocn, po4ocn, phytocn
      common /csbc_r/ zoopocn, detrocn, no3ocn, diazocn, c14ocn
      common /csbc_r/ cfc11ocn, cfc12ocn

      integer ntspas, ntspls, ntspos
      common /csbc_i/ ntspas, ntspls, ntspos

      real dtatm, dtlnd, dtocn, atatm, atlnd, atocn
      common /csbc_r/ dtlnd, dtatm, dtocn, atatm, atlnd, atocn

      real dampts, dampdz
      common /csbc_r/ dampts(nt), dampdz(nt)
#if defined llnl_plume
      common /csbc_r/ subflux(imt,jmt,nt), subz(imt,jmt)
#endif
#if defined uvic_mtlm
      integer land_map

      common /csbc_i/ land_map(imt,jmt)
#endif
