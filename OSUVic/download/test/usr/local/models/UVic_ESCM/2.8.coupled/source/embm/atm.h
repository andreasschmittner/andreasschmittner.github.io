!======================== include file "atm.h" =========================
#include "derived_options.h"
!     arrays for the energy-moisture balance model

!     note: units for heat flux are: g s-3 (or mW m-2)
!           units for fresh water flux are in g cm-2 s-1
!           downward is into the top surface (ocean, ice or land)
!           upward is into the bottom of the atmosphere
!           outward is out of the top of the atmosphere
!           inward is into the top of the atmosphere

!     at        = tracers (previous and most recent)
!     elev      = elevations (cm)
!     flux      = downward flux accumulator (1=heat, 2=fresh water, ...)
!     rh        = relative humidity (%/100)
!     tmsk      = tracer grid ocean mask (0 = land, 1 = ocean)
!     umsk      = velocity grid ocean mask (0 = land, 1 = ocean)
!     dn        = northward atmospheric tracer diffusivity
!     de        = eastward atmospheric tracer diffusivity
!     fcor      = Coriolis factor
!   heat fluxes
!     solins    = incoming short wave radiation (g s-3)
!     dnswr     = downward surface shortwave flux absorbed (g s-3)
!     uplwr     = upward surface longwave flux (g s-3)
!     upsens    = upward surface sensible heat flux (g s-3)
!     upltnt    = upward surface latent heat flux (g s-3)
!     outlwr    = outgoing atmosphere longwave flux (g s-3)
!   fresh water fluxes
!     precip    = precipitation (g cm-2 s-1)
!     evap      = evaporation (g cm-2 s-1)
!     disch     = discharge from continents (g cm-2 s-1)
!   land model
!     soilm     = soil moisture, as depth in bucket (g cm -2)
!     runoff    = water runoff from continents (g cm-2 s-1)
!     surf      = land surface temperature (C)
!   annual average temperature ("uvic_embm_running_average")
!     rtbar     = running average atmospheric temperature (C)
!     atbar     = average temperature accumulator (C)
!   wind feedback ("uvic_embm_awind")
!     tbar      = climatological average atmospheric temperature (C)
!     awx       = anomalous x component of wind (cm s-1)
!     awy       = anomalous y component of wind (cm s-1)
!     apress    = anomalous pressure (g cm-2 s-2)
!   flux adjustment ("uvic_save_flxadj")
!     flxadj    = flux adjustment
!   indicies for atmospheric tracer array
!     isat      = index for surface air temperature
!     ishum     = index for surface specific humidity
!     ico2      = index for atmospheric co2
!     mapat     = map for atmospheric tracer names

      character(10) :: mapat(nat)
      common /embm_c/ mapat

      integer isat, ishum, ico2
      common /embm_i/ isat, ishum, ico2

      real at, elev, flux, rh, tmsk, umsk, dn, de, fcor, solins, dnswr
      real uplwr, upsens, upltnt, outlwr, precip, evap, disch, soilm
      real runoff, surf, rtbar, atbar, tbar, awx, awy, apress, flxadj

      common /embm_r/ at(imt,jmt,2,nat), elev(imt,jmt)
#if defined uvic_ice_evp || defined uvic_embm_awind
      common /embm_r/ flux(imt,jmt,nat+2), rh(imt,jmt)
#else
      common /embm_r/ flux(imt,jmt,nat), rh(imt,jmt)
#endif
      common /embm_r/ tmsk(imt,jmt), umsk(imt,jmt), dn(imt,jmt,nat)
      common /embm_r/ de(imt,jmt,nat), fcor(imt,jmt), solins(imt,jmt)
      common /embm_r/ dnswr(imt,jmt), uplwr(imt,jmt), upsens(imt,jmt)
      common /embm_r/ upltnt(imt,jmt), outlwr(imt,jmt), precip(imt,jmt)
      common /embm_r/ evap(imt,jmt), disch(imt,jmt)
      common /embm_r/ soilm(imt,jmt,2), runoff(imt,jmt)
      common /embm_r/ surf(imt,jmt)
#if defined uvic_embm_running_average || defined uvic_embm_awind
      common /embm_r/ rtbar(imt,jmt), atbar(imt,jmt)
#endif
#if defined uvic_embm_awind
      common /embm_r/ tbar(imt,jmt), awx(imt,jmt), awy(imt,jmt)
      common /embm_r/ apress(imt,jmt)
#endif
#if defined uvic_save_flxadj
      common /embm_r/ flxadj(imt,jmt,2)
#endif

#if defined time_averages
!   time averaged arrays
!     ta_at         = time averaged atmospheric tracers
!     ta_solins     = time averaged incoming shortwave flux
!     ta_arswr      = time averaged atmospheric reflected shortwave flux
!     ta_dnswr      = time averaged downward shortwave flux
!     ta_absin      = time averaged absorbed shortwave coming in
!     ta_absout     = time averaged absorbed shortwave going out
!     ta_uplwr      = time averaged upward longwave flux
!     ta_upsens     = time averaged upward sensible heat flux
!     ta_upltnt     = time averaged upward latent heat flux
!     ta_outlwr     = time averaged outgoing longwave flux
!     ta_precip     = time averaged precipitation
!     ta_psno       = time averaged precipitation as snow
!     ta_evap       = time averaged evaporation
!     ta_disch      = time averaged discharge
!     ta_ws         = time averaged surface wind speed
!     ta_runoff     = time averaged runoff
!     ta_soilm      = time averaged soil moisture
!     ta_surf       = time averaged land surface temperature
!   moisture advection ("uvic_embm_adv_q")
!     ta_wx_q       = time averaged x component of wind
!     ta_wy_q       = time averaged y component of wind
!   moisture advection ("uvic_embm_adv_t")
!     ta_wx_t       = time averaged x component of wind
!     ta_wy_t       = time averaged y component of wind
!   windstress feedback ("uvic_embm_awind")
!     ta_awx        = time averaged x component of wind
!     ta_awy        = time averaged y component of wind
!     ta_rtbar      = time averaged running average temperature
!     ta_apress     = time averaged anomalous pressure
!   flux adjustment ("uvic_save_flxadj")
!     ta_flxadj     = time averaged flux adjustment
!   diffusivities ("uvic_embm_save_diff")
!     ta_dn         = time averaged dn
!     ta_de         = time averaged de
!   flux of co2 ("uvic_carbon and uvic_carbon_co2_2d")
!     ta_flxco2     = time averaged flxco2
!   uncoupled
!     ta_hflux      = time averaged heat flux
!     ta_sflux      = time averaged salt flux
!     ta_sss        = time averaged sea surface salinity
!     ta_sst        = time averaged sea surface temperature
!     ta_taux       = time averaged x component of wind stress
!     ta_tauy       = time averaged x component of wind stress

      real ta_at, ta_solins, ta_arswr, ta_dnswr, ta_absin, ta_absout
      real ta_uplwr, ta_upsens, ta_upltnt, ta_outlwr, ta_precip
      real ta_psno, ta_evap, ta_disch, ta_ws, ta_soilm, ta_runoff
      real ta_surf, ta_wx_q, ta_wy_q, ta_wx_t, ta_wy_t, ta_awx, ta_awy
      real ta_rtbar, ta_apress, ta_flxadj, ta_dn, ta_de, ta_flxco2
      real ta_hflux, ta_sflux, ta_sss, ta_sst, ta_taux, ta_tauy

      common /embm_r/ ta_at(imt,jmt,nat), ta_solins(imt,jmt)
      common /embm_r/ ta_arswr(imt,jmt), ta_dnswr(imt,jmt)
      common /embm_r/ ta_absin(imt,jmt), ta_absout(imt,jmt)
      common /embm_r/ ta_uplwr(imt,jmt), ta_upsens(imt,jmt)
      common /embm_r/ ta_upltnt(imt,jmt), ta_outlwr(imt,jmt)
      common /embm_r/ ta_precip(imt,jmt), ta_psno(imt,jmt)
      common /embm_r/ ta_evap(imt,jmt), ta_disch(imt,jmt)
      common /embm_r/ ta_ws(imt,jmt), ta_runoff(imt,jmt)
      common /embm_r/ ta_soilm(imt,jmt), ta_surf(imt,jmt)
# if defined uvic_embm_adv_q
      common /embm_r/ ta_wx_q(imt,jmt), ta_wy_q(imt,jmt)
# endif
# if defined uvic_embm_adv_t
      common /embm_r/ ta_wx_t(imt,jmt), ta_wy_t(imt,jmt)
# endif
# if defined uvic_embm_awind
      common /embm_r/ ta_awx(imt,jmt), ta_awy(imt,jmt)
      common /embm_r/ ta_rtbar(imt,jmt), ta_apress(imt,jmt)
# endif
# if defined uvic_save_flxadj
      common /embm_r/ ta_flxadj(imt,jmt,2)
# endif
# if defined uvic_embm_save_diff
      common /embm_r/ ta_dn(imt,jmt,nat), ta_de(imt,jmt,nat)
# endif
# if defined uvic_carbon && defined uvic_carbon_co2_2d
      common /embm_r/ ta_flxco2(imt,jmt)
# endif
# if !defined uvic_mom
      common /embm_r/ ta_hflux(imt,jmt), ta_sflux(imt,jmt)
      common /embm_r/ ta_sss(imt,jmt), ta_sst(imt,jmt)
      common /embm_r/ ta_taux(imt,jmt), ta_tauy(imt,jmt)
# endif
#endif
