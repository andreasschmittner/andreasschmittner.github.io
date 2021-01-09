!====================== include file "diag.h" ==========================

!     variables used for computing diagnostics:

!     tcella  = "t" cell surface area cm**2 (entire ocean)
!     ucella  = "u" cell surface area cm**2 (entire ocean)
!     tcellv  = "t" cell volume cm**3 (entire ocean)
!     ucellv  = "u" cell volume cm**3 (entire ocean)

      common /cdiag/ tcellv, ucellv, tcella(km), ucella(km)

#if defined time_step_monitor

!     ektot     = "total" kinetic energy per unit volume at "tau". units
!                 ergs/cm**3 = dyn/cm**2 = g/cm/sec**2 = 10**-7 J/cm**3.
!                 ektot is the "total" ke in the sense that it considers
!                 both the internal and external modes summed over the
!                 entire ocean volume. The contributions of
!                 vertical motions are neglected on the basis of scaling
!                 arguments (i.e., w**2 << (u**2 + v**2).
!     dtabs     = absolute value of rate of change of tracer per unit
!                 volume centered at "tau"
!     tbar      = first moment of tracer at "tau"
!     travar    = variance = second moment of tracer about mean at "tau"
!     isot1     = starting i index of section 1 for max/min overturning
!     ieot1     = ending i index of section 1 for max/min overturning
!     isot2     = starting i index of section 2 for max/min overturning
!     ieot2     = ending i index of section 2 for max/min overturning
!     jsot      = starting j index for max/min overturning
!     jeot      = ending j index for max/min overturning
!     ksot      = starting k index for max/min overturning
!     keot      = ending k index for max/min overturning
!     mrot      = regional mask region for max/min overturning
!     v_otsf    = velocity field for calculating max/min overturning
!     t_slh     = tracer fields for calculating seal level height

      common /cdiag/ ektot(0:km,jmt), dtabs(0:km,nt,jmt)
      common /cdiag/ travar(0:km,nt,jmt), tbar(0:km,nt,jmt)
      common /cdiagi/ isot1, ieot1, isot2, ieot2, jsot, jeot
      common /cdiagi/ ksot, keot, mrot
# if defined uvic_tai_otsf
      common /cdiagi/ nv_otsf
      common /cdiag/ v_otsf(jmt,km)
# endif
# if defined uvic_tai_slh
      common /cdiagi/ nt_slh
      common /cdiag/ t_slh(2:imtm1,2:jmtm1,km,2)
# endif

!     ntatio        = number of time averaged time step integrals
!     tai_ek        = average integrated kinetic energy
!     tai_t         = average integrated temperature
!     tai_s         = average integrated salinity
!     tai_tvar      = average integrated second moment of temperature
!     tai_svar      = average integrated second moment of salinity
!     tai_dt        = average integrated rate of change of temperature
!     tai_ds        = average integrated rate of change of salinity
!     tai_scan      = average scans for ocean solver
!     tai_otmax     = average maximum overturning
!     tai_otmin     = average minimum overturning
!     tai_slh       = average integrated sea level height
!     tai_hflx      = average heat flux
!     tai_sflx      = average salt flux
!     tai_dic       = average carbon
!     tai_dicflx    = average carbon flux
!     tai_alk       = average alkalinity
!     tai_o2        = average oxygen
!     tai_o2flx     = average oxygen flux
!     tai_po4       = average phosphate
!     tai_p         = average phytoplankton
!     tai_z         = average zooplankton
!     tai_d         = average detritus
!     tai_no3       = average nitrate
!     tai_di        = average diazotrophs
!     tai_c14       = average carbon 14
!     tai_dc14      = average delta carbon 14
!     tai_c14flx    = average carbon 14 flux
!     tai_cfc11     = average CFC11
!     tai_cfc11flx  = average CFC11 flux
!     tai_cfc12     = average CFC12
!     tai_cfc12flx  = average CFC12 flux

      common /cdiagi/ ntatio
      common /cdiag/ tai_ek, tai_t, tai_s, tai_tvar, tai_svar, tai_dt
      common /cdiag/ tai_ds, tai_scan, tai_otmax, tai_otmin, tai_slh
      common /cdiag/ tai_hflx, tai_sflx, tai_dic, tai_dicflx, tai_alk
      common /cdiag/ tai_o2, tai_o2flx, tai_po4, tai_p, tai_z, tai_d
      common /cdiag/ tai_no3, tai_di, tai_c14, tai_dc14, tai_c14flx
      common /cdiag/ tai_cfc11, tai_cfc11flx, tai_cfc12, tai_cfc12flx
#endif
#if defined energy_analysis

!     engint    = volume averaged internal mode energy integral
!                 components
!     engext    = volume averaged external mode energy integral
!                 components
!     buoy      = volume averaged buoyancy

!     tcerr     = maximum "t" cell continuity error
!     ucerr     = maximum "u" cell continuity error
!     itcerr    = "i" index corresponding to "tcerr"
!     jtcerr    = "jrow" index corresponding to "tcerr"
!     ktcerr    = "k" index corresponding to "tcerr"
!     iucerr    = "i" index corresponding to "ucerr"
!     jucerr    = "jrow" index corresponding to "ucerr"
!     kucerr    = "k" index corresponding to "ucerr"

!     wtbot     = maximum "adv_vbt" error at ocean bottom
!     iwtbot    = "i" index corresponding to "wtbot"
!     jwtbot    = "jrow" index corresponding to "wtbot"
!     kwtbot    = "k" index corresponding to "wtbot"
!     wubot     = maximum "adv_vbu" at ocean bottom
!     iwubot    = "i" index corresponding to "wubot"
!     jwubot    = "jrow" index corresponding to "wubot"
!     kwubot    = "k" index corresponding to "wubot"

!     wtlev     = zonally integrated adv_vbt for each level
!     wulev     = zonally integrated adv_vbu for each level

      common /cdiag/ buoy(0:km,jmt), engint(0:km,8,jmt), engext(8,jmt)
      common /cdiag/ tcerr(jmt), ucerr(jmt)
      common /cdiagi/ itcerr(jmt), jtcerr(jmt), ktcerr(jmt)
      common /cdiagi/ iucerr(jmt), jucerr(jmt), kucerr(jmt)
      common /cdiag/ wtbot(jmt), wubot(jmt)
      common /cdiagi/ iwtbot(jmt), jwtbot(jmt), kwtbot(jmt)
      common /cdiagi/ iwubot(jmt), jwubot(jmt), kwubot(jmt)
      common /cdiag/ wtlev(km,0:jmt), wulev(km,0:jmt)
#endif
#if defined term_balances
# include "termbal.h"
#endif
#if defined gyre_components

!     ttn       = northward transport of tracer components

!     ttn2      = northward transport of tracers for ocean basins
!                  (.,.,.,0)       Global
!                  (.,.,.,1:nhreg) Ocean basins
!                also,
!                  (6,.,.,.) total transport due to advection
!                  (7,.,.,.) total transport due to diffusion
!                  (8,.,.,.) total transport
# if defined isopycmix && defined gent_mcwilliams && !defined fct && !defined quicker
!                  (9,.,.,.) total transport due to isopycnal advection
# endif
      common /gyres/ ttn(8,jmt,ntmin2)
# if defined isopycmix && defined gent_mcwilliams && !defined fct && !defined quicker
      common /gyres/ ttn2(6:9,jmt,nt,0:nhreg)
# else
      common /gyres/ ttn2(6:8,jmt,nt,0:nhreg)
# endif
#endif
#if defined meridional_overturning

!     vmsf      = vertical_meridional stream function

      common /cdiag/ vmsf(jmt,km)
#endif
#if defined show_zonal_mean_of_sbc

!     zmsmf     = zonal mean surface momentum flux
!     zmstf     = zonal mean surface tracer flux
!     zmsm      = zonal mean surface momentum
!     zmst      = zonal mean surface tracers
!     zmau      = surface area weighting for "u" latitudes
!     zmat      = surface area weighting for "t" latitudes

      common /cdiag/ zmsmf(jmt,2), zmstf(jmt,nt), zmau(jmt), zmat(jmt)
      common /cdiag/ zmsm(jmt,2), zmst(jmt,nt)
#endif

#if defined tracer_yz

!     tyz(,,,1) = zonal mean of tracer T
!     tyz(,,,2) = zonal mean of d(T)/dt
!     tyz(,,,3) = zonal mean of advection of T
!     tyz(,,,4) = zonal mean of diffusion of T
!     tyz(,,,5) = zonal mean of source of T

      common /cdiag/ tyz(jmt,km,nt,5)
#endif
