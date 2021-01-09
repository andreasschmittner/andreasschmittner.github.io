!======================= include file "cembm.h" ========================

!     parameters for use in the energy balance model (also see atm.h)

!     addflux    = logical flag for adding only even mode fluxes
!     nats       = number of atmospheric time steps since mixing
!     namix      = time steps between mixing (set in atmos.in)
!     lf         = time step flag (1=>leapfrog, 2=>forward)
!     niats      = number of ice advection sub time steps
!     nivts      = time steps between recalculating ice velocities
!     nivc       = time step counter for nivts
!     ns         = number of subcycles for advection and diffusion
!     pyear      = default paleo calendar year (-/+ = BC/AD)
!     dts        = time step (2*dtatm=>leapfrog, dtatm=>forward)
!     co2ccn     = atmospheric CO2 concentration (ppmv)
!     co2flx     = atmospheric CO2 emissions flux (umol cm-2 s-1)
!     co2ccni    = initial atmospheric CO2 concentration
!     anthro     = radiative forcing by atmospheric CO2
!     co2yri     = last year of initial atmospheric CO2 concentration
!     co2yrf     = first year of final atmospheric CO2 concentration
!     co2rate    = atmospheric CO2 increase rate (linear or exponential)
!     co2for     = atmospheric CO2 forcing term (in units of heat flux)
!     c14ccn     = atmospheric C14 concentration (ppmv)
!     dc14ccn    = atmospheric dC14 concentration (permil)
!     dc14ccnn   = northern hemisphere atmospheric dC14 concentration
!     dc14ccne   = equatorial atmospheric dC14 concentration
!     dc14ccns   = southern hemisphere atmospheric dC14 concentration
!     c14prod    = atmospheric C14 production (umol cm-2 s-1)
!     cfc11ccnn  = northern hemisphere atmospheric CFC11 concentration
!     cfc11ccns  = southern hemisphere atmospheric CFC11 concentration
!     cfc12ccnn  = northern hemisphere atmospheric CFC12 concentration
!     cfc12ccns  = southern hemisphere atmospheric CFC12 concentration
!     scatter    = proportion of solar scattered by the atmosphere
!     solarconst = solar constant (g/s**3)
!     cssh       = constant used in calculation of ssh (g/g)
!     cdatm      = drag coefficient (dimensionless)
!     cpatm      = atmospheric heat capacity (cm**2/s**2/K)
!     sht        = scale height for temperature
!     shq        = scale height for specific humidity
!     shc        = scale height for carbon
!     rhoatm     = density of air at sea surface (g/cm**3)
!     esatm      = atmosphere emissivity times Stefan's constant
!     rhoocn     = representative sea surface density
!     esocn      = ocean emissivity times Stefan's constant
!     vlocn      = latent heat of vapourization of water
!     cdice      = drag coefficient (dimensionless)
!     dampice    = time scale for freezing first layer under ice (days)
!     rhoice     = ice density (g/cm**3)
!     rhosno     = snow density (g/cm**3)
!     esice      = ice emissivity times Stefan's constant
!     slice      = latent heat of sublimation of ice
!     flice      = latent heat of fusion of ice (cm**2/s**2)
!     condice    = ice conductivity (g*cm/s**3/K)
!     tsno       = air temperature for accumulating snow
!     hsno_max   = maximum snow depth
!     totaltime  = total time for long term averages
!     dtoih      = total net heat flux through top of atmosphere
!     rlapse     = lapse rate (K/cm)
!     soilmax    = soil water field capacity (cm)
!     eslnd      = land emissivity time Stefan's constant
!     pass       = atmospheric transmission coefficient (%/100)
!     cs_alb     = clear sky atmospheric albedo (%/100)
!     ice_calb   = ice coalbedo (%/100)
!     sno_calb   = snow coalbedo (%/100)
!     pcfactor   = precip - cloud correlation factor ( %/100)
!     rfactor    = lapse rate reduction factor ( %/100)
!     orbit_yr   = orbital forcing year (-/+ = BC/AD)
!     solar_yr   = solar forcing year (-/+ = BC/AD)
!     c14_yr     = c14 forcing year (-/+ = BC/AD)
!     co2_yr     = co2 forcing year (-/+ = BC/AD)
!     ice_yr     = ice sheet forcing year (-/+ = BC/AD)
!     dalt_v     = dalton number over vegetation (no vegetation model)
!     dalt_o     = dalton number over ocean
!     dalt_i     = dalton number over ice
!     rhmax      = maximum relative humidity

      logical addflux

      common /cembm_l/ addflux

      integer nats, namix, lf, niats, nivts, nivc, ns

      common /cembm_i/ nats, namix, lf, niats, nivts, nivc, ns

      real pyear, dts, co2ccn, co2flx,co2ccni, anthro, co2yri, co2yrf
      real co2rate, co2for, c14ccn, dc14ccn, dc14ccnn, dc14ccne
      real dc14ccns, c14prod, cfc11ccnn, cfc11ccns, cfc12ccnn
      real cfc12ccns, scatter, solarconst, cssh, cdatm, cpatm, sht, shq
      real shc, rhoatm, esatm, rhoocn, esocn, vlocn, cdice, dampice
      real rhoice, rhosno, esice,slice, flice, condice, tsno,  hsno_max
      real totaltime, dtoih, rlapse, soilmax, eslnd, pass, cs_alb
      real ice_calb, sno_calb, pcfactor, rfactor, orbit_yr, solar_yr
      real c14_yr, co2_yr, ice_yr, dalt_v, dalt_o, dalt_i, rhmax

      common /cembm_r/ pyear, dts, co2ccn, co2flx,co2ccni, anthro
      common /cembm_r/ co2yri,  co2yrf, co2rate, co2for, c14ccn, dc14ccn
      common /cembm_r/ dc14ccnn, dc14ccne, dc14ccns, c14prod, cfc11ccnn
      common /cembm_r/ cfc11ccns, cfc12ccnn, cfc12ccns, cssh, cdatm
      common /cembm_r/ cpatm, sht, shq, shc, rhoatm, esatm, rhoocn
      common /cembm_r/ scatter, solarconst, esocn, dampice, rhoice
      common /cembm_r/ rhosno, esice, slice, flice, condice, vlocn
      common /cembm_r/ cdice, tsno, hsno_max, totaltime, dtoih, rlapse
      common /cembm_r/ soilmax, eslnd, pass, cs_alb, ice_calb, sno_calb
      common /cembm_r/ pcfactor, rfactor,  orbit_yr, solar_yr, c14_yr
      common /cembm_r/ co2_yr, ice_yr, dalt_v, dalt_o, dalt_i, rhmax

!     ntatsa        = time step counter for time averaging
!     ntatia        = number of time averaged time step integrals
!     tai_sat       = average integrated elevated surface air temperature
!     tai_shum      = average integrated surface specific humidity
!     tai_precip    = average integrated precipitation
!     tai_evap      = average integrated evaporation
!     tai_ohice     = total integrated sea ice volume
!     tai_oaice     = total integrated sea ice area
!     tai_hsno      = total integrated snow volume
!     tai_lhice     = total integrated land ice volume
!     tai_laice     = total integrated land ice area (includes snow)
!     tai_co2ccn    = average integrated CO2 concentration
!     tai_dc14ccn    = average integrated dC14 concentration
!     tai_cfc11ccn  = average integrated CFC11 concentration
!     tai_cfc12ccn  = average integrated CFC12 concentration
!     tai_maxit     = average maximum iterations for atmospheric solver
!     tai_sst       = average sea surface temperature
!     tai_sss       = average sea surface salinity
!     tai_nsat      = average northern hemisphere surface air temperature
!     tai_ssat      = average southern hemisphere surface air temperature
!     tai_nshum     = average northern hemisphere surface specific humidity
!     tai_sshum     = average southern hemisphere surface specific humidity
!     tai_nprecip   = average northern hemisphere precipitation
!     tai_sprecip   = average southern hemisphere precipitation
!     tai_nevap     = average northern hemisphere evaporation
!     tai_sevap     = average southern hemisphere evaporation
!     tai_nohice    = total northern hemisphere sea ice volume
!     tai_sohice    = total southern hemisphere sea ice volume
!     tai_noaice    = total northern hemisphere sea ice area
!     tai_soaice    = total southern hemisphere sea ice area
!     tai_nhsno     = total northern hemisphere snow volume
!     tai_shsno     = total southern hemisphere snow volume
!     tai_nlhice    = total northern hemisphere land ice volume
!     tai_slhice    = total southern hemisphere land ice volume
!     tai_nlaice    = total northern hemisphere land ice area
!     tai_slaice    = total southern hemisphere land ice area

      integer ntatsa, ntatia

      common /cembm_i/ ntatsa, ntatia

      real tai_sat, tai_shum, tai_precip, tai_evap, tai_ohice, tai_oaice
      real tai_hsno, tai_lhice, tai_laice, tai_co2ccn, tai_dc14ccn
      real tai_cfc11ccn, tai_cfc12ccn, tai_maxit, tai_sst, tai_sss
      real tai_nsat, tai_ssat, tai_nshum, tai_sshum, tai_nprecip
      real tai_sprecip, tai_nevap, tai_sevap, tai_nohice, tai_sohice
      real tai_noaice, tai_soaice, tai_nhsno, tai_shsno, tai_nlhice
      real tai_slhice, tai_nlaice, tai_slaice

      common /cembm_r/ tai_sat, tai_shum, tai_precip, tai_evap
      common /cembm_r/ tai_ohice, tai_oaice, tai_hsno, tai_lhice
      common /cembm_r/ tai_laice, tai_co2ccn, tai_dc14ccn, tai_cfc11ccn
      common /cembm_r/ tai_cfc12ccn, tai_maxit, tai_sst, tai_sss
      common /cembm_r/ tai_nsat, tai_ssat, tai_nshum, tai_sshum
      common /cembm_r/ tai_nprecip,  tai_sprecip, tai_nevap, tai_sevap
      common /cembm_r/ tai_nohice, tai_sohice, tai_noaice, tai_soaice
      common /cembm_r/ tai_nhsno, tai_shsno, tai_nlhice, tai_slhice
      common /cembm_r/ tai_nlaice, tai_slaice
