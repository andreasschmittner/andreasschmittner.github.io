      subroutine embm_snap_def (fname, imt, jmt, nat, ncat, xt, yt
     &,                         timunit, expnam, runstamp, mapat)
!=======================================================================
!     definition routine for atmospheric time averages or snapshots

!   inputs:
!     fname        = file name
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     timunit      = time unit
!     expnam       = experiment name
!     runstamp     = run stamp
!     mapat        = tracer map

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

      integer iou, j, n, imt, jmt, nat, ncat, igs, ige, ig, jgs, jge
      integer jg, it(10), iu(10), id_time, id_xt, id_xu, id_yt, id_yu
      integer id_cat, id_cat_e, id_xt_e, id_xu_e, id_yt_e, id_yu_e

      character(*) :: fname, timunit, expnam, runstamp
      character(3) :: a3
      character(10) :: mapat(nat)

      real xt(imt), yt(jmt)
      real c0, c1, c100, c400, c1e3, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c400 = 400.
      c1e3 = 1.e3
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global attributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)
      call putatttext (iou, 'global', 'time_unit', timunit)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      call defdim ('xt', iou, ig, id_xt)
      call defdim ('yt', iou, jg, id_yt)
      call defdim ('xu', iou, ig, id_xu)
      call defdim ('yu', iou, jg, id_yu)
      call defdim ('xt_edges', iou, ig+1, id_xt_e)
      call defdim ('yt_edges', iou, jg+1, id_yt_e)
      call defdim ('xu_edges', iou, ig+1, id_xu_e)
      call defdim ('yu_edges', iou, jg+1, id_yu_e)
#  if defined uvic_ice_cpts || defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1
      call defdim ('cat', iou, ncat, id_cat)
      call defdim ('cat_edges', iou, ncat+1, id_cat_e)
#  endif

!-----------------------------------------------------------------------
!       define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &, 'time since initial condition', 'time', trim(timunit))
      call defvar ('year', iou, 1, it, c0, c0, ' ', 'F'
     &, 'year', ' ',' ')
      call defvar ('month', iou, 1, it, c0, c0, ' ', 'F'
     &, 'month', ' ',' ')
      call defvar ('day', iou, 1, it, c0, c0, ' ', 'F'
     &, 'day', ' ',' ')
      call defvar ('hour', iou, 1, it, c0, c0, ' ', 'F'
     &, 'hour', ' ',' ')
      call defvar ('minute', iou, 1, it, c0, c0, ' ', 'F'
     &, 'minute', ' ',' ')
      call defvar ('second', iou, 1, it, c0, c0, ' ', 'F'
     &, 'second', ' ',' ')
      call defvar ('period', iou, 1, it, c0, c0, ' ', 'F'
     &, 'averaging period', ' ','days')

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('xt', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude of the t grid', 'grid_longitude', 'degrees_east')
        it(1) = id_yt
      call defvar ('yt', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude of the t grid', 'grid_latitude', 'degrees_north')
        it(1) = id_xu
      call defvar ('xu', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude of the u grid', 'grid_longitude', 'degrees_east')
        it(1) = id_yu
      call defvar ('yu', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude of the u grid', 'grid_latitude', 'degrees_north')
        it(1) = id_xt_e
      call defvar ('xt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude of t grid edges', ' ', 'degrees')
        it(1) = id_yt_e
      call defvar ('yt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude of t grid edges', ' ', 'degrees')
        it(1) = id_xu_e
      call defvar ('xu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude of u grid edges', ' ', 'degrees')
        it(1) = id_yu_e
      call defvar ('yu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude of u grid edges', ' ', 'degrees')
#  if defined uvic_ice_cpts || defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1
        it(1) = id_cat
      call defvar ('cat', iou, 1, it, c0, c0, ' ', 'F'
     &, 'ice category', ' ',' ')
        it(1) = id_cat_e
      call defvar ('cat_edges', iou, 1, it, c0, c0, ' ', 'F'
     &, 'ice category', ' ',' ')
#  endif

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      call defvar ('elev', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'land elevation and ocean depth', 'surface_altitude', 'm')
      call defvar ('tmsk', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'ocean mask', ' ' ,'1')
      call defvar ('mskhr', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'horizontal region mask', ' ' ,'1')
      call defvar ('nriv', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'river mask', ' ' ,'1')
      call defvar ('tlat', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid latitude', 'latitude', 'degrees')
      call defvar ('tlon', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid longitude', 'longitude', 'degrees')
      call defvar ('ulat', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid latitude', 'latitude', 'degrees')
      call defvar ('ulon', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid longitude', 'longitude', 'degrees')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          call defvar ('slat', iou, 3, it, c100, c400, ' ', 'F'
     &,     'sea level atmospheric temperature'
     &,     'air_temperature', 'K')
        elseif (trim(mapat(n)) .eq. 'shum') then
          call defvar ('shum', iou, 3, it, -c100, c100, ' ', 'F'
     &,     'atmospheric surface specific humidity'
     &,     'specific_humidity', '1')
        elseif (trim(mapat(n)) .eq. 'co2') then
          call defvar ('co2', iou, 3, it, c0, c1e6, ' ', 'F'
     &,     'atmospheric co2', 'atmospheric_co2', 'ppm')
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call defvar ('at_'//trim(a3), iou ,3, it, -c1e6, c1e6
     &,     'F', 'unknown tracer '//trim(a3)
     &,     ' ', 'unknown_tracer_'//trim(a3), 'unknown')
        endif
      enddo
      call defvar ('sat', iou, 3, it, c100, c400, ' ', 'F'
     &, 'atmospheric surface temperature', 'air_temperature', 'K')
      call defvar ('precip', iou, 3, it, c0, c1, ' ', 'F'
     &, 'precipitation (includes snow in water equivalent)'
     &, 'precipitation_flux', 'kg m-2 s-1')
      call defvar ('psno', iou, 3, it, c0, c1, ' ', 'F'
     &, 'precipitation as snow'
     &, 'snowfall_flux', 'kg m-2 s-1')
      call defvar ('evap', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'upward evaporation plus sublimation'
     &, 'water_evaporation_flux', 'kg m-2 s-1')
      call defvar ('disch', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'river discharge', 'river_discharge_flux', 'kg m-2 s-1')
      call defvar ('outlwr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'outgoing longwave from top of atmosphere'
     &, 'toa_outgoing_longwave_flux', 'W m-2')
      call defvar ('uplwr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward longwave'
     &, 'surface_net_upward_longwave_flux', 'W m-2')
      call defvar ('upsens', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward sensible heat'
     &, 'surface_upward_sensible_heat_flux', 'W m-2')
      call defvar ('dnswr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'net surface downward short wave (absorbed)'
     &, 'surface_net_downward_shortwave_flux', 'W/m^2')
      call defvar ('upltnt', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward latent heat'
     &, 'surface_upward_latent_heat_flux', 'W m-2')
      call defvar ('solins', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'incoming solar insolation'
     &, 'toa_incoming_shortwave_flux', 'W m-2')
      call defvar ('p_alb', iou, 3, it, c0, c1, ' ', 'F'
     &, 'planetary albedo', 'planetary_albedo', '1')
      call defvar ('a_alb', iou, 3, it, c0, c1, ' ', 'F'
     &, 'atmospheric albedo', ' ', '1')
      call defvar ('s_alb', iou, 3, it, c0, c1, ' ', 'F'
     &, 'surface albedo', 'surface_albedo', '1')
      call defvar ('ws', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'surface wind speed', 'wind_speed', 'm s-1')
      call defvar ('runoff', iou, 3, it, c0, c1, ' ', 'F'
     &, 'runoff', 'runoff_flux', 'kg m-2 s-1')
# if defined uvic_embm_adv_q
      call defvar ('wx_q', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward wind for advection of humidity'
     &, 'eastward_wind', 'm s-1')
      call defvar ('wy_q', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward wind for advection of humidity'
     &, 'northward_wind', 'm s-1')
# endif
# if defined uvic_embm_adv_t
      call defvar ('wx_t', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward wind for advection of temperature'
     &, 'eastward_wind', 'm s-1')
      call defvar ('wy_t', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward wind for advection of temperature'
     &, 'northward_wind', 'm s-1')
# endif
# if defined uvic_embm_awind
      call defvar ('awx', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward wind  anomaly', 'eastward_wind_anomaly', 'm s-1')
      call defvar ('awy', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward wind anomaly', 'northward_wind_anomaly', 'm s-1')
      call defvar ('rtbar', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'running average sea level atmospheric temperature'
     &, 'air_temperature', 'K')
      call defvar ('apress', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'anomalous surface pressure', 'pressure', 'Pa')
# endif
      call defvar ('soilm', iou, 3, it, c0, c1e3, ' ', 'F'
     &, 'soil moisture', 'soil_moisture_content', 'kg m-2')
      call defvar ('surf', iou, 3, it, c100, c400, ' ', 'F'
     &, 'land surface temperature', 'surface_temperature', 'K')
# if defined uvic_save_flxadj
      call defvar ('flxadj_t', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'surface downward heat flux adjustment'
     &, ' ', 'W m-2')
      call defvar ('flxadj_s', iou,3, it, -c100, c100, ' ', 'F'
     &, 'surface downward salt flux adjustment'
     &, ' ', 'kg m-2 s-1')
# endif
# if defined uvic_embm_save_diff
      do n=1,nat
        if (n .lt. 1000) write(a3,'(i3)') n
        if (n .lt. 100) write(a3,'(i2)') n
        if (n .lt. 10) write(a3,'(i1)') n
        call defvar ('dn_'//trim(a3), iou ,3, it, -c1e20, c1e20, ' '
     &,   'F', 'northward diffusivity tracer '//trim(a3)
     &,   'northward diffusivity tracer_'//trim(a3), 'unknown')
        call defvar ('de_'//trim(a3), iou ,3, it, -c1e20, c1e20, ' '
     &,   'F', 'eastward diffusivity tracer '//trim(a3)
     &,   'eastward diffusivity tracer_'//trim(a3), 'unknown')
      enddo
# endif
# if defined uvic_embm_icedata
      call defvar ('aicel', iou, 3, it, c0, c100, ' ', 'F'
     &, 'ice sheet area fraction', 'land_ice_area_fraction', '1')
      call defvar ('hicel', iou, 3, it, c0, c1e6, ' ', 'F'
     &, 'ice sheet thickness', ' ', 'm')
# endif
# if defined uvic_carbon && defined uvic_carbon_co2_2d
      call defvar ('flxco2', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'surface downward flux of co2', ' ', 'ppmv s-1')
# endif
# if defined uvic_ice && !defined uvic_ice_cpts && !defined uvic_ice_cpts10 && !defined uvic_ice_cpts5 && !defined uvic_ice_cpts3 && !defined uvic_ice_cpts1
      call defvar ('tice', iou, 3, it, c100, c400, ' ', 'F'
     &, 'surface ice temperature', 'surface_temperature'
     &, 'K')
      call defvar ('hice', iou, 3, it, c0, c1e6, ' ', 'F'
     &, 'ice thickness', 'sea_ice_thickness', 'm')
      call defvar ('aice', iou, 3, it, c0, c100, ' ', 'F'
     &, 'ice area fraction (includes land ice area fraction)'
     &, 'sea_ice_area_fraction', '1')
      call defvar ('hsno', iou, 3, it, -0.01, c1e6, ' ', 'F'
     &, 'surface snow thickness', 'surface_snow_thickness', 'm')
# endif
# if defined uvic_ice_evp
      call defvar ('uice', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward ice velocity', 'eastward_sea_ice_velocity'
     &, 'm s-1')
      call defvar ('vice', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward ice velocity', 'northward_sea_ice_velocity'
     &, 'm s-1')
      call defvar ('xint', iou, 3, iu, -c1e20, c1e20, ' ', 'F'
     &, 'eastward ice interaction'
     &, 'downward_eastward_stress_at_sea_ice_base', 'Pa')
      call defvar ('yint', iou, 3, iu, -c1e20, c1e20, ' ', 'F'
     &, 'northward ice interaction'
     &, 'downward_northward_stress_at_sea_ice_base', 'Pa')
# endif
# if defined uvic_ice_cpts || defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1

!-----------------------------------------------------------------------
!       define time dependent 4d data (x,y,c,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_cat
      iu(3) = id_cat
      it(4) = id_time
      iu(4) = id_time
      call defvar ('ticen', iou, 4, it, c100, c400, ' ', 'F'
     &, 'multi-category ice temperature'
     &, 'surface_sea_ice_temperature', 'K')
      call defvar ('hicen', iou, 4, it, c0, c1e6, ' ', 'F'
     &, 'multi-category ice thickness', 'sea_ice_thickness', 'm')
      call defvar ('aicen', iou, 4, it, c0, c100, ' ', 'F'
     &, 'multi-category ice area', 'sea_ice_area_fraction', '1')
      call defvar ('hsnon', iou, 4, it, c0, c1e6, ' ', 'F'
     &, 'multi-category snow thickness', 'surface_snow_thickness'
     &, 'm')
# endif

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)
#endif

      return
      end

      subroutine embm_snap_out (fname, ids, ide, jds, jde, imt, jmt, nat
     &,                         ncat, xt, yt, xu, yu, dxt, dyt, dxu, dyu
     &,                         avgper, relyr, stamp, mapat, at, sat
     &,                         precip, evap, disch, outlwr, uplwr
     &,                         upsens, dnswr, upltnt, solins, p_alb
     &,                         a_alb, s_alb, elev, psno, ws, runoff
#if defined uvic_embm_adv_q
     &,                         wx_q, wy_q
#endif
#if defined uvic_embm_adv_t
     &,                         wx_t, wy_t
#endif
#if defined uvic_embm_awind
     &,                         awx, awy, rtbar, apress
#endif
     &,                         soilm, surf
#if defined uvic_ice_cpts || defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1
     &,                         ticen, hicen, aicen, hsnon
#elif defined uvic_ice
     &,                         tice, hice, aice, hsno
#endif
#if defined uvic_ice_evp
     &,                         uice, vice, xint, yint
#endif
     &,                         tlat, tlon, ulat, ulon
#if defined uvic_save_flxadj
     &,                         flxadj_t, flxadj_s
#endif
#if defined uvic_embm_save_diff
     &,                         dn, de
#endif
#if defined uvic_embm_icedata
     &,                         aicel, hicel
#endif
#if defined uvic_carbon && defined uvic_carbon_co2_2d
     &,                         flxco2
#endif
     &,                         tmsk, mskhr, nriv, ntrec)
!=======================================================================
!     output routine for atmospheric time averages or snapshots

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!     imt, jmt ...       = global array dimensions
!     xt, yt ...         = global axes
!     dxt, dyt ...       = grid widths
!     avgper             = length of averaging period
!     relyr              = time in years
!     stamp              = time stamp
!     at, ...            = data to be written

!   outputs:
!     ntrec              = number of time record in file

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

      integer iou, j, k, ln, n, ntrec, imt, jmt, nat, ncat
      integer ids, ide, jds, jde, igs, ige, ig, jgs, jge, jg
      integer ils, ile, jls, jle, ib(10), ic(10)
      integer mskhr(ids:ide,jds:jde), nriv(ids:ide,jds:jde)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      character(*) :: fname, stamp
      character(3) :: a3
      character(10) :: mapat(nat)

      real xt(imt), xu(imt), yt(jmt), yu(jmt)
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt)
      real avgper, at(ids:ide,jds:jde,nat)
      real sat(ids:ide,jds:jde), precip(ids:ide,jds:jde)
      real evap(ids:ide,jds:jde), disch(ids:ide,jds:jde)
      real outlwr(ids:ide,jds:jde), uplwr(ids:ide,jds:jde)
      real upsens(ids:ide,jds:jde), dnswr(ids:ide,jds:jde)
      real upltnt(ids:ide,jds:jde), solins(ids:ide,jds:jde)
      real p_alb(ids:ide,jds:jde), a_alb(ids:ide,jds:jde)
      real s_alb(ids:ide,jds:jde), elev(ids:ide,jds:jde)
      real psno(ids:ide,jds:jde), ws(ids:ide,jds:jde)
      real runoff(ids:ide,jds:jde)
# if defined uvic_embm_adv_q
      real wx_q(ids:ide,jds:jde), wy_q(ids:ide,jds:jde)
# endif
# if defined uvic_embm_adv_t
      real wx_t(ids:ide,jds:jde), wy_t(ids:ide,jds:jde)
# endif
# if defined uvic_embm_awind
      real awx(ids:ide,jds:jde), awy(ids:ide,jds:jde)
      real rtbar(ids:ide,jds:jde), apress(ids:ide,jds:jde)
# endif
      real soilm(ids:ide,jds:jde), surf(ids:ide,jds:jde)
# if defined uvic_ice_cpts || defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1
      real cat(ncat), cat_e(ncat+1)
      real ticen(ids:ide,jds:jde,ncat), hicen(ids:ide,jds:jde,ncat)
      real aicen(ids:ide,jds:jde,ncat), hsnon(ids:ide,jds:jde,ncat)
# elif defined uvic_ice
      real tice(ids:ide,jds:jde), hice(ids:ide,jds:jde)
      real aice(ids:ide,jds:jde), hsno(ids:ide,jds:jde)
# endif
# if defined uvic_ice_evp
      real uice(ids:ide,jds:jde), vice(ids:ide,jds:jde)
      real xint(ids:ide,jds:jde), yint(ids:ide,jds:jde)
# endif
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real ulat(ids:ide,jds:jde), ulon(ids:ide,jds:jde)
# if defined uvic_save_flxadj
      real flxadj_t(ids:ide,jds:jde), flxadj_s(ids:ide,jds:jde)
# endif
# if defined uvic_embm_save_diff
      real dn(ids:ide,jds:jde,nat), de(ids:ide,jds:jde,nat)
# endif
#if defined uvic_embm_icedata
      real aicel(ids:ide,jds:jde), hicel(ids:ide,jds:jde)
#endif
# if defined uvic_carbon && defined uvic_carbon_co2_2d
      real flxco2(ids:ide,jds:jde)
# endif
      real tmsk(ids:ide,jds:jde)
      real relyr, tmp, xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1)
      real yu_e(jmt+1), tmpmask(ids:ide,jds:jde)
      real c0, c1, c10, c100, c1e3, C2K, p1
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      c0 = 0.
      c1 = 1.
      c10 = 10.
      c100 = 100.
      c1e3 = 1.e3
      C2K = 273.15
      p1 = 0.1

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, relyr, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpijm(ils:ile,jls:jle) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      call putvars ('time', iou, ntrec, relyr, c1, c0)
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = nyear
      call putvars ('year', iou, ntrec, tmp, c1, c0)
      tmp = nmonth
      call putvars ('month', iou, ntrec, tmp, c1, c0)
      tmp = nday
      call putvars ('day', iou, ntrec, tmp, c1, c0)
      tmp = nhour
      call putvars ('hour', iou, ntrec, tmp, c1, c0)
      tmp = nmin
      call putvars ('minute', iou, ntrec, tmp, c1, c0)
      tmp = nsec
      call putvars ('second', iou, ntrec, tmp, c1, c0)
      call putvars ('period', iou, ntrec, avgper, c1, c0)

      if (ntrec .eq. 1) then

!-----------------------------------------------------------------------
!       write 1d data (x, y or z)
!-----------------------------------------------------------------------
        allocate ( tmpi(igs:ige) )
        allocate ( tmpj(jgs:jge) )
        allocate ( tmpie(igs:ige+1) )
        allocate ( tmpje(jgs:jge+1) )

        ib(1) = 1
        ic(1) = ig
        tmpi(igs:ige) = xt(igs:ige)
        call putvara ('xt', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = xu(igs:ige)
        call putvara ('xu', iou, ig, ib, ic, tmpi, c1, c0)

        ic(1) = jg
        tmpj(jgs:jge) = yt(jgs:jge)
        call putvara ('yt', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = yu(jgs:jge)
        call putvara ('yu', iou, jg, ib, ic, tmpj, c1, c0)

        ic(1) = ig + 1
        call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xt_e(igs:ige+1)
        call putvara ('xt_edges', iou, ig+1, ib, ic, tmpie, c1, c0)
        call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xu_e(igs:ige+1)
        call putvara ('xu_edges', iou, ig+1, ib, ic, tmpie, c1, c0)

        ic(1) = jg + 1
        call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
        call putvara ('yt_edges', iou, jg+1, ib, ic, tmpje, c1, c0)
        call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
        call putvara ('yu_edges', iou, jg+1, ib, ic, tmpje, c1, c0)

        deallocate ( tmpi )
        deallocate ( tmpj )
        deallocate ( tmpie )
        deallocate ( tmpje )
#  if defined uvic_ice_cpts || defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1

        do n=1,ncat
          cat(n) = float(n)
          cat_e(n) = float(n) - 0.5
        enddo
        cat_e(ncat+1) = cat(ncat) + 0.5
        ic(1) = ncat
        call putvara ('cat', iou, ncat, ib, ic, cat, 1., 0.)
        ic(1) = ncat+1
        call putvara ('cat_edges', iou, ncat+1, ib, ic, cat_e, 1., 0.)
#  endif

!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
        tmpij(ils:ile,jls:jle) = elev(ils:ile,jls:jle)
        call putvara ('elev', iou, ln, ib, ic, tmpij, c100, c0)
        tmpij(ils:ile,jls:jle) = tmsk(ils:ile,jls:jle)
        call putvara ('tmsk', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = mskhr(ils:ile,jls:jle)
        call putvara ('mskhr', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = nriv(ils:ile,jls:jle)
        call putvara ('nriv', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlat(ils:ile,jls:jle)
        call putvara ('tlat', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlon(ils:ile,jls:jle)
        call putvara ('tlon', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulat(ils:ile,jls:jle)
        call putvara ('ulat', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulon(ils:ile,jls:jle)
        call putvara ('ulon', iou, ln, ib, ic, tmpij, c1, c0)

      endif

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('slat', iou, ln, ib, ic, tmpij, c1, -C2K)
        elseif (trim(mapat(n)) .eq. 'shum') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('shum', iou, ln, ib, ic, tmpij, c1, c0)
        elseif (trim(mapat(n)) .eq. 'co2') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('co2', iou, ln, ib, ic, tmpij, c1, c0)
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('at_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        endif
      enddo
      tmpij(ils:ile,jls:jle) = sat(ils:ile,jls:jle)
      call putvara('sat', iou, ln, ib, ic, tmpij, c1, -C2K)
      tmpij(ils:ile,jls:jle) = precip(ils:ile,jls:jle)
      call putvara ('precip', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = psno(ils:ile,jls:jle)
      call putvara ('psno', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = evap(ils:ile,jls:jle)
      call putvara ('evap', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = disch(ils:ile,jls:jle)
      call putvara ('disch', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = outlwr(ils:ile,jls:jle)
      call putvara ('outlwr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = uplwr(ils:ile,jls:jle)
      call putvara ('uplwr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = upsens(ils:ile,jls:jle)
      call putvara ('upsens', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = dnswr(ils:ile,jls:jle)
      call putvara ('dnswr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = upltnt(ils:ile,jls:jle)
      call putvara ('upltnt', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = solins(ils:ile,jls:jle)
      call putvara ('solins', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpmask(:,:) = 1.
      where (p_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = p_alb(ils:ile,jls:jle)
      call putvaramsk ('p_alb', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
      tmpmask(:,:) = 1.
      where (a_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = a_alb(ils:ile,jls:jle)
      call putvaramsk ('a_alb', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
      tmpmask(:,:) = 1.
      where (s_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = s_alb(ils:ile,jls:jle)
      call putvaramsk ('s_alb', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
      tmpij(ils:ile,jls:jle) = ws(ils:ile,jls:jle)
      call putvara ('ws', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = runoff(ils:ile,jls:jle)
      call putvara ('runoff', iou, ln, ib, ic, tmpij, p1, c0)
# if defined uvic_embm_adv_q
      tmpij(ils:ile,jls:jle) = wx_q(ils:ile,jls:jle)
      call putvara ('wx_q', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = wy_q(ils:ile,jls:jle)
      call putvara ('wy_q', iou, ln, ib, ic, tmpij, c100, c0)
# endif
# if defined uvic_embm_adv_t
      tmpij(ils:ile,jls:jle) = wx_t(ils:ile,jls:jle)
      call putvara ('wx_t', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = wy_t(ils:ile,jls:jle)
      call putvara ('wy_t', iou, ln, ib, ic, tmpij, c100, c0)
# endif
# if defined uvic_embm_awind
      tmpij(ils:ile,jls:jle) = awx(ils:ile,jls:jle)
      call putvara ('awx', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) =  awy(ils:ile,jls:jle)
      call putvara ('awy', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = rtbar(ils:ile,jls:jle)
      call putvara ('rtbar', iou, ln, ib, ic, tmpij, c1, -C2K)
      tmpij(ils:ile,jls:jle) =  apress(ils:ile,jls:jle)
      call putvara ('apress', iou, ln, ib, ic, tmpij, c10, c0)
# endif
      tmpmask(:,:) = 0.
      where (tmsk(:,:) .lt. 0.5) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle)
      call putvaramsk ('soilm', iou, ln, ib, ic, tmpij, tmpijm, p1, c0)
      tmpij(ils:ile,jls:jle) = surf(ils:ile,jls:jle)
      call putvaramsk ('surf', iou, ln, ib, ic, tmpij, tmpijm, c1, -C2K)
# if defined uvic_ice && !defined uvic_ice_cpts && !defined uvic_ice_cpts10 && !defined uvic_ice_cpts5 && !defined uvic_ice_cpts3 && !defined uvic_ice_cpts1
      tmpmask(:,:) = 0.
      where (hice(:,:) .gt. 0) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = tice(ils:ile,jls:jle)
      call putvaramsk ('tice', iou, ln, ib, ic, tmpij, tmpijm, c1, -C2K)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle)
      call putvara ('hice', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = aice(ils:ile,jls:jle)
      call putvara ('aice', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsno(ils:ile,jls:jle)
      call putvara ('hsno', iou, ln, ib, ic, tmpij, c100, c0)
# endif
# if defined uvic_ice_evp
      tmpij(ils:ile,jls:jle) = uice(ils:ile,jls:jle)
      call putvara ('uice', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = vice(ils:ile,jls:jle)
      call putvara ('vice', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = xint(ils:ile,jls:jle)
      call putvara ('xint', iou, ln, ib, ic, tmpij, c10, c0)
      tmpij(ils:ile,jls:jle) = yint(ils:ile,jls:jle)
      call putvara ('yint', iou, ln, ib, ic, tmpij, c10, c0)
# endif
# if defined uvic_save_flxadj
      tmpij(ils:ile,jls:jle) = flxadj_t(ils:ile,jls:jle)
      call putvara('flxadj_t', iou, ln, ib, ic, tmpij, cal2J, c0)
      tmpij(ils:ile,jls:jle) = flxadj_s(ils:ile,jls:jle)
      call putvara('flxadj_s', iou, ln, ib, ic, tmpij, p1, c0)
# endif
# if defined uvic_embm_save_diff
      do n=1,nat
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpij(ils:ile,jls:jle) = dn(ils:ile,jls:jle,n)
        call putvara('dn_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = de(ils:ile,jls:jle,n)
        call putvara('de_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
      enddo
# endif
# if defined uvic_embm_icedata
      tmpij(ils:ile,jls:jle) = aicel(ils:ile,jls:jle)
      call putvara('aicel', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hicel(ils:ile,jls:jle)
      call putvara('hicel', iou, ln, ib, ic, tmpij, c100, c0)
# endif
# if defined uvic_carbon && defined uvic_carbon_co2_2d
      tmpij(ils:ile,jls:jle) = flxco2(ils:ile,jls:jle)
      call putvara('flxco2', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined uvic_ice_cpts || defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1

!-----------------------------------------------------------------------
!     write 4d data (x,y,c,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ic(3) = 1
      ib(4) = ntrec
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpmask(:,:) = 0.
      where (tmsk(:,:) .ge. 0.5) tmpmask(:,:) = 1.
      do k=1,ncat
        ib(3) = k
        tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = ticen(ils:ile,jls:jle,k)
        call putvaramsk ('ticen', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, -C2K)
        tmpij(ils:ile,jls:jle) = hicen(ils:ile,jls:jle,k)
        call putvaramsk ('hicen', iou, ln, ib, ic, tmpij, tmpijm
     &,   c100, c0)
        tmpij(ils:ile,jls:jle) = aicen(ils:ile,jls:jle,k)
        call putvaramsk ('aicen', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
        tmpij(ils:ile,jls:jle) = hsnon(ils:ile,jls:jle,k)
        call putvaramsk ('hsnon', iou, ln, ib, ic, tmpij, tmpijm
     &,   c100, c0)
      enddo
# endif

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)
#endif

      return
      end
