      subroutine embm_tsi_def (fname, timunit, expnam, runstamp)
!=======================================================================
!     output routine for atmospheric time step integrals

!   inputs:
!     fname      = file name
!     timunit    = time unit
!     expnam     = experiment name
!     runstamp   = run stamp

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

      character(*) :: fname, timunit, expnam, runstamp

      integer id(1), id_time, iou

      real c0, c1, c100, c400, c1e3, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c400 = 400.
      c1e3 = 1.e3
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

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
      id(1) = id_time

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      call defvar ('time', iou, 1, id, c0, c0, 'T', 'D'
     &,   'time since initial condition', 'time', trim(timunit))
      call defvar ('year', iou, 1, id, c0, c0, ' ', 'F'
     &,   'year', ' ',' ')
      call defvar ('month', iou, 1, id, c0, c0, ' ', 'F'
     &,   'month', ' ',' ')
      call defvar ('day', iou, 1, id, c0, c0, ' ', 'F'
     &,   'day', ' ',' ')
      call defvar ('hour', iou, 1, id, c0, c0, ' ', 'F'
     &,   'hour', ' ',' ')
      call defvar ('minute', iou, 1, id, c0, c0, ' ', 'F'
     &,   'minute', ' ',' ')
      call defvar ('second', iou, 1, id, c0, c0, ' ', 'F'
     &,   'second', ' ',' ')
      call defvar ('period', iou, 1, id, c0, c0, ' ', 'F'
     &,   'averaging period', ' ','days')
      call defvar ('sat', iou, 1, id, c100, c400, ' ', 'F'
     &,   'global average surface air temperature', ' ', 'K')
      call defvar ('shum', iou, 1, id, -c100, c100, ' ', 'F'
     &,   'global average surface specific humidity', ' ', '1')
      call defvar ('precip', iou, 1, id, -c100, c100, ' '
     &,   'F', 'global average precipitation', ' ','kg m-2 s-1')
      call defvar ('evap', iou, 1, id, -c100, c100, ' '
     &,   'F', 'global average evaporation', ' ','kg m-2 s-1')
      call defvar ('co2ccn', iou, 1, id, c0, c1e3, ' '
     &,   'F', 'global average CO2 concentration', ' ','ppm')
      call defvar ('maxit', iou, 1, id, c0, c1e3, ' '
     &,   'F', 'maximum atmosphere solver iterations', ' ','')
      call defvar ('sst', iou, 1, id, -c100, c400, ' '
     &,   'F', 'global average sea surface temperature', ' ','K')
      call defvar ('sss', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface salinity', ' ','psu')
# if defined uvic_carbon_14
      call defvar ('dc14ccn', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &,   'global average delta C14 concentration', ' ','permil')
# endif
#  if defined uvic_cfc || defined uvic_cfc11
      call defvar ('cfc11ccn', iou, 1, id, -c1e3, c1e3, ' '
     &,   'F', 'global average CFC11 concentration', ' ','ppm')
# endif
# if defined uvic_cfc || defined uvic_cfc12
      call defvar ('cfc12ccn', iou, 1, id, -c1e3, c1e3, ' '
     &,   'F', 'global average CFC12 concentration', ' ','ppm')
# endif
# if defined uvic_ice
      call defvar ('snow_volume', iou, 1, id, c0, c1e20, ' '
     &,   'F', 'global snow volume', ' ', 'm3')
      call defvar ('sea_ice_volume', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global sea ice volume', ' ', 'm3')
      call defvar ('sea_ice_area', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global sea ice area', ' ', 'm2')
#   if defined uvic_icedata_transient || defined uvic_ism
      call defvar ('land_ice_volume', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'global land ice volume', ' ', 'm3')
      call defvar ('land_ice_area', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global land ice area', ' ', 'm2')
#  endif
# endif
# if defined uvic_tai_ns
      call defvar ('n_sat', iou, 1, id, c100, c400, ' '
     &,   'F', 'NH average surface air temperature', ' ', 'K')
      call defvar ('s_sat', iou, 1, id, c100, c400, ' '
     &,   'F', 'SH average surface air temperature', ' ', 'K')
      call defvar ('n_shum', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average surface specific humidity', ' ', '1')
      call defvar ('s_shum', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average surface specific humidity', ' ', '1')
      call defvar ('n_precip', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average precipitation', ' ','kg m-2 s-1')
      call defvar ('s_precip', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average precipitation', ' ','kg m-2 s-1')
      call defvar ('n_evap', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average evaporation', ' ','kg m-2 s-1')
      call defvar ('s_evap', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average evaporation', ' ','kg m-2 s-1')
#  if defined uvic_ice
      call defvar ('n_snow_volume', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH snow volume', ' ', 'm3')
      call defvar ('s_snow_volume', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH snow volume', ' ', 'm3')
      call defvar ('n_sea_ice_volume', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH sea ice volume', ' ', 'm3')
      call defvar ('s_sea_ice_volume', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH sea ice volume', ' ', 'm3')
      call defvar ('n_sea_ice_area', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH sea ice area', ' ', 'm2')
      call defvar ('s_sea_ice_area', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH sea ice area', ' ', 'm2')
#   if defined uvic_icedata_transient || defined uvic_ism
      call defvar ('n_land_ice_volume', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'NH land ice volume', ' ', 'm3')
      call defvar ('s_land_ice_volume', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'SH land ice volume', ' ', 'm3')
      call defvar ('n_land_ice_area', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'NH land ice area', ' ', 'm2')
      call defvar ('s_land_ice_area', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'SH land ice area', ' ', 'm2')
#   endif
#  endif
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

      subroutine embm_tsi_out (fname, avgper, relyr, stamp, sat, shum
     &,                        precip, evap, v_oice, a_oice, v_snow
     &,                        v_lice, a_lice, co2ccn, dc14ccn, cfc11ccn
     &,                        cfc12ccn, scan, sst, sss, nsat, ssat
     &,                        nshum, sshum, nprecip, sprecip, nevap
     &,                        sevap, v_noice, v_soice, a_noice, a_soice
     &,                        v_nsnow, v_ssnow, v_nlice, v_slice
     &,                        a_nlice, a_slice, ntrec)
!=======================================================================
!     output routine for atmospheric time step integrals

!   inputs:
!     fname    = file name
!     avgper   = length of averaging period
!     relyr    = time in years
!     stamp    = time stamp
!     sat, ... = data to be written

!   outputs:
!     ntrec    = number of time record in file

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real sat, shum, precip, evap, v_oice, a_oice, v_snow, v_lice
      real a_lice, co2ccn, scan, sst, sss, dc14ccn, cfc11ccn, cfc12ccn
      real nsat, ssat, nshum, sshum, nprecip, sprecip, nevap, sevap
      real v_noice, v_soice, a_noice, a_soice,  v_nsnow, v_ssnow
      real v_nlice, v_slice, a_nlice, a_slice,  avgper, relyr, tmp
      real c0, c1, c100, c1e4, c1e6, C2K, p1, p001, p035, cal2J

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c1e4 = 1.e4
      c1e6 = 1.e6
      C2K = 273.15
      p1 = 0.1
      p001 = 0.001
      p035 = 0.035
      cal2J = 2.389e-05

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, relyr, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

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
      call putvars ('sat', iou, ntrec, sat, c1, -C2K)
      call putvars ('shum', iou, ntrec, shum, c1, c0)
      call putvars ('precip', iou, ntrec, precip, p1, c0)
      call putvars ('evap', iou, ntrec, evap, p1, c0)
      call putvars ('co2ccn', iou, ntrec, co2ccn, c1, c0)
      call putvars ('maxit', iou, ntrec, scan, c1, c0)
      call putvars ('sst', iou, ntrec, sst, c1, -C2K)
      call putvars ('sss', iou, ntrec, sss, p001, -p035)
#  if defined uvic_carbon_14
      call putvars ('dc14ccn', iou, ntrec, dc14ccn, c1, c0)
#  endif
#  if defined uvic_cfc || defined uvic_cfc11
      call putvars ('cfc11ccn', iou, ntrec, cfc11ccn, c1, c0)
#  endif
#  if defined uvic_cfc || defined uvic_cfc12
      call putvars ('cfc12ccn', iou, ntrec, cfc12ccn, c1, c0)
#  endif
#  if defined uvic_ice
      call putvars ('sea_ice_volume', iou, ntrec, v_oice, c1e6, c0)
      call putvars ('sea_ice_area', iou, ntrec, a_oice, c1e4, c0)
      call putvars ('snow_volume', iou, ntrec, v_snow, c1e6, c0)
#   if defined uvic_icedata_transient || defined uvic_ism
      call putvars ('land_ice_volume', iou, ntrec, v_lice, c1e6, c0)
      call putvars ('land_ice_area', iou, ntrec, a_lice, c1e4, c0)
#   endif
#  endif
#  if defined uvic_tai_ns
      call putvars ('n_sat', iou, ntrec, nsat, c1, -C2K)
      call putvars ('s_sat', iou, ntrec, ssat, c1, -C2K)
      call putvars ('n_shum', iou, ntrec, nshum, c1, c0)
      call putvars ('s_shum', iou, ntrec, sshum, c1, c0)
      call putvars ('n_precip', iou, ntrec, nprecip, p1, c0)
      call putvars ('s_precip', iou, ntrec, sprecip, p1, c0)
      call putvars ('n_evap', iou, ntrec, nevap, p1, c0)
      call putvars ('s_evap', iou, ntrec, sevap, p1, c0)
#   if defined uvic_ice
      call putvars ('n_sea_ice_volume', iou, ntrec, v_noice, c1e6, c0)
      call putvars ('s_sea_ice_volume', iou, ntrec, v_soice, c1e6, c0)
      call putvars ('n_sea_ice_area', iou, ntrec, a_noice, c1e4, c0)
      call putvars ('s_sea_ice_area', iou, ntrec, a_soice, c1e4, c0)
      call putvars ('n_snow_volume', iou, ntrec, v_nsnow, c1e6, c0)
      call putvars ('s_snow_volume', iou, ntrec, v_ssnow, c1e6, c0)
#    if defined uvic_icedata_transient || defined uvic_ism
      call putvars ('n_land_ice_volume', iou, ntrec, v_nlice, c1e6, c0)
      call putvars ('s_land_ice_volume', iou, ntrec, v_slice, c1e6, c0)
      call putvars ('n_land_ice_area', iou, ntrec, a_nlice, c1e4, c0)
      call putvars ('s_land_ice_area', iou, ntrec, a_slice, c1e4, c0)
#    endif
#   endif
#  endif

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)
#endif

      return
      end
