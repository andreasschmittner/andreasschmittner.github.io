      subroutine mom_tsi_def (fname, timunit, expnam, runstamp)
!=======================================================================
!     output routine for ocean time step integrals

!   inputs:
!     fname      = file name
!     timunit    = time unit
!     expnam     = experiment name
!     runstamp   = run stamp

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_mom

      implicit none

      character(*) :: fname, timunit, expnam, runstamp

      integer id(1), id_time, iou, ntrec

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
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global atributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)
      call putatttext (iou, 'global', 'time_unit', timunit)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      id = id_time

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      call defvar ('time', iou, 1, id, c0, c0, 'T', 'D'
     &, 'time since initial condition', 'time', trim(timunit))
      call defvar ('year', iou, 1, id, c0, c0, ' ', 'F'
     &, 'year', ' ',' ')
      call defvar ('month', iou, 1, id, c0, c0, ' ', 'F'
     &, 'month', ' ',' ')
      call defvar ('day', iou, 1, id, c0, c0, ' ', 'F'
     &, 'day', ' ',' ')
      call defvar ('hour', iou, 1, id, c0, c0, ' ', 'F'
     &, 'hour', ' ',' ')
      call defvar ('minute', iou, 1, id, c0, c0, ' ', 'F'
     &, 'minute', ' ',' ')
      call defvar ('second', iou, 1, id, c0, c0, ' ', 'F'
     &, 'second', ' ',' ')
      call defvar ('period', iou, 1, id, c0, c0, ' ', 'F'
     &, 'averaging period', ' ','days')
      call defvar ('ke', iou, 1, id, c0, c1e6, ' ', 'F'
     &, 'Kinetic energy per unit volume', ' ', 'J m-3')
      call defvar ('tbar', iou, 1, id, c100, c400, ' ', 'F'
     &, 'Global average T', ' ', 'K')
      call defvar ('sbar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average S', ' ', '1e-3')
      call defvar ('tvar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'Variance of T', ' ', 'K2')
      call defvar ('svar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'Variance of S', ' ', '1e-6')
      call defvar ('abs_chg_t', iou, 1, id, c0, c1e20, ' '
     &, 'F', '|dT/dt|', ' ', 'K s-1')
      call defvar ('abs_chg_s', iou, 1, id, c0, c1e20, ' '
     &, 'F', '|dS/dt|', ' ', '1e-3 s-1')
      call defvar ('scans', iou, 1, id, c0, c1e6, ' ', 'F'
     &, 'elliptic solver iterations', ' ', '1')
      call defvar ('hflx', iou, 1, id, -c1e20, c1e20, ' '
     &,   'F', 'global average ocean heat flux', ' ','W m-2')
      call defvar ('sflx', iou, 1, id, -c1e20, c1e20, ' '
     &,   'F', 'global average ocean salt flux', ' ','kg m-2 s-1')
# if defined uvic_carbon
      call defvar ('dicbar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average dic', ' ', 'mol m-3')
      call defvar ('dicflx', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average carbon flux', ' ', 'mol m-2 s-1')
#  if defined uvic_carbon_14
      call defvar ('c14bar', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'Global average c14', ' ', 'mol m-3')
      call defvar ('dc14bar', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'Global average delta c14', ' ', 'permil')
      call defvar ('c14flx', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'Global average c14 flux', ' ', 'mol m-2 s-1')
#  endif
# endif
# if defined uvic_alk
      call defvar ('alkbar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average alkalinity', ' ', 'mol m-3')
# endif
# if defined uvic_o2
      call defvar ('o2bar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average oxygen', ' ', 'mol m-3')
      call defvar ('o2flx', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average flux of oxygen', ' ', 'mol m-2 s-1')
# endif
# if defined uvic_npzd
      call defvar ('po4bar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average phosphate', ' ', 'mol m-3')
      call defvar ('pbar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average phytoplankton', ' ', 'mol N m-3')
      call defvar ('zbar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average zooplankton', ' ', 'mol N m-3')
      call defvar ('dbar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average detritus', ' ', 'mol N m-3')
# endif
# if defined uvic_nitrogen
      call defvar ('no3bar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average nitrate', ' ', 'mol m-3')
      call defvar ('diazbar', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average diazotrophs', ' ', 'mol N m-3')
# endif
# if defined uvic_cfc || defined uvic_cfc11
      call defvar ('cfc11bar', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'Global average CFC11', ' ', 'mol m-3')
      call defvar ('cfc11flx', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'Global average CFC11 flux', ' ', 'mol m-2 s-1')
# endif
# if defined uvic_cfc || defined uvic_cfc12
      call defvar ('cfc12bar', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'Global average CFC12', ' ', 'mol m-3')
      call defvar ('cfc12flx', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'Global average CFC12 flux', ' ', 'mol m-2 s-1')
# endif
# if defined uvic_tai_otsf
      call defvar ('motsf_max', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'maximum meridional overturning streamfunction', ' '
     &, 'm3 s-1')
      call defvar ('motsf_min', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'minimum meridional overturning streamfunction', ' '
     &, 'm3 s-1')
# endif
# if defined uvic_tai_slh
      call defvar ('slh', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'relative sea level height', ' ', 'm')
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

      subroutine mom_tsi_out (fname, avgper, relyr, stamp, ektot, tbar1
     &,                       tbar2, travar1, travar2, dtabs1, dtabs2
     &,                       scan, hflx, sflx, tbar_dic, dicflx
     &,                       tbar_alk, tbar_o2, o2flx, tbar_po4, tbar_p
     &,                       tbar_z, tbar_d, tbar_no3, tbar_di
     &,                       tbar_c14, tbar_dc14, c14flx, tbar_cfc11
     &,                       cfc11flx, tbar_cfc12, cfc12flx, otmax
     &,                       otmin, slh, ntrec)
!=======================================================================
!     output routine for ocean time step integrals

!   inputs:
!     fname      = file name
!     avgper     = length of averaging period
!     relyr      = time in years
!     stamp      = time stamp
!     ektot, ... = data to be written

!   outputs:
!     ntrec      = number of time record in file

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_mom

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real ektot, tbar1, tbar2, travar1, travar2, dtabs1, dtabs2, scan
      real hflx, sflx, tbar_dic, dicflx, tbar_alk, tbar_o2, o2flx
      real tbar_po4, tbar_p, tbar_z, tbar_d, tbar_no3, tbar_di, tbar_c14
      real tbar_dc14, c14flx, tbar_cfc11, cfc11flx, tbar_cfc12, cfc12flx
      real otmax, otmin, slh, avgper, relyr, tmp, c0, c1, c100, c1e3
      real c1e6, C2K, p001, p035, p1, cal2J

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c1e3 = 1.e3
      c1e6 = 1.e6
      C2K = 273.15
      p001 = 0.001
      p035 = 0.035
      p1 = 0.1
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
      call putvars ('ke', iou, ntrec, ektot, c100, c0)
      call putvars ('tbar', iou, ntrec, tbar1, c1, -C2K)
      call putvars ('sbar', iou, ntrec, tbar2, p001, -p035)
      tmp = 0.0
      if (travar1 .gt. 0.0) tmp = (sqrt(travar1)+273.15)**2
      call putvars ('tvar', iou, ntrec, tmp, c1, c0)
      tmp = 0.0
      if (travar2 .gt. 0.0) tmp = (sqrt(travar2)*1000.+35.)**2
      call putvars ('svar', iou, ntrec, tmp, c1, c0)
      tmp = dtabs1
      call putvars ('abs_chg_t', iou, ntrec, dtabs1, c1, c0)
      tmp = dtabs2
      call putvars ('abs_chg_s', iou, ntrec, dtabs2, p001, c0)
      call putvars ('scans', iou, ntrec, scan, c1, c0)
      call putvars ('hflx', iou, ntrec, hflx, cal2J, c0)
      call putvars ('sflx', iou, ntrec, sflx, p1, c0)
# if defined uvic_carbon
      call putvars ('dicbar', iou, ntrec, tbar_dic, c1, c0)
      call putvars ('dicflx', iou, ntrec, dicflx, c100, c0)
#  if defined uvic_carbon_14
      call putvars ('c14bar', iou, ntrec, tbar_c14, c1, c0)
      call putvars ('dc14bar', iou, ntrec, tbar_dc14, c1, c0)
      call putvars ('c14flx', iou, ntrec, c14flx, c100, c0)
#  endif
# endif
# if defined uvic_alk
      call putvars ('alkbar', iou, ntrec, tbar_alk, c1, c0)
# endif
# if defined uvic_o2
      call putvars ('o2bar', iou, ntrec, tbar_o2, c1, c0)
      call putvars ('o2flx', iou, ntrec, o2flx, c100, c0)
# endif
# if defined uvic_npzd
      call putvars ('po4bar', iou, ntrec, tbar_po4, c1e3, c0)
      call putvars ('pbar', iou, ntrec, tbar_p, c1e3, c0)
      call putvars ('zbar', iou, ntrec, tbar_z, c1e3, c0)
      call putvars ('dbar', iou, ntrec, tbar_d, c1e3, c0)
# endif
# if defined uvic_nitrogen
      call putvars ('no3bar', iou, ntrec, tbar_no3, c1e3, c0)
      call putvars ('diazbar', iou, ntrec, tbar_di, c1e3, c0)
# endif
# if defined uvic_cfc || defined uvic_cfc11
      call putvars ('cfc11bar', iou, ntrec, tbar_cfc11, c1, c0)
      call putvars ('cfc11flx', iou, ntrec, cfc11flx, c100, c0)
# endif
# if defined uvic_cfc || defined uvic_cfc12
      call putvars ('cfc12bar', iou, ntrec, tbar_cfc12, c1, c0)
      call putvars ('cfc12flx', iou, ntrec, cfc12flx, c100, c0)
# endif
# if defined uvic_tai_otsf
      call putvars ('motsf_max', iou, ntrec, otmax, c1e6, c0)
      call putvars ('motsf_min', iou, ntrec, otmin, c1e6, c0)
# endif
# if defined uvic_tai_slh
      call putvars ('slh', iou, ntrec, slh, c100, c0)
# endif

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)
#endif

      return
      end
