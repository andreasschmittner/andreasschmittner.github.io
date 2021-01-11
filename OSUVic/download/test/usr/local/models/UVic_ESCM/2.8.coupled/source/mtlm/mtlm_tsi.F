      subroutine mtlm_tsi_def (fname, timunit, expnam, runstamp)
!=======================================================================
!     output routine for land time step integrals

!   inputs:
!     fname      = file name
!     timunit    = time unit
!     expnam     = experiment name
!     runstamp   = run stamp

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_mtlm

      implicit none

      character(*) :: fname, timunit, expnam, runstamp

      integer id(1), id_time, iou, ntrec

      real c0, c1, c1e20

      c0 = 0.
      c1 = 1.
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
      call defvar ('CS', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon', ' ', 'kg')
      call defvar ('RESP_S', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration', ' ', 'kg s-1')
      call defvar ('CV', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon', ' ', 'kg')
      call defvar ('NPP', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity', ' ', 'kg s-1')
      call defvar ('GPP', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total gross primary productivity', ' ', 'kg s-1')
      call defvar ('HT', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average vegetation height', ' ', 'm')
      call defvar ('LAI', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average leaf area index', ' ', '1')
      call defvar ('LYING_SNOW', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total snow mass', ' ', 'kg')
      call defvar ('TSOIL', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average soil temperature', ' ', 'K')
      call defvar ('TSTAR', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average skin temperature', ' ', 'K')
      call defvar ('M', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average soil moisture', ' ', 'Kg m-2')
      call defvar ('ET', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average evapotranspiration', ' ', 'kg m-2 s-1')

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)
#endif

      return
      end

      subroutine mtlm_tsi_out (fname, avgper, relyr, stamp, CS, RESP_S
     &,                        CV, NPP, GPP, HT, LAI, LYING_SNOW, TSOIL
     &,                        TSTAR, M, ET, ntrec)
!=======================================================================
!     output routine for land time step integrals

!   inputs:
!     fname   = file name
!     avgper  = length of averaging period
!     relyr   = time in years
!     stamp   = time stamp
!     CS, ... = data to be written

!   outputs:
!     ntrec   = number of time record in file

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_mtlm

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real CS, RESP_S, CV, NPP, GPP, HT, LAI, LYING_SNOW, TSOIL
      real TSTAR, M, ET, avgper, relyr, tmp, c0, c1

      c0 = 0.
      c1 = 1.

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
      call putvars ('CS', iou, ntrec, CS, c1, c0)
      call putvars ('RESP_S', iou, ntrec, RESP_S, c1, c0)
      call putvars ('CV', iou, ntrec, CV, c1, c0)
      call putvars ('NPP', iou, ntrec, NPP, c1, c0)
      call putvars ('GPP', iou, ntrec, GPP, c1, c0)
      call putvars ('HT', iou, ntrec, HT, c1, c0)
      call putvars ('LAI', iou, ntrec, LAI, c1, c0)
      call putvars ('LYING_SNOW', iou, ntrec, LYING_SNOW, c1, c0)
      call putvars ('TSOIL', iou, ntrec, TSOIL, c1, c0)
      call putvars ('TSTAR', iou, ntrec, TSTAR, c1, c0)
      call putvars ('M', iou, ntrec, M, c1, c0)
      call putvars ('ET', iou, ntrec, ET, c1, c0)

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)
#endif

      return
      end
