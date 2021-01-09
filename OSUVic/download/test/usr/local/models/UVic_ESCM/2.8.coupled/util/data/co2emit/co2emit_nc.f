      program co2emit_nc

!=======================================================================
!     creates ice data file co2emit.nc
!=======================================================================

      implicit none

      integer id_time, iou, n, ntrec, m
      parameter (ntrec=260)
      real time(ntrec), co2emit(ntrec,3), data(ntrec) 

!=======================================================================
!     define netcdf file
!=======================================================================

      open (10,file='co2emit.txt')
      do n=1,ntrec
        m = n
        read (10,'(f7.2,f7.4,f9.6,f9.6)') time(m), co2emit(m,1)
     &,   co2emit(m,2), co2emit(m,3)
      enddo

      call opennew ("../co2emit.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call defvar ('co2emit_fuel', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'co2 emissions from fossil fuels', '', 'GT yr-1')
      call defvar ('co2emit_land', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'co2 emissions from land changes', '', 'GT yr-1')
      call defvar ('co2emit', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'co2 emissions', '', 'GT yr-1')
      call enddef (iou)
      do n=1,ntrec
        call putvars ('time', iou, n, time(n), 1., 0.)
        call putvars ('co2emit_fuel', iou,  n, co2emit(n,1), 1., 0.)
        call putvars ('co2emit_land', iou,  n, co2emit(n,2), 1., 0.)
        call putvars ('co2emit', iou,  n, co2emit(n,3), 1., 0.)
      enddo
      call closefile (iou)

      end
