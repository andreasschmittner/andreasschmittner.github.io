      program cfcccn_nc

!=======================================================================
!     creates ice data file cfc.nc
!=======================================================================

      implicit none

      integer id_time, iou, n, nr, m
      parameter (nr=68)
      real time(nr), cfc11n(nr), cfc12n(nr), cfc11s(nr), cfc12s(nr)

!=======================================================================
!     define netcdf file
!=======================================================================

      open (10,file='cfcccn.txt')
      do n=1,nr
        m = n
        read (10,'(5f8.2)') time(m), cfc11n(m), cfc12n(m), cfc11s(m)
     &,   cfc12s(m)        
      enddo
      close(10)

      call opennew ("../cfcccn.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call defvar ('cfc11ccnn', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'northern cfc11 concentration', '', 'ppt')
      call defvar ('cfc11ccns', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'sorthern cfc11 concentration', '', 'ppt')
      call defvar ('cfc12ccnn', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'northern cfc12 concentration', '', 'ppt')
      call defvar ('cfc12ccns', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'sorthern cfc12 concentration', '', 'ppt')
      call enddef (iou)
      do n=1,nr
        call putvars ('time', iou, n, time(n), 1., 0.)
        call putvars ('cfc11ccnn', iou,  n, cfc11n(n), 1., 0.)
        call putvars ('cfc11ccns', iou,  n, cfc11s(n), 1., 0.)
        call putvars ('cfc12ccnn', iou,  n, cfc12n(n), 1., 0.)
        call putvars ('cfc12ccns', iou,  n, cfc12s(n), 1., 0.)
      enddo
      call closefile (iou)

      end
