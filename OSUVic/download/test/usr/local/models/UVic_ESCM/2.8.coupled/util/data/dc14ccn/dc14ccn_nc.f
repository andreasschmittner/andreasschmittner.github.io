      program dc14ccn_nc

!=======================================================================
!     creates ice data file dc14ccn.nc
!=======================================================================

      implicit none

      integer id_time, iou, n, nr, m
      parameter (nr=233)
      real time(nr), dc14ccnn(nr), dc14ccne(nr), dc14ccns(nr)

!=======================================================================
!     define netcdf file
!=======================================================================

      open (10,file='dc14ccn.txt')
      do n=1,nr
        m = n
        read (10,'(f7.1, 3f10.3)')
     $    time(m), dc14ccnn(m), dc14ccne(m), dc14ccns(m)
      enddo
      close(10)

      call opennew ("../dc14ccn.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call defvar ('dc14ccnn', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'northern c14 concentration', '', 'permil')
      call defvar ('dc14ccne', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'equatorial c14 concentration', '', 'permil')
      call defvar ('dc14ccns', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'sorthern c14 concentration', '', 'permil')
      call enddef (iou)
      do n=1,nr
        call putvars ('time', iou, n, time(n), 1., 0.)
        call putvars ('dc14ccnn', iou,  n,  dc14ccnn(n), 1., 0.)
        call putvars ('dc14ccne', iou,  n,  dc14ccne(n), 1., 0.)
        call putvars ('dc14ccns', iou,  n,  dc14ccns(n), 1., 0.)
      enddo
      call closefile (iou)

      end
