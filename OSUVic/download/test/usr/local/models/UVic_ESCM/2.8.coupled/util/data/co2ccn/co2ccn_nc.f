      program co2ccn_nc

!=======================================================================
!     creates ice data file co2ccn.nc
!=======================================================================

      implicit none

      integer id_time, iou, n, ntrec, m
      parameter (ntrec=589)
      real time(ntrec), co2(ntrec), data(ntrec) 

!=======================================================================
!     define netcdf file
!=======================================================================

      open (10,file='co2ccn.txt')
      do n=1,ntrec
        m = ntrec + 1 - n
        read (10,'(f9.1,f8.3,f4.1)') time(m), co2(m), data(m)
        if (data(m) .eq. 0) then
          data(m) = 1 
        elseif (data(m) .eq. 1) then
          data(m) = 2 
        elseif (data(m) .eq. 2) then
          data(m) = 3 
        elseif (data(m) .eq. 3) then
          data(m) = 4 
        elseif (data(m) .eq. 4) then
          data(m) = 5 
        endif
      enddo

      call opennew ("../co2ccn.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call defvar ('co2ccn', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'co2 concentration', '', 'ppmv')
      call defvar ('data_set', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'data set', '', '1')
      call putatttext (iou, 'data_set', 'data_set_1'
     &, 'http://cdiac.esd.ornl.gov/ftp/maunaloa-co2/maunaloa.co2')
      call putatttext (iou, 'data_set', 'data_set_2'
     &, 'http://cdiac.ornl.gov/ftp/trends/co2/lawdome.combined.dat')
      call putatttext (iou, 'data_set', 'data_set_3'
     &, 'http://www.ngdc.noaa.gov/paleo/taylor/taylor.html')
      call putatttext (iou, 'data_set', 'data_set_4'
     &, 'http://www.ngdc.noaa.gov/paleo/taylor/taylor-latequat.html')
      call putatttext (iou, 'data_set', 'data_set_5'
     &, 'http://www.ngdc.noaa.gov/paleo/taylor/taylor-glacial.html')
      call putatttext (iou, 'data_set', 'data_set_6'
     &, 'http://cdiac.esd.ornl.gov/trends/co2/vostok.htm')
      call enddef (iou)
      do n=1,ntrec
        call putvars ('time', iou, n, time(n), 1., 0.)
        call putvars ('co2ccn', iou,  n, co2(n), 1., 0.)
        call putvars ('data_set', iou,  n, data(n), 1., 0.)
      enddo
      call closefile (iou)

      end
