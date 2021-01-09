      program solar_nc

!=======================================================================
!     creates ice data file solar.nc
!=======================================================================

      implicit none

      integer id_time, iou, n, ntrec, m
      parameter (ntrec=485)
      real time(ntrec), solar(ntrec), data(ntrec) 

!=======================================================================
!     define netcdf file
!=======================================================================

      open (10,file='solar.txt')
      do n=1,ntrec
        m = ntrec + 1 - n
        read (10,'(f8.2,f10.4,f4.1)') time(n), solar(n), data(n)
      enddo

      call opennew ("../solar.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call defvar ('solar_constant', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'solar constant', '', 'W m-2')
      call defvar ('data_set', iou, 1, (/id_time/), 0., 1.
     &, ' ', 'F', 'data set', '', '1')
      call putatttext (iou, 'data_set', 'data_set_1.1'
     &, 'ftp://ftp.ngdc.noaa.gov/paleo/climate_forcing/')
      call putatttext (iou, 'data_set', 'data_set_1.2'
     &, '      solar_variability/bard_irradiance.txt')
      call putatttext (iou, 'data_set', 'data_set_2.1'
     &, 'ftp://ftp.ngdc.noaa.gov/paleo/contributions_by_author/')
      call putatttext (iou, 'data_set', 'data_set_2.2'
     &, '      lean1995/irradiance_data.txt')
      call enddef (iou)
      do n=1,ntrec
        call putvars ('time', iou, n, time(n), 1., 0.)
        call putvars ('solar_constant', iou,  n, solar(n), 1., 0.)
        call putvars ('data_set', iou,  n, data(n), 1., 0.)
      enddo
      call closefile (iou)

      end
