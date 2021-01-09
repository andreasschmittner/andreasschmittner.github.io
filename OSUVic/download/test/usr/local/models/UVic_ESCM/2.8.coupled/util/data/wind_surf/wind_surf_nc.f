      program wind_surf_nc

!=======================================================================
!     creates wind data files wind_surf_mth.nc and wind_surf_ann.nc
!=======================================================================

      implicit none

      real, allocatable :: avar(:,:,:), var(:,:,:,:)
      real, allocatable :: xt(:), yt(:)

      integer id, jd
      parameter (id=144, jd=73)
      real data(id,jd,2), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, k, n, id_xt, id_yt, id_time
      real psi, theta, phi, time, daymon(12)

      logical exists
      
      character(120) :: path

      data daymon(1:6)  / 15.5,  45.,   74.5, 105.,  135.5, 166. /
      data daymon(7:12) /196.5, 227.5, 258.,  288.5, 319.,  349.5/

!=======================================================================
!     set path and read path file if it exists
!=======================================================================

      path = '/usr/local/models/UVic_ESCM/data_source'
      inquire (file='../path', exist=exists)
      if (exists) then
        open (10,file='../path')
        read (10,'(a)') path
      endif
      path = trim(path)//'/ncep/surface/'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate (xt(imt))
      allocate (yt(jmt))
      allocate (var(imt,jmt,3,12))
      allocate (avar(imt,jmt,3))
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read data grid
!=======================================================================

      call openfile (trim(path)//'wspd.mon.ltm.nc', iou)
      call getvara ('lon', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('lat', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!     flip latitude
      yd(1:jd) = yd(jd:1:-1)
      call closefile (iou)


      avar(:,:,1) = 0.
      do n=1,12

!=======================================================================
!     read wind speed data
!=======================================================================

        call openfile (trim(path)//'wspd.mon.ltm.nc', iou)
        call getvara ('wspd', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,1), 1., 0.)
!       flip data in latitude
        data(:,1:jd,1) = data(:,jd:1:-1,1)
        call closefile (iou)
        
!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_sclr (data(:,:,1), xd, yd, id, jd, var(:,:,1,n)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)


!=======================================================================
!       set cyclic boundary condition
!=======================================================================

        var(1,:,1,n) = var(imt-1,:,1,n)
        var(imt,:,1,n) = var(2,:,1,n)

        avar(:,:,1) = avar(:,:,1) + var(:,:,1,n)
      enddo
      avar(:,:,1) = avar(:,:,1)/12.

!=======================================================================
!     write monthly netcdf wind speed data
!=======================================================================

      call opennew ("../wind_surf_mth.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('ws', iou, 3, (/id_xt,id_yt,id_time/), 0., 1.e6
     &, ' ', 'F', 'surface wind speed', 'wind_speed', 'm s-1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        do i=1,imt
          do j=1,jmt
            if ( var(i,j,1,n).lt.0.) print*, 'var: ',i,j,n,var(i,j,1,n)
          enddo
        enddo          
        call putvara ('ws', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(:,:,1,n), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf wind speed data
!=======================================================================

      call opennew ("../wind_surf_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('ws', iou, 2, (/id_xt,id_yt/), 0., 1.e6
     &, ' ', 'F', 'surface wind speed', 'wind_speed', 'm s-1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('ws', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar(:,:,1), 1., 0.)
      call closefile (iou)

      do n=1,12

!=======================================================================
!       read wind components
!=======================================================================

        call openfile (trim(path)//'uwnd.mon.ltm.nc', iou)
        call getvara ('uwnd', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,1), 1., 0.)
!       flip data in latitude
        data(:,1:jd,1) = data(:,jd:1:-1,1)
        call closefile (iou)
        call openfile (trim(path)//'vwnd.mon.ltm.nc', iou)
        call getvara ('vwnd', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,2), 1., 0.)
!       flip data in latitude
        data(:,1:jd,2) = data(:,jd:1:-1,2)
        call closefile (iou)

!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_vctr (data(:,:,1:2), xd, yd, id, jd
     &,   var(:,:,1:2,n), xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

      enddo

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:,:,:) = var(imt-1,:,:,:)
      var(imt,:,:,:) = var(2,:,:,:)

!=======================================================================
!     calculate wind angle from components
!=======================================================================

      avar(:,:,:) = 0.
      do n=1,12
        do j=1,jmt
          do i=1,imt
            avar(i,j,1) = avar(i,j,1) + var(i,j,1,n) 
            avar(i,j,2) = avar(i,j,2) + var(i,j,2,n) 
            var(i,j,3,n) = sqrt(var(i,j,1,n)**2+(var(i,j,2,n)**2))
            var(i,j,3,n) = min(1.,max(-1.,var(i,j,1,n)/var(i,j,3,n)))
            var(i,j,3,n) = acos(var(i,j,3,n))
            if (var(i,j,2,n) .lt. 0.) var(i,j,3,n) = -var(i,j,3,n)
          enddo
        enddo
      enddo
      avar(:,:,1:2) = avar(:,:,1:2)/12.
      do j=1,jmt
        do i=1,imt
          avar(i,j,3) = sqrt(avar(i,j,1)**2+(avar(i,j,2)**2))
          avar(i,j,3) = min(1.,max(-1.,avar(i,j,1)/avar(i,j,3)))
          avar(i,j,3) = acos(avar(i,j,3))
          if (avar(i,j,2) .lt. 0.) avar(i,j,3) = -avar(i,j,3)
        enddo
      enddo

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call openfile ("../wind_surf_mth.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('wa', iou, 3, (/id_xt,id_yt,id_time/), -4.
     &,  4., ' ', 'F', 'wind angle referenced to the east'
     &, ' ', 'radians')
      call enddef (iou)
      do n=1,12
        call putvara ('wa', iou, imt*jmt, (/1,1,n/)
     &,              (/imt,jmt,1/), var(:,:,3,n), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call openfile ("../wind_surf_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('wa', iou, 2, (/id_xt,id_yt/), -4.
     &,  4., ' ', 'F', 'wind angle referenced to the east'
     &, ' ', 'radians')
      call enddef (iou)
      call putvara ('wa', iou, imt*jmt, (/1,1/)
     &,              (/imt,jmt/), avar(:,:,3), 1., 0.)
      call closefile (iou)

      end
