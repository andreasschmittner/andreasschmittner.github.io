      program wind_stress_nc

!=======================================================================
!     creates wind data files wind_stress_mth.nc and wind_stress_ann.nc
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: avar(:,:,:), var(:,:,:,:), vmask(:,:)
      real, allocatable :: xu(:), yu(:)

      integer id, jd
      parameter (id=192, jd=94)
      real data(id,jd,2), xd(id), yd(jd)
      integer maskd(id,jd)

      integer i, imt, iou, j, jmt, n, id_xu, id_yu, id_time
      real psi, theta, phi, daymon(12)

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
      path = trim(path)//'/ncep/surface_gauss/'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xu', iou, imt)
      call getdimlen ('yu', iou, jmt)
      allocate ( xu(imt) )
      allocate ( yu(jmt) )
      allocate ( var(imt,jmt,2,12) )
      allocate ( vmask(imt,jmt) )
      allocate ( maskv(imt,jmt) )
      allocate ( avar(imt,jmt,2) )
      call getvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call getvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read kmt data
!=======================================================================

      call openfile ("../kmt.nc", iou)
      call getvara ('kmt', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &,  vmask, 1., 0.)
      do j=1,jmt-1
        do i=1,imt-1
         maskv(i,j) = int(min(vmask(i,j),vmask(i+1,j)
     &,               vmask(i,j+1),vmask(i+1,j+1)))
        enddo
        maskv(imt,j) = maskv(2,j)
      enddo
      maskv(:,jmt) = 0.
      call closefile (iou)

!=======================================================================
!     read data mask
!=======================================================================

      call openfile (trim(path)//'lsmask.19294.nc', iou)
      call getvara ('lat', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!     flip latitude data
      yd(1:jd) = yd(jd:1:-1)
      call getvara ('lon', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('lsmask', iou, id*jd, (/1,1,1/), (/id,jd,1/)
     &,  data(:,:,1), 1., 0.)
!     flip dmask in latitude and make into ocean mask (0=land, 1=ocean)
      maskd(:,1:jd) = int(data(:,jd:1:-1,1) + 1)
      call closefile (iou)

      avar(:,:,:) = 0.0
      do n=1,12
      
!=======================================================================
!       read wind stress data
!=======================================================================

        call openfile (trim(path)//'uflx.mon.ltm.nc', iou)
        call getvara ('uflx', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,1), 1., 0.)
!       flip data in latitude and change sign for stress into ocean
        data(:,1:jd,1) = -data(:,jd:1:-1,1)
        call closefile (iou)

        call openfile (trim(path)//'vflx.mon.ltm.nc', iou)
        call getvara ('vflx', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,2), 1., 0.)
!       flip data in latitude and change sign for stress into ocean
        data(:,1:jd,2) = -data(:,jd:1:-1,2)
        call closefile (iou)

!=======================================================================
!       mask, rotate and interpolate data
!=======================================================================

        call set_land (data(:,:,1), 2.e20, maskd, id, jd, 1)
        call set_land (data(:,:,2), 2.e20, maskd, id, jd, 1)
        call rot_intrp_vctr (data(:,:,:), xd, yd, id, jd, var(:,:,:,n)
     &,   xu, yu, imt, jmt, phi, theta, psi, -1.e20, 0)
        call extrap (var(:,:,1,n), -1.e10, xu, maskv, imt, jmt, 1)
        call extrap (var(:,:,2,n), -1.e10, xu, maskv, imt, jmt, 1)
        call set_land (var(:,:,1,n), 0., maskv, imt, jmt, 1)
        call set_land (var(:,:,2,n), 0., maskv, imt, jmt, 1)

!=======================================================================
!       set cyclic boundary condition
!=======================================================================

        var(1,:,:,n) = var(imt-1,:,:,n)
        var(imt,:,:,n) = var(2,:,:,n)

        avar(:,:,:) = avar(:,:,:) + var(:,:,:,n)
      enddo
      avar(:,:,:) = avar(:,:,:)/12.

!      vmask(:,:) = 0.
!      where (maskv(:,:) .gt. 0) vmask(:,:) = 1.
       vmask(:,:) = 1.

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call opennew ("../wind_stress_mth.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defdim ('xu', iou, imt, id_xu)
      call defdim ('yu', iou, jmt, id_yu)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call defvar ('xu', iou, 1, (/id_xu/), 0., 0., 'X', 'F'
     &, 'longitude of the u grid', 'longitude', 'degrees_east')
      call defvar ('yu', iou, 1, (/id_yu/), 0., 0., 'Y', 'F'
     &, 'latitude of the u grid', 'latitude', 'degrees_north')
      call defvar ('taux', iou, 3, (/id_xu,id_yu,id_time/), -1.e4, 1.e4,
     &  ' ', 'F', 'surface eastward momentum flux'
     &, 'surface_downward_eastward_stress', 'Pa')
      call defvar ('tauy', iou, 3, (/id_xu,id_yu,id_time/), -1.e4, 1.e4,
     &  ' ', 'F', 'surface northward momentum flux'
     &, 'surface_downward_northward_stress', 'Pa')
      call enddef (iou)
      call putvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call putvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        call putvaramsk ('taux', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(1,1,1,n), vmask, 1., 0.)
        call putvaramsk ('tauy', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(:,:,2,n), vmask, 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../wind_stress_ann.nc", iou)
      call redef (iou)
      call defdim ('xu', iou, imt, id_xu)
      call defdim ('yu', iou, jmt, id_yu)
      call defvar ('xu', iou, 1, (/id_xu/), 0., 0., 'X', 'F'
     &, 'longitude of the u grid', 'longitude', 'degrees_east')
      call defvar ('yu', iou, 1, (/id_yu/), 0., 0., 'Y', 'F'
     &, 'latitude of the u grid', 'latitude', 'degrees_north')
      call defvar ('taux', iou, 2, (/id_xu,id_yu/), -1.e4, 1.e4,
     &  ' ', 'F', 'surface eastward momentum flux'
     &, 'surface_downward_eastward_stress', 'Pa')
      call defvar ('tauy', iou, 2, (/id_xu,id_yu/), -1.e4, 1.e4,
     &  ' ', 'F', 'surface northward momentum flux'
     &, 'surface_downward_northward_stress', 'Pa')
      call enddef (iou)
      call putvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call putvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      call putvaramsk ('taux', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar(:,:,1), vmask, 1., 0.)
      call putvaramsk ('tauy', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar(:,:,2), vmask, 1., 0.)
      call closefile (iou)

      end
