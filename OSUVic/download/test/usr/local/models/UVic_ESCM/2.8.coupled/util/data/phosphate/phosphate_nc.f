      program phosphate_nc

!=======================================================================
!     creates netcdf files of annual data and monthly surface data
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:,:), tmp(:,:,:), vmask(:,:,:)
      real, allocatable :: xt(:), yt(:), zt(:)

      integer id, jd, kd, kdm
      parameter (id=360, jd=180, kd=33, kdm=7)
      real data(id,jd,kd), xd(id), yd(jd), zd(kd)

      integer i, imt, iou, j, jmt, k, km, n, id_xt, id_yt
      integer id_zt, id_time, ka, kb, m
      real psi, theta, phi, wta, wtb, vmax, vmin, missing_value

      logical exists
      
      character(120) :: path, name, name2, long_name, stand_name, units

!=======================================================================
!     set names and values
!=======================================================================

      name = "phosphate"
      long_name = "phosphate"
      stand_name = "phosphate"
      vmax = 4.2
      vmin = 0.
      units = "mol m-3"
      missing_value = -99.
      name2 = name

!=======================================================================
!     set path and read path file if it exists
!=======================================================================

      path = '/usr/local/models/UVic_ESCM/data_source'
      inquire (file='../path', exist=exists)
      if (exists) then
        open (10,file='../path')
        read (10,'(a)') path
      endif
      path = trim(path)//'/WOA/'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      call getdimlen ('zt', iou, km)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( zt(km) )
      allocate ( var(imt,jmt,km) )
      allocate ( tmp(imt,jmt,2) )
      allocate ( vmask(imt,jmt,km) )
      allocate ( maskv(imt,jmt) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvara ('zt', iou, km, (/1/), (/km/), zt, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read kmt data
!=======================================================================

      call openfile ("../kmt.nc", iou)
      call getvara ('kmt', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, vmask(:,:,1), 1., 0.)
      maskv(:,:) = nint(vmask(:,:,1))
      call closefile (iou)

!=======================================================================
!     read annual data
!=======================================================================

      call openfile (trim(path)//trim(name)//'_annual.nc', iou)
      call getvara ('X', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('Y', iou, jd, (/1/), (/jd/), yd, 1., 0.)
      call getvara ('Z', iou, kd, (/1/), (/kd/), zd, 1., 0.)
      call getvara (name2, iou, id*jd*kd, (/1,1,1,1/)
     &, (/id,jd,kd,1/), data(:,:,:), 1., 0.)
      where (data(:,:,:) .le. missing_value) data(:,:,:) = 2e20

!=======================================================================
!     rotate and interpolate data for all levels
!=======================================================================

      do k=1,km
        ka = 1
        do m=2,kd
          if (zd(m) .lt. zt(k)) ka = m
        enddo
        kb = min(ka+1,kd)
        ka = max(kb-1,1)
        wtb = min(1.0,max(0.0,(zt(k) - zd(ka))/(zd(kb) - zd(ka))))
        call rot_intrp_sclr (data(:,:,ka), xd, yd, id, jd, tmp(:,:,1)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call rot_intrp_sclr (data(:,:,kb), xd, yd, id, jd, tmp(:,:,2)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call intrp_vert (tmp(:,:,1), tmp(:,:,2), wtb, -1e20
     &,   var(:,:,k), imt, jmt)
      enddo
      call extrap (var(:,:,:), -1.e20, xt, maskv(:,:), imt, jmt, km)

!=======================================================================
!       set cyclic boundary condition
!=======================================================================

      var(1,:,:) = var(imt-1,:,:)
      var(imt,:,:) = var(2,:,:)
      
!=======================================================================
!     create 3D mask
!=======================================================================

      vmask(:,:,:) = 0
      do k=1,km
        where (maskv(:,:) .ge. k) vmask(:,:,k) = 1
      enddo

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../"//trim(name)//"_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defdim ('zt', iou, km, id_zt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('zt', iou, 1, (/id_zt/), 0., 0., 'Z', 'F'
     &, 'depth', 'depth', 'depth')
      call defvar (trim(name), iou, 3, (/id_xt,id_yt,id_zt/)
     &, vmin, vmax,' ', 'F', trim(long_name), trim(stand_name)
     &, trim(units))
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('zt', iou, jmt, (/1/), (/km/), zt, 1., 0.)
      call putvaramsk (trim(name), iou, imt*jmt*km, (/1,1,1/)
     &, (/imt,jmt,km/), var(:,:,:), vmask, 1000., 0.)
      call closefile (iou)
      
      do n=1,12

!=======================================================================
!     read monthly data
!=======================================================================

        call openfile (trim(path)//trim(name)//'_monthly.nc', iou)
        call getvara ('X', iou, id, (/1/), (/id/), xd, 1., 0.)
        call getvara ('Y', iou, jd, (/1/), (/jd/), yd, 1., 0.)
        call getvara (name2, iou, id*jd*kdm, (/1,1,1,n/)
     &,   (/id,jd,kdm,1/), data(:,:,1:kdm), 1., 0.)
        where (data(:,:,:) .le. missing_value) data(:,:,:) = 2e20

!=======================================================================
!     rotate and interpolate monthly model surface data
!=======================================================================

        k = 1
        ka = 1
        do m=2,kd
          if (zd(m) .lt. zt(k)) ka = m
        enddo
        kb = min(ka+1,kd)
        ka = max(kb-1,1)
        wtb = min(1.0,max(0.0,(zt(k) - zd(ka))/(zd(kb) - zd(ka))))
        call rot_intrp_sclr (data(:,:,ka), xd, yd, id, jd, tmp(:,:,1)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call rot_intrp_sclr (data(:,:,kb), xd, yd, id, jd, tmp(:,:,2)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call intrp_vert (tmp(:,:,1), tmp(:,:,2), wtb, -1e20
     &,   var(:,:,n), imt, jmt)
        call extrap (var(:,:,n), -1.e10, xt, maskv, imt, jmt, 1)

      enddo

!=======================================================================
!     write monthly netcdf surface data
!=======================================================================

      call opennew ("../"//trim(name)//"_mth.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'months')
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      write(name2, '(f8.2)') zt(1)
      long_name = trim(long_name)//' at '//trim(name2)//' m'
      stand_name = trim(stand_name)//' at '//trim(name2)//' m'
      call defvar (trim(name), iou, 3, (/id_xt,id_yt,id_time/)
     &, vmin, vmax,' ', 'F', trim(long_name), trim(stand_name)
     &, trim(units))
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, real(n), 1., 0.)
        call putvaramsk (trim(name), iou, imt*jmt, (/1,1,n/)
     &,   (/imt,jmt,1/), var(:,:,n), vmask(:,:,1), 1000., 0.)
      enddo
      call closefile (iou)

      end
