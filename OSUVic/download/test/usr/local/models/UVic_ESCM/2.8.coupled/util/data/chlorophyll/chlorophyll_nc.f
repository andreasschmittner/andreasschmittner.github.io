      program chlorophyll_nc

!=======================================================================
!     creates netcdf files of annual surface data
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:,:), tmp(:,:,:), vmask(:,:,:)
      real, allocatable :: xt(:), yt(:), zt(:)

      integer id, jd, kd, kdm
      parameter (id=360, jd=180, kd=7, kdm=7)
      real data(id,jd,kd), xd(id), yd(jd), zd(kd)

      integer i, imt, iou, j, jmt, k, km, n, id_xt, id_yt
      integer id_zt, id_time, ka, kb, m
      real psi, theta, phi, wta, wtb, vmax, vmin
      real missing_value

      logical exists
      
      character(120) :: path, name, name2, long_name, stand_name, units

!=======================================================================
!     set names and values
!=======================================================================

      name = "chlorophyll"
      long_name = "chlorophyll"
      stand_name = "chlorophyll"
      vmax = 9.
      vmin = 0.
      units = "um l-1"
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
!     rotate and interpolate annual model surface data
!=======================================================================

      k=1
      ka = 1
      do m=2,kd
        if (zd(m) .lt. zt(k)) ka = m
      enddo
      kb = min(ka+1,kd)
      ka = max(kb-1,1)
      wtb = min(1.0,max(0.0,(zt(k) - zd(ka))/(zd(kb) - zd(ka))))
      call rot_intrp_sclr (data(:,:,ka), xd, yd, id, jd, tmp(:,:,1)
     &, xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
      call rot_intrp_sclr (data(:,:,kb), xd, yd, id, jd, tmp(:,:,2)
     &, xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
      call intrp_vert (tmp(:,:,1), tmp(:,:,2), wtb, -1e20
     &, var(:,:,k), imt, jmt)
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
!     write annual netcdf surface data
!=======================================================================

      call opennew ("../"//trim(name)//"_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      write(name2, '(f8.2)') zt(1)
      long_name = trim(long_name)//' at '//trim(name2)//' m'
      stand_name = trim(stand_name)//' at '//trim(name2)//' m'
      call defvar (trim(name), iou, 2, (/id_xt,id_yt/)
     &, vmin, vmax,' ', 'F', trim(long_name), trim(stand_name)
     &, trim(units))
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvaramsk (trim(name), iou, imt*jmt*km, (/1,1/)
     &, (/imt,jmt/), var(:,:,:), vmask, 1., 0.)
      call closefile (iou)
      
      end
