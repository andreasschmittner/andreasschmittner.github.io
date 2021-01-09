      program veg_class_nc

!=======================================================================
!     creates current vegetation class file
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:), vmask(:,:), lat_t(:,:)
      real, allocatable :: xt(:), yt(:)
      real nvt
      parameter (nvt=7)

      integer id, jd
      parameter (id=360, jd=180)
      real data(id,jd), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, k, n, id_xt, id_yt
      real psi, theta, phi

      logical exists
      
      character(120) :: path
      
!=======================================================================
!     set path and read path file if it exists
!=======================================================================

      path = '/usr/local/models/UVic_ESCM/data_source'
      inquire (file='../path', exist=exists)
      if (exists) then
        open (10,file='../path')
        read (10,'(a)') path
      endif
      path = trim(path)//'/veg/'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( var(imt,jmt) )
      allocate ( vmask(imt,jmt) )
      allocate ( maskv(imt,jmt) )
      allocate ( lat_t(imt,jmt) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call getvara ('lat_t', iou, imt*jmt, (/1,1/), (/imt,jmt/), lat_t
     &, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read kmt data
!=======================================================================

      call openfile ("../kmt.nc", iou)
      call getvara ('kmt', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, vmask(:,:), 1., 0.)
      maskv(:,:) = nint(vmask(:,:))
      call closefile (iou)

!=======================================================================
!     read data
!=======================================================================

      call openfile (trim(path)//'veg_class.nc', iou)
      call getvara ('X', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('Y', iou, jd, (/1/), (/jd/), yd, 1., 0.)
      call getvara ('class', iou, id*jd, (/1,1,1,1/), (/id,jd,1,1/)
     &, data(:,:), 1., 0.)
!     flip data in latitude
      yd(1:jd) = yd(jd:1:-1)
      data(:,1:jd) = data(:,jd:1:-1)
      call closefile (iou)

!=======================================================================
!     reduce vegetation classes to nvt
!=======================================================================

      do j=1,jd
        do i=1,id
          if (data(i,j) .eq. 1) data(i,j) = 1.
          if (data(i,j) .eq. 2) data(i,j) = 2.
          if (data(i,j) .eq. 3) data(i,j) = 2.
          if (data(i,j) .eq. 4) data(i,j) = 2.
          if (data(i,j) .eq. 5) data(i,j) = 2.
          if (data(i,j) .eq. 6) data(i,j) = 3.
          if (data(i,j) .eq. 7) data(i,j) = 3.
          if (data(i,j) .eq. 8) data(i,j) = 3.
          if (data(i,j) .eq. 9) data(i,j) = 4.
          if (data(i,j) .eq. 10) data(i,j) = 5.
          if (data(i,j) .eq. 11) data(i,j) = 6.
          if (data(i,j) .eq. 12) data(i,j) = 3.
          if (data(i,j) .eq. 13) data(i,j) = 7.
          if (data(i,j) .eq. 14) data(i,j) = 3.
          if (data(i,j) .eq. 15) data(i,j) = 3.
        enddo
      enddo
      
      where (data(:,:) .lt. 1 .or. data(:,:) .gt. nvt) data(:,:) = 1.e20

!=======================================================================
!       rotate and interpolate data
!=======================================================================

      call rot_intrp_sclr (data(:,:), xd, yd, id, jd, var(:,:)
     &, xt, yt, imt, jmt, phi, theta, psi, -1.e20, 1)
      call extrap2 (var(:,:), -1.e10, xt, imt, jmt)

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:) = var(imt-1,:)
      var(imt,:) = var(2,:)
      
      print*, 'Warning: setting vegetation class to 7 over Arctica'
      where (lat_t(:,:) .ge. 85 .and. maskv(:,:) .eq. 0) var(:,:) = 7

      vmask(:,:) = 0.
      where (maskv(:,:) .eq. 0) vmask(:,:) = 1.

!=======================================================================
!     write netcdf veg data
!=======================================================================

      call opennew ("../veg_class.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('veg', iou, 2, (/id_xt,id_yt/), 0., nvt
     &, ' ', 'F', 'present day vegetation class', ' ', '1')
      call putatttext (iou, 'veg', 'type1', 'tropical forest')
      call putatttext (iou, 'veg', 'type2', 'temperate/boreal forest')
      call putatttext (iou, 'veg', 'type3', 'grass')
      call putatttext (iou, 'veg', 'type4', 'shrub')
      call putatttext (iou, 'veg', 'type5', 'tundra')
      call putatttext (iou, 'veg', 'type6', 'desert')
      call putatttext (iou, 'veg', 'type7', 'ice')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvaramsk ('veg', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var, vmask, 1., 0.)
      call closefile (iou)

      end
