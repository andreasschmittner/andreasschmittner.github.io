      program dtr

!=======================================================================
!     creates diurnal temperature range data
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: avar(:,:), var(:,:,:), vmask(:,:), lat_t(:,:)
      real, allocatable :: xt(:), yt(:)

      integer id, jd
      parameter (id=720, jd=360)
      real data(id,jd), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, k, n, id_xt, id_yt, id_time
      real psi, theta, phi, rad, alat
      
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
      path = trim(path)//'/dtr/'
      rad = acos(-1.)/180.

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( var(imt,jmt,12) )
      allocate ( avar(imt,jmt) )
      allocate ( maskv(imt,jmt) )
      allocate ( vmask(imt,jmt) )
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
     &, var(:,:,1), 1., 0.)
      maskv(:,:) = nint(var(:,:,1))
      call closefile (iou)

!=======================================================================
!     read data
!=======================================================================

      call openfile (trim(path)//'dtr.nc', iou)
      call getvara ('X', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('Y', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!     flip latitude
      yd(1:jd) = yd(jd:1:-1)
      call closefile (iou)

      do n=1,12

!=======================================================================
!       read data
!=======================================================================

        call openfile (trim(path)//'dtr.nc', iou)
        call getvara ('dtr', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:), 1., 0.)
!       flip data in latitude
        data(:,1:jd) = data(:,jd:1:-1)
        call closefile (iou)

        do j=1,jd
          do i=1,id
            if (data(i,j) .lt. -999.) data(i,j) = 2.e20
          enddo
        enddo

!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_sclr (data(:,:), xd, yd, id, jd, var(:,:,n)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call extrap (var(:,:,n), -1.e10, xt, maskv, imt, jmt, 1)

      enddo

      avar(:,:) = 0.
      print*, 'Warning: setting dtr to 20 above 85 N and below 60 S'
      do n=1,12
        where (lat_t(:,:) .ge. 85) var(:,:,n) = 20.
        where (lat_t(:,:) .le. -60) var(:,:,n) = 20.
      
        avar(:,:) = avar(:,:) + var(:,:,n)
      enddo
      avar(:,:) = avar(:,:)/12

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:,:) = var(imt-1,:,:)
      var(imt,:,:) = var(2,:,:)

      vmask(:,:) = 0.
      where (maskv(:,:) .eq. 0) vmask(:,:) = 1.

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call opennew ("../dtr_mth.nc", iou)
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
      call defvar ('dtr', iou, 3, (/id_xt,id_yt,id_time/), 0., 40.
     &, ' ', 'F', 'diurnal temperature range', '', 'K')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, real(n), 1., 0.)
        call putvaramsk ('dtr', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(1,1,n), vmask, 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../dtr_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('dtr', iou, 2, (/id_xt,id_yt,id_time/), 0., 40.
     &, ' ', 'F', 'diurnal temperature range', '', 'K')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvaramsk ('dtr', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar, vmask, 1., 0.)
      call closefile (iou)

      end
