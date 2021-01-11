      program precip_nc

!=======================================================================
!     creates precip files precip_mth.nc and precip_ann.nc
!=======================================================================

      implicit none

      real, allocatable :: avar(:,:), var(:,:,:), xt(:), yt(:)

      integer id, jd
      parameter (id=192, jd=94)
      real data(id,jd), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, n, id_xt, id_yt, id_time
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
      path = trim(path)//'/ncep/surface_gauss/'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( avar(imt,jmt) )
      allocate ( var(imt,jmt,12) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

      avar(:,:) = 0.
      do n=1,12

!=======================================================================
!     read data
!=======================================================================

        call openfile (trim(path)//'prate.mon.ltm.nc', iou)
        call getvara ('lon', iou, id, (/1/), (/id/), xd, 1., 0.)
        call getvara ('lat', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!       flip latitude
        yd(1:jd) = yd(jd:1:-1)
        call getvara ('prate', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:), 1., 0.)
!       flip data in latitude
        data(:,1:jd) = data(:,jd:1:-1)
        call closefile (iou)
      
!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_sclr (data, xd, yd, id, jd, var(:,:,n)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

!=======================================================================
!     set cyclic boundary condition
!=======================================================================
     
        var(1,:,n) = var(imt-1,:,n)
        var(imt,:,n) = var(2,:,n)

        avar(:,:) = avar(:,:) + var(:,:,n)

      enddo
      avar(:,:) = avar(:,:)/12.
     
!=======================================================================
!       write monthly netcdf data
!=======================================================================

      call opennew ("../precip_mth.nc", iou)
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
      call defvar ('precip', iou, 3, (/id_xt,id_yt,id_time/), -1., 1.
     &, ' ', 'F', 'precipitation rate', ' ' ,'kg m-2 s-1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, real(n), 1., 0.)
        call putvara ('precip', iou, imt*jmt, (/1,1,n/),(/imt,jmt,1/)
     &, var(:,:,n), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!       write annual netcdf data
!=======================================================================

      call opennew ("../precip_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('precip', iou, 2, (/id_xt,id_yt/), -1., 1.
     &, ' ', 'F', 'precipitation rate', ' ' ,'kg m-2 s-1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('precip', iou, imt*jmt, (/1,1/),(/imt,jmt/)
     &, avar(:,:), 1., 0.)
      call closefile (iou)

      end
