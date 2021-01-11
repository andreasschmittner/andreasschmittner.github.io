      program ice4g_nc

!=======================================================================
!     creates ice data file ice4g.nc
!=======================================================================

      implicit none

      real, allocatable :: var(:,:,:)
      real, allocatable :: xt(:), yt(:)

      integer id, jd
      parameter (id=360, jd=180)
      real data(id,jd,3), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, k, n, id_xt, id_yt, id_time
      real psi, theta, phi, time

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
      path = trim(path)//'/ice/'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( var(imt,jmt,2) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read present day topography data
!=======================================================================

      call openfile (trim(path)//'top.nc', iou)
      call getvara ('X', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('Y', iou, jd, (/1/), (/jd/), yd, 1., 0.)
      call getvara ('topography', iou, id*jd, (/1,1,22/),
     &  (/id,jd,1/), data(:,:,3), 1., 0.)
!     flip latitude
      yd(1:jd) = yd(jd:1:-1)
!     flip data in latitude
      data(:,1:jd,3) = data(:,jd:1:-1,3)
      call closefile (iou)

!=======================================================================
!     define netcdf file
!=======================================================================

      call opennew ("../ice4g.nc", iou)
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
      call defvar ('hicel', iou, 3, (/id_xt,id_yt,id_time/), 0., 1.e10
     &, ' ', 'F', 'elevation difference from present day', '', 'm')
      call defvar ('aicel', iou, 3, (/id_xt,id_yt,id_time/), 0., 1.
     &, ' ', 'F', 'mask of ice area', '', '')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call closefile (iou)

      do n=1,22

!=======================================================================
!       read data
!=======================================================================

        call openfile (trim(path)//'top.nc', iou)
        call getvara ('topography', iou, id*jd, (/1,1,n/)
     &,   (/id,jd,1/), data(:,:,1), 1., 0.)
        call closefile (iou)
        call openfile (trim(path)//'ice.nc', iou)
        call getvara ('ice', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,2), 1., 0.)
        call closefile (iou)
!       flip data in latitude
        data(:,1:jd,1:2) = data(:,jd:1:-1,1:2)

!=======================================================================
!       subtract off present day topography
!=======================================================================

        data(:,:,1) = (data(:,:,1) - data(:,:,3))*data(:,:,2)
        where (data(:,:,1) < 0.) data(:,:,1) = 0.

!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_sclr (data(:,:,1), xd, yd, id, jd, var(:,:,1)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

        call rot_intrp_sclr (data(:,:,2), xd, yd, id, jd, var(:,:,2)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

!=======================================================================
!       set cyclic boundary condition
!=======================================================================

        var(1,:,:) = var(imt-1,:,:)
        var(imt,:,:) = var(2,:,:)

!=======================================================================
!       write  netcdf data
!=======================================================================

        time = float(n - 20)*1000.
        call openfile ("../ice4g.nc", iou)
        call putvars ('time', iou, n, time, 1., 0.)
        call putvara ('hicel', iou, imt*jmt, (/1,1,n/)
     &,              (/imt,jmt,1/), var(:,:,1), 1., 0.)
        call putvara ('aicel', iou, imt*jmt, (/1,1,n/)
     &,              (/imt,jmt,1/), var(:,:,2), 1., 0.)
        call closefile (iou)

      enddo

      end
