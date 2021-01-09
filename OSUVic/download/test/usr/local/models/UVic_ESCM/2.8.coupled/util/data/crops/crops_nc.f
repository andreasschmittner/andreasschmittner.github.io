      program crops_nc

!=======================================================================
!     creates croplands data file crops.nc
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:)
      real, allocatable :: xt(:), yt(:)
      real ntr
      parameter (ntr=293)

      integer id, jd
      parameter (id=720, jd=360)
      real data(id,jd), xd(id), yd(jd)

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
      allocate ( maskv(imt,jmt) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read kmt data
!=======================================================================

!     read kmt for masking
      call openfile ("../kmt.nc", iou)
      call getvara ('kmt', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:), 1., 0.)
      maskv(:,:) = nint(var(:,:))
      call closefile (iou)

!=======================================================================
!     read crop dataset grid
!=======================================================================

      call openfile (trim(path)//'crop_1700-1992_v1.0.nc', iou)
      call getvara ('longitude', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('latitude', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!     flip data in latitude
      yd(1:jd) = yd(jd:1:-1)
      call closefile (iou)

!=======================================================================
!     define croplands file
!=======================================================================

      call opennew ("../crops.nc", iou)
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
      call defvar ('farea', iou, 3, (/id_xt,id_yt,id_time/), 0., 1.
     &, ' ', 'F', 'crop area', '', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call closefile (iou)

      do n=1, ntr
      
!=======================================================================
!       read crop dataset data
!=======================================================================

        call openfile (trim(path)//'crop_1700-1992_v1.0.nc', iou)
        call getvars ('time', iou, n, time, 1., 0.)
        call getvara ('farea', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:), 1., 0.)
!       flip data in latitude
        data(:,1:jd) = data(:,jd:1:-1)
        call closefile (iou)

!=======================================================================
!       rotate and interpolate data
!=======================================================================


        call rot_intrp_sclr (data(:,:), xd, yd, id, jd, var(:,:)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 1)
        where (maskv(:,:) .ne. 0 .or. var(:,:) .gt. 1.) var(:,:) = 0     

!=======================================================================
!       set cyclic boundary condition
!=======================================================================

        var(1,:) = var(imt-1,:)
        var(imt,:) = var(2,:)

!=======================================================================
!       write netcdf croplands data
!=======================================================================

        call openfile ("../crops.nc", iou)
        call putvars ('time', iou, n, time, 1., 0.)
        call putvara ('farea', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(:,:), 1., 0.)
        call closefile (iou)

      enddo

      end
