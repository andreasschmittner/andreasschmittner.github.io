      program elev_nc

!=======================================================================
!     creates elevation data file elev.nc
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:), xt(:), yt(:)

      integer id, jd
      parameter (id=360, jd=180)
      real data(id,jd), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, n, id_xt, id_yt
      real psi, theta, phi, total

      logical exists
      
      character(120) :: path
      
      print*, "Warning: This may not reproduce original elevation data"
      
!=======================================================================
!     set path and read path file if it exists
!=======================================================================

      path = '/usr/local/models/UVic_ESCM/data_source'
      inquire (file='../path', exist=exists)
      if (exists) then
        open (10,file='../path')
        read (10,'(a)') path
      endif
      path = trim(path)//'/topog/'

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

      call openfile ("../kmt.nc", iou)
      call getvara ('kmt', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:), 1., 0.)
      maskv(:,:) = nint(var(:,:))
      call closefile (iou)

!=======================================================================
!     read elevation data
!=======================================================================

      call openfile (trim(path)//'elev.1-deg.nc', iou)
      call getvara ('lon', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('lat', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!     flip latitude
      yd(1:jd) = yd(jd:1:-1)
      call getvara ('data', iou, id*jd, (/1,1,1/), (/id,jd,1/)
     &, data(:,:), 1., 0.)
!     flip data in latitude
      data(:,1:jd) = data(:,jd:1:-1)
      call closefile (iou)
        
!=======================================================================
!     rotate and interpolate data
!=======================================================================

      call rot_intrp_sclr (data, xd, yd, id, jd, var
     &, xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
     
!=======================================================================
!     zero elevation over water and smooth data
!=======================================================================

      where (maskv(:,:) .ne. 0) var(:,:) = 0.
      where (var(:,:) .le. 0) var(:,:) = 0.

      call smooth (var, imt, jmt, 4)

      where (maskv(:,:) .ne. 0) var(:,:) = 0.
      where (var(:,:) .le. 0) var(:,:) = 0.

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:) = var(imt-1,:)
      var(imt,:) = var(2,:)

!=======================================================================
!     write netcdf elev data
!=======================================================================

      call opennew ("../elev.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('elev', iou, 2, (/id_xt,id_yt/), -1.e4, 1.e4, ' '
     &, 'F', 'land elevation and ocean depth', 'surface_altitude', 'm')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('elev', iou, imt*jmt, (/1,1/)
     &,            (/imt,jmt/), var, 1., 0.)
      call closefile (iou)

      end


      subroutine smooth (var, imt, jmt, num)
!=======================================================================
!     smooth data (5 point smoothing)
!=======================================================================
      
      implicit none
      
      integer i, imt, j, jmt, n, num
      
      real var(imt,jmt)
      real, allocatable :: tmp(:,:)

      allocate ( tmp(imt,jmt) )
      
      print*, 'Warning: applying smoothing'

      var(1,:) = var(imt-1,:)
      var(imt,:) = var(2,:)
      var(:,1) = var(:,2)
      var(:,jmt) = var(:,jmt-1)        
      do n=1,num
        tmp(:,:) = var(:,:)
        do j=2,jmt-1
          do i=2,imt-1
            var(i,j) = 0.5*tmp(i,j) + 0.125*(tmp(i+1,j) + tmp(i-1,j)
     &               + tmp(i,j+1) + tmp(i,j-1))
          enddo
        enddo
        var(1,:) = var(imt-1,:)
        var(imt,:) = var(2,:)
        var(:,1) = var(:,2)
        var(:,jmt) = var(:,jmt-1)        
      enddo
     
      return
      end            
