      program tmsk_nc

!=======================================================================
!     creates tracer grid mask data file tmsk.nc
!=======================================================================

      implicit none

      real, allocatable :: var(:,:,:), xt(:), yt(:), lat_t(:,:)

      integer i, imt, iou, j, jmt, n, id_xt, id_yt

      real psi, theta, phi, year

      logical exists
      
!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( var(imt,jmt,2) )
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
      where (var(:,:,1) .gt. 0.) var(:,:,1) = 1.

!=======================================================================
!     read ice data
!=======================================================================
   
      inquire (file='../ice4g.nc', exist=exists)
      if (exists) then
        year = 2000.
        print*, 'Warning: year for ice shelves = ',year
        call opentime ('../ice4g.nc', year, n, iou)
        call getvara ('aicel', iou, imt*jmt, (/1,1,n/)
     &,   (/imt,jmt,1/), var(:,:,2), 1., 0.)
        print*, 'Warning: adding ice shelves around Antarctica'
        where (var(:,:,2) .gt. 0.5 .and. lat_t(:,:) .lt. -55.) 
     &    var(:,:,1) = 0.
      endif

!=======================================================================
!     fill isolated bays
!=======================================================================

      call fill_bays (var(:,:,1), imt, jmt)

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:,1) = var(imt-1,:,1)
      var(imt,:,1) = var(2,:,1)

!=======================================================================
!     write netcdf tmsk data
!=======================================================================

      call opennew ("../tmsk.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('tmsk', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'I', 'ocean mask', ' ' ,'')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('tmsk', iou, imt*jmt, (/1,1/)
     &,            (/imt,jmt/), var(:,:,1), 1., 0.)
      call closefile (iou)

      end


      subroutine fill_bays (var, imt, jmt)
!=======================================================================
!     fill isolated bays
!=======================================================================
      
      implicit none
      
      integer i, imt, j, jmt
      
      real var(imt,jmt)
      real, allocatable :: tmp(:,:)

      allocate ( tmp(imt,jmt) )

      var(1,:) = var(imt-1,:)
      var(imt,:) = var(2,:)
      var(:,1) = var(:,2)
      var(:,jmt) = var(:,jmt-1)         
      do j=1,jmt-1
        do i=1,imt-1
          tmp(i,j) = min(var(i,j), var(i+1,j), var(i,j+1), var(i+1,j+1))
        enddo
      enddo
      tmp(1,:) = tmp(imt-1,:)
      tmp(imt,:) = tmp(2,:)
      tmp(:,1) = tmp(:,2)
      tmp(:,jmt) = tmp(:,jmt-1)      
      do j=2,jmt
        do i=2,imt
          var(i,j) = max(tmp(i,j), tmp(i-1,j), tmp(i,j-1), tmp(i-1,j-1))
        enddo
      enddo
      var(1,:) = var(imt-1,:)
      var(imt,:) = var(2,:)
      var(:,1) = var(:,2)
      var(:,jmt) = var(:,jmt-1)
      
      return
      end            


