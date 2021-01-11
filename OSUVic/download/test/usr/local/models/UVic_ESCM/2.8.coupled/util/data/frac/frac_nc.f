      program frac_nc

!=======================================================================
!     creates vegetation fractions
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:,:), vmask(:,:), lat_t(:,:)
      real, allocatable :: xt(:), yt(:)
      integer nvt
      parameter (nvt=9)

      integer id, jd
      parameter (id=97, jd=73)
      real data(id,jd,nvt), xd(id), yd(jd), dataf(id,jd)

      integer id2, jd2
      parameter (id2=360, jd2=180)
      real data2(id2,jd2), xd2(id2), yd2(jd2)

      integer i, imt, iou, j, jmt, k, n, id_xt, id_yt
      real time, psi, theta, phi
      real x, y, f1, f2, f3, f4, f5, f6, f7, f8, f9

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
      allocate ( var(imt,jmt,nvt) )
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
      maskv(:,:) = 0
      where (vmask(:,:) .eq. 0) maskv(:,:) = 1
      vmask(:,:) = maskv(:,:)
      call closefile (iou)

!=======================================================================
!    read data current vegetation file
!=======================================================================

      print*, trim(path)//'frac_igbp.nc'
      call openfile (trim(path)//'frac_igbp.nc', iou)
      call getvara ('longitude', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('latitude', iou, jd, (/1/), (/jd/), yd, 1., 0.)
      call getvara ('BTREE', iou, id*jd, (/1,1/), (/id,jd/)
     &, data(:,:,1), 1., 0.)
      call getvara ('NTREE', iou, id*jd, (/1,1/), (/id,jd/)
     &, data(:,:,2), 1., 0.)
      call getvara ('C3G', iou, id*jd, (/1,1/), (/id,jd/)
     &, data(:,:,3), 1., 0.)
      call getvara ('C4G', iou, id*jd, (/1,1/), (/id,jd/)
     &, data(:,:,4), 1., 0.)
      call getvara ('SHRUB', iou, id*jd, (/1,1/), (/id,jd/)
     &, data(:,:,5), 1., 0.)
      call getvara ('SOIL', iou, id*jd, (/1,1/), (/id,jd/)
     &, data(:,:,6), 1., 0.)
      call getvara ('URBAN_and_WATER', iou, id*jd, (/1,1/)
     &, (/id,jd/), data(:,:,7), 1., 0.)
      call getvara ('ICE', iou, id*jd, (/1,1/), (/id,jd/)
     &, data(:,:,9), 1., 0.)
      call closefile (iou)

!     combine urban, water and ice into "other"
      data(:,:,7) = data(:,:,7) + data(:,:,9)

!=======================================================================
!     rotate and interpolate data
!=======================================================================

      do n=1,nvt
        call rot_intrp_sclr (data(:,:,n), xd, yd, id, jd, var(:,:,n)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 1)
        call extrap2 (var(:,:,n), -1.e10, xt, imt, jmt)
      enddo
      
!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:,:) = var(imt-1,:,:)
      var(imt,:,:) = var(2,:,:)

!=======================================================================
!     write netcdf veg data
!=======================================================================

      call opennew ("../frac.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('BTREE', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'F', 'Broadleaf tree vegetation fraction', '', '1')
      call defvar ('NTREE', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'F', 'Needleleaf tree vegetation fraction', '', '1')
      call defvar ('C3G', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'F', 'C3 grass vegetation fraction', '', '1')
      call defvar ('C4G', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'F', 'C4 grass vegetation fraction', '', '1')
      call defvar ('SHRUB', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'F', 'Shrub vegetation fraction', '', '1')
      call defvar ('SOIL', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'F', 'Soil vegetation fraction', '', '1')
      call defvar ('OTHER', iou, 2, (/id_xt,id_yt/), 0., 1.
     &, ' ', 'F', 'Other vegetation fraction', '', '1')
!      call defvar ('URBAN', iou, 2, (/id_xt,id_yt/), 0., 1.
!     &, ' ', 'F', 'Urban vegetation fraction', '', '1')
!      call defvar ('WATER', iou, 2, (/id_xt,id_yt/), 0., 1.
!     &, ' ', 'F', 'Water vegetation fraction', '', '1')
!      call defvar ('ICE', iou, 2, (/id_xt,id_yt/), 0., 1.
!     &, ' ', 'F', 'Ice vegetation fraction', '', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvaramsk ('BTREE', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,1), vmask, 1., 0.)
      call putvaramsk ('NTREE', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,2), vmask, 1., 0.)
      call putvaramsk ('C3G', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,3), vmask, 1., 0.)
      call putvaramsk ('C4G', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,4), vmask, 1., 0.)
      call putvaramsk ('SHRUB', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,5), vmask, 1., 0.)
      call putvaramsk ('SOIL', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,6), vmask, 1., 0.)
      call putvaramsk ('OTHER', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,7), vmask, 1., 0.)
!      call putvaramsk ('WATER', iou, imt*jmt, (/1,1/), (/imt,jmt/)
!     &, var(:,:,7), vmask, 1., 0.)
!      call putvaramsk ('SOIL', iou, imt*jmt, (/1,1/), (/imt,jmt/)
!     &, var(:,:,8), vmask, 1., 0.)
!      call putvaramsk ('ICE', iou, imt*jmt, (/1,1/), (/imt,jmt/)
!     &, var(:,:,9), vmask, 1., 0.)
      call closefile (iou)

      end
