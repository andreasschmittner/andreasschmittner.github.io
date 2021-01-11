      program grid_nc

!=======================================================================
!     creates grid data file grid.nc
!=======================================================================

      implicit none

      integer imt, jmt, km
      parameter (imt=102, jmt=102, km=19)

      real var(imt,jmt,5), psi, theta, phi, ddzw
      real xt(imt), xu(imt), yt(jmt), yu(jmt), zt(km), zw(km)
      real dxt(imt), dyt(jmt), dxu(imt), dyu(jmt), dzt(km), dzw(km)
      real psi, theta, phi

      integer id_xt, id_yt, id_zt, id_xu, id_yu, id_zw
      integer id_dxt, id_dyt, id_dzt, id_dxu, id_dyu, id_dzw
      integer id_phi, id_theta, id_psi
      integer iou, i, j, k, n

!=======================================================================
!     define rotation in radians
!=======================================================================

!     rotate
!      phi   = -2.268928028
!      theta = -0.226892803
!      psi   = 3.141592654
!     unrotate
!      psi   = 2.268928028
!      theta = 0.226892803
!      phi   = -3.141592654
!     no rotation
      psi   = 0.
      theta = 0.
      phi   = 0.

!=======================================================================
!     define grid
!=======================================================================

!     parabolic depth distribution
      zw(1) = 50.  ! bottom depth of the first "t" grid box (m)
      zw(km) = 5400.  ! bottom depth of the last "t" grid box (m)
      ddzw = 30.  ! second derivative (26. gives ~5400m for km=19)
!      ddzw = 2.*(zw(km)/float(km) - zw(1))/float(km-1)
      do k=1,km
        dzt(k) = zw(1) + (float(k) - 1.0)*ddzw
        dzw(k) = zw(1) + (float(k) - 0.5)*ddzw
      enddo
      zt(1) = zw(1) - dzw(1)/2.   ! for a "u" centred grid
!      zt(1) = zw(1) - dzt(1)/2.   ! for a "t" centred grid
      do k=2,km
        zt(k) = zt(k-1) + dzw(k-1)
        zw(k) = zw(k-1) + dzt(k)
      enddo
      
!     constant 3.6 (longitude) by 1.8 (latitude) degree grid
      xt(1) = -1.8
      dxt(:) = 3.6
      yt(1) = -90.9
      dyt(:) = 1.8

      dxt(1) = dxt(imt-1)
      dxt(imt) = dxt(2)
      do i=1,imt-1
        dxu(i) = (dxt(i+1) + dxt(i))/2.
      enddo
      dxu(1) = dxu(imt-1)
      dxu(imt) = dxu(2)
      do i=2,imt
        xt(i) = xt(i-1) + dxu(i)
      enddo
      xt(1) = xt(imt-1)-360.
      xt(imt) = xt(2)+360.
      do i=1,imt
        xu(i) = xt(i) + dxt(i)/2.
      enddo
      xu(1) = xu(imt-1)-360.
      xu(imt) = xu(2)+360.

      do j=1,jmt-1
        dyu(j) = (dyt(j+1) + dyt(j))/2.
      enddo
      dyu(jmt) = dyu(jmt-1)
      do j=2,jmt
        yt(j) = yt(j-1) + dyu(j)
      enddo
      do j=1,jmt
        yu(j) = yt(j) + dyt(j)/2.
      enddo

!=======================================================================
!     define latitude and longitude
!=======================================================================

      do j=1,jmt
        var(:,j,5) = yt(j)
      enddo
      call rot_intrp_sclr (var(:,:,5), xt, yt, imt, jmt, var(:,:,1)
     &,  xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

      do j=1,jmt
        var(:,j,5) = yu(j)
      enddo
      call rot_intrp_sclr (var(:,:,5), xt, yt, imt, jmt, var(:,:,2)
     &,  xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

      do i=1,imt
        var(i,:,5) = xt(i)
      enddo
      call rot_intrp_sclr (var(:,:,5), xt, yt, imt, jmt, var(:,:,3)
     &,  xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

      do i=1,imt
        var(i,:,5) = xu(i)
      enddo
      call rot_intrp_sclr (var(:,:,5), xt, yt, imt, jmt, var(:,:,4)
     &,  xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)

      var(1,:,:) = var(imt-1,:,:)
      var(imt,:,:) = var(2,:,:)
      var(:,jmt,:) = var(:,jmt-1,:)
      var(:,1,:) = var(:,2,:)      

!=======================================================================
!     write grid data
!=======================================================================

!     write netcdf data
      call opennew ("../grid.nc", iou)
      call redef (iou)
      call defdim ('psi', iou, 1, id_psi)
      call defdim ('theta', iou, 1, id_theta)
      call defdim ('phi', iou, 1, id_phi)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defdim ('zt', iou, km, id_zt)
      call defdim ('xu', iou, imt, id_xu)
      call defdim ('yu', iou, jmt, id_yu)
      call defdim ('zw', iou, km, id_zw)
      call defdim ('xt', iou, imt, id_dxt)
      call defdim ('yt', iou, jmt, id_dyt)
      call defdim ('zt', iou, km, id_dzt)
      call defdim ('xu', iou, imt, id_dxu)
      call defdim ('yu', iou, jmt, id_dyu)
      call defdim ('zw', iou, km, id_dzw)
      call defvar ('psi', iou, 1, (/id_psi/), 0., 0., ' ', 'F'
     &, 'rotation angle psi', ' ', 'radians')
      call defvar ('theta', iou, 1, (/id_theta/), 0., 0., ' ', 'F'
     &, 'rotation angle theta', ' ', 'radians')
      call defvar ('phi', iou, 1, (/id_phi/), 0., 0., ' ', 'F'
     &, 'rotation angle phi', ' ', 'radians')
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('zt', iou, 1, (/id_zt/), 0., 0., 'Z', 'F'
     &, 'depth of the t grid', 'depth', 'm')
      call defvar ('xu', iou, 1, (/id_xu/), 0., 0., 'X', 'F'
     &, 'longitude of the u grid', 'longitude', 'degrees_east')
      call defvar ('yu', iou, 1, (/id_yu/), 0., 0., 'Y', 'F'
     &, 'latitude of the u grid', 'latitude', 'degrees_north')
      call defvar ('zw', iou, 1, (/id_zw/), 0., 0., 'Z', 'F'
     &, 'depth of the w grid', 'depth', 'm')
      call defvar ('dxt', iou, 1, (/id_dxt/), 0., 0., ' ', 'F'
     &, 't grid width in x', '', 'degrees')
      call defvar ('dyt', iou, 1, (/id_dyt/), 0., 0., ' ', 'F'
     &, 't grid hieght in y', '', 'degrees')
      call defvar ('dzt', iou, 1, (/id_dzt/), 0., 0., ' ', 'F'
     &, 't grid thickness in z', '', 'm')
      call defvar ('dxu', iou, 1, (/id_dxu/), 0., 0., ' ', 'F'
     &, 'u grid width in x', '', 'degrees')
      call defvar ('dyu', iou, 1, (/id_dyu/), 0., 0., ' ', 'F'
     &, 'u grid hieght in y', '', 'degrees')
      call defvar ('dzw', iou, 1, (/id_dzw/), 0., 0., ' ', 'F'
     &, 'w grid thickness in z', '', 'm')
      call defvar ('lat_t', iou, 2, (/id_xt,id_yt/), -1.e2, 1.e2
     &, '', 'F', 'latitude of the t grid', '', 'degrees')
      call defvar ('lat_u', iou, 2, (/id_xu,id_yu/), -1.e2, 1.e2
     &, '', 'F', 'latitude of the u grid', '', 'degrees')
      call defvar ('lon_t', iou, 2, (/id_xt,id_yt/), -4.e2, 4.e2
     &, '', 'F', 'longitude of the t grid', '', 'degrees')
      call defvar ('lon_u', iou, 2, (/id_xu,id_yu/), -4.e2, 4.e2
     &, '', 'F', 'longitude of the u grid', '', 'degrees')
      call enddef (iou)
      call putvars ('psi', iou, 1, psi, 1., 0.)
      call putvars ('theta', iou, 1, theta, 1., 0.)
      call putvars ('phi', iou, 1, phi, 1., 0.)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('zt', iou, km, (/1/), (/km/), zt, 1., 0.)
      call putvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call putvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      call putvara ('zw', iou, km, (/1/), (/km/), zw, 1., 0.)
      call putvara ('dxt', iou, imt, (/1/), (/imt/), dxt, 1., 0.)
      call putvara ('dyt', iou, jmt, (/1/), (/jmt/), dyt, 1., 0.)
      call putvara ('dzt', iou, km, (/1/), (/km/), dzt, 1., 0.)
      call putvara ('dxu', iou, imt, (/1/), (/imt/), dxu, 1., 0.)
      call putvara ('dyu', iou, jmt, (/1/), (/jmt/), dyu, 1., 0.)
      call putvara ('dzw', iou, km, (/1/), (/km/), dzw, 1., 0.)
      call putvara ('lat_t', iou, imt*jmt, (/1,1/)
     &,               (/imt,jmt/), var(1,1,1), 1., 0.)
      call putvara ('lat_u', iou, imt*jmt, (/1,1/)
     &,               (/imt,jmt/), var(1,1,2), 1., 0.)
      call putvara ('lon_t', iou, imt*jmt, (/1,1/)
     &,               (/imt,jmt/), var(1,1,3), 1., 0.)
      call putvara ('lon_u', iou, imt*jmt, (/1,1/)
     &,               (/imt,jmt/), var(1,1,4), 1., 0.)
      call closefile (iou)

      end
