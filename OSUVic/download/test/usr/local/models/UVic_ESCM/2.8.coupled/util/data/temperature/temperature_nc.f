      program temperature_nc

!=======================================================================
!     creates netcdf files of annual data and monthly surface data
!=======================================================================

      implicit none

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:,:,:), tmp(:,:,:), vmask(:,:,:)
      real, allocatable :: xt(:), yt(:), zt(:)

      integer id, jd, kd, kdm
      parameter (id=360, jd=180, kd=33, kdm=7)
      real data(id,jd,kd), xd(id), yd(jd), zd(kd)

      integer i, imt, iou, j, jmt, k, km, n, id_xt, id_yt
      integer id_zt, id_time, ka, kb, m
      real psi, theta, phi, wta, wtb, vmax, vmin
      real missing_value, daymon(12)
      real (kind=8) :: t, s, p, theta0

      logical exists
      
      character(120) :: path, name, name2, long_name, stand_name, units

      data daymon(1:6)  / 15.5,  45.,   74.5, 105.,  135.5, 166. /
      data daymon(7:12) /196.5, 227.5, 258.,  288.5, 319.,  349.5/

!=======================================================================
!     set names and values
!=======================================================================

      name = "temperature"
      vmax = 310.
      vmin = 270.
      units = "K"
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
      allocate ( var(imt,jmt,km,2) )
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
!     rotate and interpolate data for all levels
!=======================================================================

      do k=1,km
        ka = 1
        do m=2,kd
          if (zd(m) .lt. zt(k)) ka = m
        enddo
        kb = min(ka+1,kd)
        ka = max(kb-1,1)
        wtb = min(1.0,max(0.0,(zt(k) - zd(ka))/(zd(kb) - zd(ka))))
        call rot_intrp_sclr (data(:,:,ka), xd, yd, id, jd, tmp(:,:,1)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call rot_intrp_sclr (data(:,:,kb), xd, yd, id, jd, tmp(:,:,2)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call intrp_vert (tmp(:,:,1), tmp(:,:,2), wtb, -1e20
     &,   var(:,:,k,1), imt, jmt)
      enddo
      call extrap (var(:,:,:,1), -1.e20, xt, maskv(:,:), imt, jmt, km)

!=======================================================================
!       set cyclic boundary condition
!=======================================================================

      var(1,:,:,1) = var(imt-1,:,:,1)
      var(imt,:,:,1) = var(2,:,:,1)
      
!=======================================================================
!     calculate potential temperature from in situ
!=======================================================================

      inquire (file='../salinity_ann.nc', exist=exists)
      if (exists) then
        call openfile ("../salinity_ann.nc", iou)
        call getvara ('salinity', iou, imt*jmt*km, (/1,1,1/)
     &,   (/imt,jmt,km/), var(:,:,:,2), 1., 0.)
        call closefile (iou)

        do k=1,km
          do j=1,jmt
            do i=1,imt
              t = var(i,j,k,1)
              s = var(i,j,k,2)
              p = zt(k)
              call potem (t, s, p, theta0)
              var(i,j,k,1) = theta0
            enddo
          enddo
        enddo
        long_name = "potential temperature"
        stand_name = "sea_water_potential_temperature"
      else
        print*, 'Warning: salinity_ann.nc missing, in situ temperature'
        long_name = "in situ temperature"
        stand_name = "sea_water_temperature"
      endif

!=======================================================================
!     create 3D mask
!=======================================================================

      vmask(:,:,:) = 0
      do k=1,km
        where (maskv(:,:) .ge. k) vmask(:,:,k) = 1
      enddo

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../"//trim(name)//"_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defdim ('zt', iou, km, id_zt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('zt', iou, 1, (/id_zt/), 0., 0., 'Z', 'F'
     &, 'depth', 'depth', 'm')
      call defvar (trim(name), iou, 3, (/id_xt,id_yt,id_zt/)
     &, vmin, vmax,' ', 'F', trim(long_name), trim(stand_name)
     &, trim(units))
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('zt', iou, jmt, (/1/), (/km/), zt, 1., 0.)
      call putvaramsk (trim(name), iou, imt*jmt*km, (/1,1,1/)
     &, (/imt,jmt,km/), var(:,:,:,1), vmask, 1., -273.15)
      call closefile (iou)
      
      do n=1,12

!=======================================================================
!     read monthly data
!=======================================================================

        call openfile (trim(path)//trim(name)//'_monthly.nc', iou)
        call getvara ('X', iou, id, (/1/), (/id/), xd, 1., 0.)
        call getvara ('Y', iou, jd, (/1/), (/jd/), yd, 1., 0.)
        call getvara (name2, iou, id*jd*kdm, (/1,1,1,n/)
     &,   (/id,jd,kdm,1/), data(:,:,1:kdm), 1., 0.)
        where (data(:,:,:) .le. missing_value) data(:,:,:) = 2e20

!=======================================================================
!     rotate and interpolate monthly model surface data
!=======================================================================

        k = 1
        ka = 1
        do m=2,kd
          if (zd(m) .lt. zt(k)) ka = m
        enddo
        kb = min(ka+1,kd)
        ka = max(kb-1,1)
        wtb = min(1.0,max(0.0,(zt(k) - zd(ka))/(zd(kb) - zd(ka))))
        call rot_intrp_sclr (data(:,:,ka), xd, yd, id, jd, tmp(:,:,1)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call rot_intrp_sclr (data(:,:,kb), xd, yd, id, jd, tmp(:,:,2)
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        call intrp_vert (tmp(:,:,1), tmp(:,:,2), wtb, -1e20
     &,   var(:,:,n,1), imt, jmt)
        call extrap (var(:,:,n,1), -1.e10, xt, maskv, imt, jmt, 1)

      enddo

!=======================================================================
!     calculate potential temperature from in situ
!=======================================================================

      inquire (file='../sss_mth.nc', exist=exists)
      if (exists) then
        call openfile ("../sss_mth.nc", iou)
        call getvara ('sss', iou, imt*jmt*12, (/1,1,1/)
     &,   (/imt,jmt,12/), var(:,:,1:12,2), 1., 0.)
        call closefile (iou)

        k=1
        do j=1,jmt
          do i=1,imt
            t = var(i,j,k,1)
            s = var(i,j,k,2)
            p = zt(k)
            call potem (t, s, p, theta0)
            var(i,j,k,1) = theta0
          enddo
        enddo
        long_name = "potential temperature"
        stand_name = "potential_temperature"
      else
        print*, 'Warning: sss_mth.nc missing, in situ temperature'
        long_name = "in situ temperature"
        stand_name = "in_situ_temperature"
      endif
       
!=======================================================================
!     write monthly netcdf surface data
!=======================================================================

      call opennew ("../sst_mth.nc", iou)
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
      write(name2, '(f8.2)') zt(1)
      long_name = trim(long_name)//' at '//trim(name2)//' m'
      stand_name = trim(stand_name)//' at '//trim(name2)//' m'
      call defvar ('sst', iou, 3, (/id_xt,id_yt,id_time/)
     &, vmin, vmax,' ', 'F', 'sea surface temperature'
     &, 'sea_surface_temperature', trim(units))
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        call putvaramsk ('sst', iou, imt*jmt, (/1,1,n/)
     &,   (/imt,jmt,1/), var(:,:,n,1), vmask(:,:,1), 1., -273.15)
      enddo
      call closefile (iou)

      end
      
      subroutine potem (t, s, p, theta)

!=======================================================================
!     this subroutine calculates potential temperature as a function
!     of in-situ temperature, salinity, and pressure.

!     input [units]:
!       in-situ temperature (t): [degrees centigrade]
!       salinity (s): [per mil]
!       pressure (p): [decibars, approx. as meters of depth]
!     output [units]:
!       potential temperature (theta): [degrees centigrade]

!     references:
!        based on Fofonoff and Froese (1958) as shown in ...
!        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
!          New York, 1962, page 17, table iv.

!-----------------------------------------------------------------------

      implicit real(kind=8) (a-h,o-z)

!=======================================================================

      b1    = -1.60d-5*p
      b2    = 1.014d-5*p*t
      t2    = t*t
      t3    = t2*t
      b3    = -1.27d-7*p*t2
      b4    = 2.7d-9*p*t3
      b5    = 1.322d-6*p*s
      b6    = -2.62d-8*p*s*t
      s2    = s*s
      p2    = p*p
      b7    = 4.1d-9*p*s2
      b8    = 9.14d-9*p2
      b9    = -2.77d-10*p2*t
      b10   = 9.5d-13*p2*t2
      b11   = -1.557d-13*p2*p
      potmp = b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11
      theta = t-potmp

      return
      end
