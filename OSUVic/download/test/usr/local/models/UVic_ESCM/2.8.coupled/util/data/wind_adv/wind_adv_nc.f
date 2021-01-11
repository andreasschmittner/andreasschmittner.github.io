      program wind_adv_nc

!=======================================================================
!     creates wind data files wind_adv_mth.nc and wind_adv_ann.nc
!=======================================================================

      implicit none

      real, allocatable :: avar(:,:,:), var(:,:,:,:)
      real, allocatable :: xu(:), yu(:)

      integer id, jd, kd
      parameter (id=144, jd=73, kd=17)
      real data(id,jd,kd,5), xd(id), yd(jd)
      real z(kd), dz(kd), twt

      integer i, imt, iou, j, jmt, k, m, n, id_xu, id_yu, id_time
      real psi, theta, phi, wt, daymon(12), fwindq, fwindt, nf

      logical exists
      
      character(120) :: path
      
      data daymon(1:6)  / 15.5,  45.,   74.5, 105.,  135.5, 166. /
      data daymon(7:12) /196.5, 227.5, 258.,  288.5, 319.,  349.5/
      fwindq = 0.5
      fwindt = 0.1
      nf = 8
      
!=======================================================================
!     set path and read path file if it exists
!=======================================================================

      path = '/usr/local/models/UVic_ESCM/data_source'
      inquire (file='../path', exist=exists)
      if (exists) then
        open (10,file='../path')
        read (10,'(a)') path
      endif
      path = trim(path)//'/ncep/pressure/'
      
!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xu', iou, imt)
      call getdimlen ('yu', iou, jmt)
      allocate ( xu(imt) )
      allocate ( yu(jmt) )
      allocate ( var(imt,jmt,2,12) )
      allocate ( avar(imt,jmt,2) )
      call getvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call getvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read data grid
!=======================================================================

      call openfile (trim(path)//'uwnd.mon.ltm.nc', iou)
      call getvara ('lon', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('lat', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!     flip latitude
      yd(1:jd) = yd(jd:1:-1)
      call closefile (iou)

      avar(:,:,:) = 0.0
      do n=1,12

!=======================================================================
!       read data
!=======================================================================

        call openfile (trim(path)//'uwnd.mon.ltm.nc', iou)
        call getvara ('uwnd', iou, id*jd*kd, (/1,1,1,n/), (/id,jd,kd,1/)
     &,   data(:,:,:,1), 1., 0.)
!       flip data in latitude
        data(:,1:jd,:,1) = data(:,jd:1:-1,:,1)
        call closefile (iou)

        call openfile (trim(path)//'vwnd.mon.ltm.nc', iou)
        call getvara ('vwnd', iou, id*jd*kd, (/1,1,1,n/), (/id,jd,kd,1/)
     &,   data(:,:,:,2), 1., 0.)
!       flip data in latitude
        data(:,1:jd,:,2) = data(:,jd:1:-1,:,2)
        call closefile (iou)

        call openfile (trim(path)//'hgt.mon.ltm.nc', iou)
        call getvara ('hgt', iou, id*jd*kd, (/1,1,1,n/), (/id,jd,kd,1/)
     &,   data(:,:,:,3), 1., 0.)
!       flip data in latitude
        data(:,1:jd,:,3) = data(:,jd:1:-1,:,3)
        call closefile (iou)

!       only 8 levels of humidity instead of 17
        call openfile (trim(path)//'shum.mon.ltm.nc', iou)
        call getvara ('shum', iou, id*jd*8, (/1,1,1,n/), (/id,jd,8,1/)
     &,   data(:,:,1:8,4), 1., 0.)
        data(:,:,9:kd,4) = 0.
!       flip data in latitude
        data(:,1:jd,:,4) = data(:,jd:1:-1,:,4)
        call closefile (iou)

!=======================================================================
!       calculate advecting wind from weighted data
!=======================================================================

        do j=1,jd
          do i=1,id

            dz(1) = (data(i,j,2,3) - data(i,j,1,3))/2. + data(i,j,1,3)
            do k=2,kd-1
              dz(k) = (data(i,j,k+1,3) - data(i,j,k,3))/2.
     &              + (data(i,j,k,3) - data(i,j,k-1,3))/2.
            enddo
            dz(kd) = (data(i,j,kd,3) - data(i,j,kd-1,3))/2.

            twt = data(i,j,kd,4)*dz(kd)
            data(i,j,kd,1) = data(i,j,kd,1)*data(i,j,kd,4)*dz(kd)
            do k=1,kd-1
              data(i,j,kd,1) = data(i,j,kd,1)
     &                       + data(i,j,k,1)*data(i,j,k,4)*dz(k)
              twt = twt + data(i,j,k,4)*dz(k)
            enddo
            if (twt .ne. 0.) data(i,j,kd,1) = data(i,j,kd,1)/twt

            twt = data(i,j,kd,4)*dz(kd)
            data(i,j,kd,2) = data(i,j,kd,2)*data(i,j,kd,4)*dz(kd)
            do k=1,kd-1
              data(i,j,kd,2) = data(i,j,kd,2)
     &                       + data(i,j,k,2)*data(i,j,k,4)*dz(k)
              twt = twt + data(i,j,k,4)*dz(k)
            enddo
            if (twt .ne. 0.) data(i,j,kd,2) = data(i,j,kd,2)/twt

          enddo
        enddo

!=======================================================================
!       set wt to 1. for surface winds 0. for height integrated
!=======================================================================

        wt = 0.
        data(:,:,1,1) = data(:,:,1,1)*wt + data(:,:,kd,1)*(1.-wt)
        data(:,:,1,2) = data(:,:,1,2)*wt + data(:,:,kd,2)*(1.-wt)

!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_vctr (data(:,:,1,1:2), xd, yd, id, jd
     &,   var(:,:,1:2,n), xu, yu, imt, jmt, phi, theta, psi, -1.e20, 0)

!=======================================================================
!       set cyclic boundary condition
!=======================================================================

        var(1,:,:,n) = var(imt-1,:,:,n)
        var(imt,:,:,n) = var(2,:,:,n)

        avar(:,:,:) = avar(:,:,:) + var(:,:,:,n)
      enddo
      avar(:,:,:) = avar(:,:,:)/12.


      if (fwindq .ne. 1.) then
        print*, 'Warning: multiplying wind by factor: ', fwindq
        do i=1,imt
          do j=1,jmt
            avar(i,j,:) = avar(i,j,:)*fwindq
            var(i,j,:,:) = var(i,j,:,:)*fwindq
          enddo
        enddo
      endif
      
      
      if (nf .ge. 1.) then
        print*, 'Warning: Soothing wind'
        do m=1,nf
          do n=1,12
            call filter (var(1,1,1,n), imt, jmt)
            call filter (var(1,1,2,n), imt, jmt)
          enddo
          call filter (avar(1,1,1), imt, jmt)
          call filter (avar(1,1,2), imt, jmt)
        enddo
      endif

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call opennew ("../wind_adv_mth.nc", iou)
      call redef (iou)
      call defdim ('time', iou, 0, id_time)
      call defdim ('xu', iou, imt, id_xu)
      call defdim ('yu', iou, jmt, id_yu)
      call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'
     &, 'time', 'time', 'common_year since 1-1-1 00:00:0.0')
      call putatttext (iou, 'time', 'calendar', '365_day')
      call defvar ('xu', iou, 1, (/id_xu/), 0., 0., 'X', 'F'
     &, 'longitude of the u grid', 'longitude', 'degrees_east')
      call defvar ('yu', iou, 1, (/id_yu/), 0., 0., 'Y', 'F'
     &, 'latitude of the u grid', 'latitude', 'degrees_north')
      call defvar ('wx_q', iou, 3, (/id_xu,id_yu,id_time/), -1000.
     &,  1000., ' ', 'F', 'eastward wind for advection of humidity'
     &, 'eastward_wind', 'm s-1')
      call defvar ('wy_q', iou, 3, (/id_xu,id_yu,id_time/), -1000.
     &,  1000., ' ', 'F', 'northward wind for advection of humidity'
     &, 'northward_wind', 'm s-1')
      call defvar ('wx_t', iou, 3, (/id_xu,id_yu,id_time/), -1000.
     &,  1000., ' ', 'F', 'eastward wind for advection of temperature'
     &, 'eastward_wind', 'm s-1')
      call defvar ('wy_t', iou, 3, (/id_xu,id_yu,id_time/), -1000.
     &,  1000., ' ', 'F', 'northward wind for advection of temperature'
     &, 'northward_wind', 'm s-1')
      call enddef (iou)
      call putvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call putvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        call putvara ('wx_q', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(1,1,1,n), 1., 0.)
        call putvara ('wy_q', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(:,:,2,n), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../wind_adv_ann.nc", iou)
      call redef (iou)
      call defdim ('xu', iou, imt, id_xu)
      call defdim ('yu', iou, jmt, id_yu)
      call defvar ('xu', iou, 1, (/id_xu/), 0., 0., 'X', 'F'
     &, 'longitude of the u grid', 'longitude', 'degrees_east')
      call defvar ('yu', iou, 1, (/id_yu/), 0., 0., 'Y', 'F'
     &, 'latitude of the u grid', 'latitude', 'degrees_north')
      call defvar ('wx_q', iou, 2, (/id_xu,id_yu/), -1000.
     &, 1000., ' ', 'F', 'eastward wind for advection of humidity'
     &, 'eastward_wind', 'm s-1')
      call defvar ('wy_q', iou, 2, (/id_xu,id_yu/), -1000.
     &, 1000., ' ', 'F', 'northward wind for advection of humidity'
     &, 'northward_wind', 'm s-1')
      call defvar ('wx_t', iou, 2, (/id_xu,id_yu/), -1000.
     &, 1000., ' ', 'F', 'eastward wind for advection of temperature'
     &, 'eastward_wind', 'm s-1')
      call defvar ('wy_t', iou, 2, (/id_xu,id_yu/), -1000.
     &, 1000., ' ', 'F', 'northward wind for advection of temperature'
     &, 'northward_wind', 'm s-1')
      call enddef (iou)
      call putvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call putvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      call putvara ('wx_q', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar(:,:,1), 1., 0.)
      call putvara ('wy_q', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar(:,:,2), 1., 0.)
      call closefile (iou)

      avar(:,:,:) = 0.0
      do n=1,12

!=======================================================================
!       read data
!=======================================================================

        call openfile (trim(path)//'uwnd.mon.ltm.nc', iou)
        call getvara ('uwnd', iou, id*jd*kd, (/1,1,1,n/), (/id,jd,kd,1/)
     &,   data(:,:,:,1), 1., 0.)
!       flip data in latitude
        data(:,1:jd,:,1) = data(:,jd:1:-1,:,1)
        call closefile (iou)

        call openfile (trim(path)//'vwnd.mon.ltm.nc', iou)
        call getvara ('vwnd', iou, id*jd*kd, (/1,1,1,n/), (/id,jd,kd,1/)
     &,   data(:,:,:,2), 1., 0.)
!       flip data in latitude
        data(:,1:jd,:,2) = data(:,jd:1:-1,:,2)
        call closefile (iou)

        call openfile (trim(path)//'hgt.mon.ltm.nc', iou)
        call getvara ('hgt', iou, id*jd*kd, (/1,1,1,n/), (/id,jd,kd,1/)
     &,   data(:,:,:,3), 1., 0.)
!       flip data in latitude
        data(:,1:jd,:,3) = data(:,jd:1:-1,:,3)
        call closefile (iou)

        call openfile (trim(path)//'air.mon.ltm.nc', iou)
        call getvara ('air', iou, id*jd*kd, (/1,1,1,n/), (/id,jd,kd,1/)
     &,   data(:,:,:,4), 1., 0.)
        data(:,:,:,4) = data(:,:,:,4) + 273.15
!       flip data in latitude
        data(:,1:jd,:,4) = data(:,jd:1:-1,:,4)
        call closefile (iou)

!=======================================================================
!       calculate advecting wind from weighted data
!=======================================================================

        do j=1,jd
          do i=1,id

            dz(1) = (data(i,j,2,3) - data(i,j,1,3))/2. + data(i,j,1,3)
            do k=2,kd-1
              dz(k) = (data(i,j,k+1,3) - data(i,j,k,3))/2.
     &              + (data(i,j,k,3) - data(i,j,k-1,3))/2.
            enddo
            dz(kd) = (data(i,j,kd,3) - data(i,j,kd-1,3))/2.

            twt = data(i,j,kd,4)*dz(kd)
            data(i,j,kd,1) = data(i,j,kd,1)*data(i,j,kd,4)*dz(kd)
            do k=1,kd-1
              data(i,j,kd,1) = data(i,j,kd,1)
     &                       + data(i,j,k,1)*data(i,j,k,4)*dz(k)
              twt = twt + data(i,j,k,4)*dz(k)
            enddo
            if (twt .ne. 0.) data(i,j,kd,1) = data(i,j,kd,1)/twt

            twt = data(i,j,kd,4)*dz(kd)
            data(i,j,kd,2) = data(i,j,kd,2)*data(i,j,kd,4)*dz(kd)
            do k=1,kd-1
              data(i,j,kd,2) = data(i,j,kd,2)
     &                       + data(i,j,k,2)*data(i,j,k,4)*dz(k)
              twt = twt + data(i,j,k,4)*dz(k)
            enddo
            if (twt .ne. 0.) data(i,j,kd,2) = data(i,j,kd,2)/twt

          enddo
        enddo

!=======================================================================
!       set wt to 1. for surface winds 0. for height integrated
!=======================================================================

        wt = 0.
        data(:,:,1,1) = data(:,:,1,1)*wt + data(:,:,kd,1)*(1.-wt)
        data(:,:,1,2) = data(:,:,1,2)*wt + data(:,:,kd,2)*(1.-wt)

!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_vctr (data(:,:,1,1:2), xd, yd, id, jd
     &,   var(:,:,1:2,n), xu, yu, imt, jmt, phi, theta, psi, -1.e20, 0)

!=======================================================================
!       set cyclic boundary condition
!=======================================================================

        var(1,:,:,n) = var(imt-1,:,:,n)
        var(imt,:,:,n) = var(2,:,:,n)

        avar(:,:,:) = avar(:,:,:) + var(:,:,:,n)
      enddo
      avar(:,:,:) = avar(:,:,:)/12.


      if (fwindt .ne. 1.) then
        print*, 'Warning: multiplying wind by factor: ', fwindt
        do i=1,imt
          do j=1,jmt
            avar(i,j,:) = avar(i,j,:)*fwindt
            var(i,j,:,:) = var(i,j,:,:)*fwindt
          enddo
        enddo
      endif
      
      
      if (nf .ge. 1.) then
        print*, 'Warning: Soothing wind'
        do m=1,nf
          do n=1,12
            call filter (var(1,1,1,n), imt, jmt)
            call filter (var(1,1,2,n), imt, jmt)
          enddo
          call filter (avar(1,1,1), imt, jmt)
          call filter (avar(1,1,2), imt, jmt)
        enddo
      endif

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call openfile ("../wind_adv_mth.nc", iou)
      do n=1,12
        call putvara ('wx_t', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(:,:,1,n), 1., 0.)
        call putvara ('wy_t', iou, imt*jmt, (/1,1,n/), (/imt,jmt,1/)
     &,   var(:,:,2,n), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call openfile ("../wind_adv_ann.nc", iou)
      call putvara ('wx_t', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar(:,:,1), 1., 0.)
      call putvara ('wy_t', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, avar(:,:,2), 1., 0.)
      call closefile (iou)

      end

      subroutine filter (data, imt, jmt)

      implicit none

      integer i, j, imt, jmt
      real data(imt,jmt), tmp(imt,jmt), wt

      data(1,:) = data(imt-1,:)
      data(imt,:) = data(2,:)
      data(:,1) = data(:,2)
      data(:,jmt) = data(:,jmt-1)
      tmp(:,:) = data(:,:)
      wt = 0.5
      do j=2,jmt-1
        do i=2,imt-1
          data(i,j) = wt*tmp(i,j) 
     &      +  (tmp(i-1,j+1) + tmp(i,j+1) + tmp(i+1,j+1)  
     &      +  tmp(i-1,j)                 + tmp(i+1,j)
     &      +  tmp(i-1,j-1)  + tmp(i,j-1) + tmp(i+1,j-1))*wt/8.
        enddo
      enddo
      data(1,:) = data(imt-1,:)
      data(imt,:) = data(2,:)
      data(:,1) = data(:,2)
      data(:,jmt) = data(:,jmt-1)

      return
      end       
