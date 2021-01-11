      program a_calb_nc

!=======================================================================
!     creates atmospheric coalbedo data file a_calb.nc
!     also creates a_alb.nc, s_alb.nc and p_alb.nc from data
!=======================================================================

      implicit none

      real, allocatable :: avar(:,:,:), var(:,:,:,:)
      real, allocatable :: xt(:), yt(:)

      integer id, jd
      parameter (id=144, jd=72)
      real data(id,jd,12,3), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, k, l, m, n, id_xt, id_yt, id_time, nf
      real psi, theta, phi, daymon(12), pass, ff
      real rad, radian, sj
      
      logical exists
      
      character(120) :: path
      
      data daymon(1:6)  / 15.5,  45.,   74.5, 105.,  135.5, 166. /
      data daymon(7:12) /196.5, 227.5, 258.,  288.5, 319.,  349.5/

      rad = acos(-1.)/180.
      radian = 360./(2.* 4.0*atan(1.0))
!     assuming absorbtion on incoming and refected, if just on incoming
!     use the sqrt(1-scatter) for pass
      pass = 1. - 0.23
      nf = 2
      ff = 1.015

!=======================================================================
!     set path and read path file if it exists
!=======================================================================

      path = '/usr/local/models/UVic_ESCM/data_source'
      inquire (file='../path', exist=exists)
      if (exists) then
        open (10,file='../path')
        read (10,'(a)') path
      endif

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( var(imt,jmt,12,3) )
      allocate ( avar(imt,jmt,3) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read data
!=======================================================================

      call openfile (trim(path)//'/ERBE/surface_albedo.nc', iou)
      call getvara ('X', iou, id, (/1/), (/id/), xd, 1., 0.)
      call getvara ('Y', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!     flip latitude
      yd(1:jd) = yd(jd:1:-1)
      call closefile (iou)

      do n=1,12
        call openfile (trim(path)//'/ERBE/total_albedo.nc', iou)
        call getvara ('albedo', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,n,2), 1., 0.)
!       flip data in latitude
        data(:,1:jd,n,2) = data(:,jd:1:-1,n,2)
        call closefile (iou)
      enddo
      
      data(:,:,:,2) = data(:,:,:,2)*ff

      do n=1,12
        call openfile (trim(path)//'/ERBE/surface_albedo.nc', iou)
        call getvara ('albedo', iou, id*jd, (/1,1,n/), (/id,jd,1/)
     &,   data(:,:,n,3), 1., 0.)
        data(:,1:jd,n,3) = data(:,jd:1:-1,n,3)
        call closefile (iou)
      enddo
      data(:,:,:,:) = data(:,:,:,:)/100.
      
      do k=2,3

!=======================================================================
!       get nearest albedo in time to get values for dark periods
!=======================================================================

        do n=1,12
          do j=1,jd
            do i=1,id
              if (data(i,j,n,k).gt.1. .or. data(i,j,n,k).lt.0.) then
                l = 100
                do m=1,12
                  if (data(i,j,m,k).le.1 .and. data(i,j,m,k).ge.0) then
                    if (abs(m-n) .lt. l) then
                      l = abs(m-n)
                      data(i,j,n,k) = data(i,j,m,k)
                    elseif (abs(m-n) .eq. l) then
                      data(i,j,n,k) = 0.5*(data(i,j,n,k)+data(i,j,m,k))
                    endif
                  endif
                enddo
              endif
              if (data(i,j,n,k).gt.1. .or. data(i,j,n,k).lt.0.) then
                data(i,j,n,k) = 2.e20
              endif
            enddo
          enddo
        enddo

!=======================================================================
!       rotate and interpolate data
!=======================================================================

        do n=1,12
          call rot_intrp_sclr (data(:,:,n,k), xd, yd, id, jd
     &,     var(:,:,n,k), xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
        enddo
          
      enddo
       
!=======================================================================
!     calculate atmospheric albedo data
!=======================================================================

      do n=1,12
        do j=1,jmt
          do i=1,imt
            var(i,j,n,1) = (var(i,j,n,2) - var(i,j,n,3)*pass*pass)
     &                     /(1. - var(i,j,n,3)*pass*pass)
          enddo
        enddo
      enddo

      if (nf .gt. 0) then
        print*, 'Warning: applying smoothing to atmospheric albedo'
        do n=1,12
          call smooth (var(:,:,n,1), imt, jmt, nf)
        enddo
      endif

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:,:,:) = var(imt-1,:,:,:)
      var(imt,:,:,:) = var(2,:,:,:)
      avar(:,:,:) = 0.
      do k=1,3
        do n=1,12
          avar(:,:,k) = avar(:,:,k) + var(:,:,n,k)
        enddo
      enddo
      avar(:,:,:) = avar(:,:,:)/12

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call opennew ("../p_alb_mth.nc", iou)
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
      call defvar ('p_alb', iou, 3, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'planetary albedo'
     &, 'planetary_albedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        call putvara ('p_alb', iou, imt*jmt, (/1,1,n/)
     &,              (/imt,jmt,1/), var(1,1,n,2), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../p_alb_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('p_alb', iou, 2, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'planetary albedo'
     &, 'planetary_coalbedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('p_alb', iou, imt*jmt, (/1,1/)
     &,              (/imt,jmt/), avar(1,1,2), 1., 0.)
      call closefile (iou)

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call opennew ("../a_alb_mth.nc", iou)
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
      call defvar ('a_alb', iou, 3, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'atmospheric albedo'
     &, 'atmospheric_albedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        call putvara ('a_alb', iou, imt*jmt, (/1,1,n/)
     &,              (/imt,jmt,1/), var(1,1,n,1), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../a_alb_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('a_alb', iou, 2, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'atmospheric albedo'
     &, 'atmospheric_albedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('a_alb', iou, imt*jmt, (/1,1/)
     &,              (/imt,jmt/), avar(1,1,1), 1., 0.)
      call closefile (iou)


!=======================================================================
!     write monthly netcdf data
!=======================================================================

      call opennew ("../s_alb_mth.nc", iou)
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
      call defvar ('s_alb', iou, 3, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'surface albedo'
     &, 'surface_albedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        call putvara ('s_alb', iou, imt*jmt, (/1,1,n/)
     &,              (/imt,jmt,1/), var(1,1,n,3), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      call opennew ("../s_alb_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('s_alb', iou, 2, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'surface albedo'
     &, 'surface_albedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('s_alb', iou, imt*jmt, (/1,1/)
     &,              (/imt,jmt/), avar(1,1,3), 1., 0.)
      call closefile (iou)

!=======================================================================
!     write monthly netcdf data
!=======================================================================

      var(:,:,:,1) = 1. - var(:,:,:,1)
      
      call opennew ("../a_calb_mth.nc", iou)
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
      call defvar ('a_calb', iou, 3, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'atmospheric coalbedo'
     &, 'atmospheric_coalbedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      do n=1,12
        call putvars ('time', iou, n, daymon(n)/365., 1., 0.)
        call putvara ('a_calb', iou, imt*jmt, (/1,1,n/)
     &,              (/imt,jmt,1/), var(1,1,n,1), 1., 0.)
      enddo
      call closefile (iou)

!=======================================================================
!     write annual netcdf data
!=======================================================================

      avar(:,:,1) = 1. - avar(:,:,1)

      call opennew ("../a_calb_ann.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('a_calb', iou, 2, (/id_xt,id_yt,id_time/)
     &, 0., 1., ' ', 'F', 'atmospheric coalbedo'
     &, 'atmospheric_coalbedo', '1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('a_calb', iou, imt*jmt, (/1,1/)
     &,              (/imt,jmt/), avar(1,1,1), 1., 0.)
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
