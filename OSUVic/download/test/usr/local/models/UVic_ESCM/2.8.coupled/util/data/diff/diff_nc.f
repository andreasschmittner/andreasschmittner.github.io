      program diff_nc

!=======================================================================
!     creates diffusion file diff.nc
!=======================================================================

      implicit none

      integer nat
      parameter (nat=2)

      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:,:,:), xt(:), yt(:), xu(:), yu(:)
      real, allocatable :: vtmp(:,:,:)

      integer id, jd
      parameter (id=361, jd=181)
      real data(id,jd,2,nat,2), xd(id), yd(jd), dtmp(id,jd,2)

      integer i, imt, iou, j, jmt, k, m, n, nf
      integer id_xt, id_yt, id_xu, id_yu
      real psi, theta, phi, rad, syd, s2yd, fy

      real (kind=8) :: dphir, dthetar, dpsir, rlt, rln, a
      real (kind=8) :: x, y, temp1, temp2
      
      character(3) :: a3

      rad = acos(-1.)/180.
      nf = 0

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( yt(jmt) )
      allocate ( xu(imt) )
      allocate ( yu(jmt) )
      allocate ( var(imt,jmt,2,2) )
      allocate ( vtmp(imt,jmt,2) )
      allocate ( maskv(imt,jmt) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call getvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

!=======================================================================
!     read kmt data
!=======================================================================

      call openfile ("../kmt.nc", iou)
      call getvara ('kmt', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &, var(:,:,1,1), 1., 0.)
      maskv(:,:) = nint(var(:,:,1,1))
      call closefile (iou)

!=======================================================================
!     set diffusion on a "data" grid
!=======================================================================

      data(:,:,:,:,:) = 1.e10
      do j=1,jd
        yd(j) = float(j-1) - 90.
        syd = abs(sin(yd(j)*rad))
        s2yd = abs(sin(2.*yd(j)*rad))
        fy = abs(yd(j))/90.
        do i=1,id
          xd(i) = float(i-1)

!         diffusion north of the Equator
          if (yd(j) .ge. 0.) then

!           heat e-w
            data(i,j,1,1,1) = 3.3e10*(1.4 - 1.1*syd**2)
            data(i,j,1,1,2) = 3.3e10*(1.4 - 1.1*syd**2)

!           heat n-s
            data(i,j,2,1,1) = 3.3e10*(1.4 - 1.1*syd**2)
            data(i,j,2,1,2) = 3.3e10*(1.4 - 1.1*syd**2)

!           moisture e-w
            data(i,j,1,2,1) = 2.e9 + 1.5e10*fy**3 + 3.e11*s2yd**6
            data(i,j,1,2,2) = 2.e9 + 1.5e10*fy**3 + 3.e11*s2yd**6

!           moisture n-s
            data(i,j,2,2,1) = 7.e9 + 0.e9*fy**3 + 5.e9*s2yd**4
            data(i,j,2,2,2) = 7.e9 + 0.e9*fy**3 + 5.e9*s2yd**4

          endif
            
!         diffusion south of the Equator
          if (yd(j) .lt. 0.) then

!           heat e-w
            data(i,j,1,1,1) = 3.3e10*(1.4 - 1.3*syd**2)
            data(i,j,1,1,2) = 3.3e10*(1.4 - 1.3*syd**2)
            
!           heat n-s
            data(i,j,2,1,1) = 3.3e10*(1.4 - 1.3*syd**2)
            data(i,j,2,1,2) = 3.3e10*(1.4 - 1.3*syd**2)

!           moisture e-w over ocean
            data(i,j,1,2,1) = 2.e9 + 1.5e10*fy**3 + 3.e11*s2yd**6
            data(i,j,1,2,2) = 2.e9 + 1.5e10*fy**3 + 3.e11*s2yd**6

!           moisture n-s over ocean
            data(i,j,2,2,1) = 7.e9 + 1.5e10*fy**3 + 2.5e10*s2yd**4
            data(i,j,2,2,2) = 7.e9 + 1.5e10*fy**3 + 2.5e10*s2yd**4

          endif

        enddo
      enddo

!=======================================================================
!     define diffusion file
!=======================================================================

      call opennew ("../diff.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defdim ('xu', iou, imt, id_xu)
      call defdim ('yu', iou, jmt, id_yu)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('xu', iou, 1, (/id_xu/), 0., 0., 'X', 'F'
     &, 'longitude of the u grid', 'longitude', 'degrees_east')
      call defvar ('yu', iou, 1, (/id_yu/), 0., 0., 'Y', 'F'
     &, 'latitude of the u grid', 'latitude', 'degrees_north')
      do n=1,nat
        if (n .lt. 1000) write(a3,'(i3)') n
        if (n .lt. 100) write(a3,'(i2)') n
        if (n .lt. 10) write(a3,'(i1)') n
        call defvar ('dn_'//trim(a3), iou ,2, (/id_xt,id_yu/), 0., 1.e15
     &, ' ', 'F', 'northward diffusion for tracer '//trim(a3), '', '')
        call defvar ('de_'//trim(a3), iou ,2, (/id_xu,id_yt/), 0., 1.e15
     &, ' ', 'F', 'eastward diffusion for tracer '//trim(a3), '', '')
      enddo
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('xu', iou, imt, (/1/), (/imt/), xu, 1., 0.)
      call putvara ('yu', iou, jmt, (/1/), (/jmt/), yu, 1., 0.)
      
      do n=1,nat
        if (n .lt. 1000) write(a3,'(i3)') n
        if (n .lt. 100) write(a3,'(i2)') n
        if (n .lt. 10) write(a3,'(i1)') n

!=======================================================================
!       rotate and interpolate diffusion
!=======================================================================

        var(:,:,:,:) = 0.        


        do m=1,2

          if (n .eq. 2) then
          
!           use unrotated components for humidity
            call rot_intrp_sclr (data(1,1,1,n,m), xd, yd, id, jd
     &,                          var(1,1,1,m), xt, yt, imt, jmt
     &,                          phi, theta, psi, -1.e20, 0)
            call rot_intrp_sclr (data(1,1,2,n,m), xd, yd, id, jd
     &,                          var(1,1,2,m), xt, yt, imt, jmt
     &,                          phi, theta, psi, -1.e20, 0)
          else

!           vector method
            do k=1,4
              dtmp(:,:,:) = data(:,:,1:2,n,m)
              if (k .eq. 2 .or. k .eq. 4) dtmp(:,:,1) = -dtmp(:,:,1)
              if (k .eq. 3 .or. k .eq. 4) dtmp(:,:,2) = -dtmp(:,:,2)
              call rot_intrp_vctr (dtmp, xd, yd, id, jd, vtmp, xu, yt
     &,         imt, jmt, phi, theta, psi, -1.e20, 0)
              var(:,:,:,m) = var(:,:,:,m) + 0.25*abs(vtmp(:,:,:))
            enddo
            vtmp(:,:,:) = var(:,:,:,m)

!           ellipse method
!            call rot_intrp_sclr (data(1,1,1,n,m), xd, yd, id, jd
!     &,                          var(1,1,1,m), xt, yt, imt, jmt
!     &,                          phi, theta, psi, -1.e20, 0)
!            call rot_intrp_sclr (data(1,1,2,n,m), xd, yd, id, jd
!     &,                          var(1,1,2,m), xt, yt, imt, jmt
!     &,                          phi, theta, psi, -1.e20, 0)
!            dphir = phi
!            dthetar = theta
!            dpsir = psi
!            do j=1,jmt
!              do i=1,imt
!                rlt = yt(j)
!                rln = xt(i)
!                call drot_angle (rlt, rln, dphir, dthetar, dpsir, a)
!                x = var(i,j,1,m)! + 0.1*var(i,j,2,m)
!                y = var(i,j,2,m)! + 0.1*var(i,j,1,m)
!                temp1 = (x**2)*(y**2)/(y**2 + (x**2)*(tan(a))**2)
!                temp2 = (y**2)*(1-(y**2)/(y**2 + (x**2)*((tan(a))**2)))
!                var(i,j,1,m) = sqrt(temp1 + temp2)
!                a = a + 90.*rad
!                temp1 = (x**2)*(y**2)/(y**2 + (x**2)*(tan(a))**2)
!                temp2 = (y**2)*(1-(y**2)/(y**2 + (x**2)*((tan(a))**2)))
!                var(i,j,2,m) = sqrt(temp1 + temp2)
!              enddo
!            enddo

          endif
        enddo
        
        do j=1,jmt-1
          do i=1,imt-1
!           use land diffusion over land points
            if (maskv(i,j) .eq. 0) var(i,j,1,1) = var(i,j,1,2)
            if (maskv(i+1,j) .eq. 0) var(i,j,1,1) = var(i,j,1,2)
            if (maskv(i,j) .eq. 0) var(i,j,2,1) = var(i,j,2,2)
            if (maskv(i,j+1) .eq. 0) var(i,j,2,1) = var(i,j,2,2)
          enddo
        enddo
        
        if (nf .ge. 1.) then
          print*, 'Warning: Soothing diffusion'
          do m=1,nf
            call filter (var(1,1,1,1), imt, jmt)
            call filter (var(1,1,2,1), imt, jmt)
          enddo
        endif
        
        var(1,:,:,:) = var(imt-1,:,:,:)
        var(imt,:,:,:) = var(2,:,:,:)
        var(:,1,:,:) = var(:,2,:,:)
        var(:,jmt,:,:) = var(:,jmt-1,:,:)

!=======================================================================
!       write diffusion
!=======================================================================

        call putvara ('de_'//trim(a3), iou, imt*jmt, (/1,1/)
     &, (/imt,jmt/), var(:,:,1,1), 1., 0.)
        call putvara ('dn_'//trim(a3), iou, imt*jmt, (/1,1/)
     &, (/imt,jmt/), var(:,:,2,1), 1., 0.)

      enddo
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
