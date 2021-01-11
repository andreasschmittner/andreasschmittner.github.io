      program kmt_nc

!=======================================================================
!     creates tracer grid ocean depth mask data file kmt.nc
!=======================================================================

      implicit none

      character(1), allocatable :: line(:)
      real, allocatable :: var(:,:), xt(:), yt(:), zw(:)

      integer id, jd
      parameter (id=360, jd=180)
      real data(id,jd), xd(id), yd(jd)

      integer i, imt, iou, j, jmt, k, km, n, id_xt, id_yt
      real psi, theta, phi, total, dz
      
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
      path = trim(path)//'/topog/'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      call getdimlen ('zw', iou, km)
      allocate ( xt(imt) )
      allocate ( line(imt) )
      allocate ( yt(jmt) )
      allocate ( var(imt,jmt) )
      allocate ( zw(km) )
      call getvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call getvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call getvara ('zw', iou, km, (/1/), (/km/), zw, 1., 0.)
      call getvars ('psi', iou, 1, psi, 1., 0.)
      call getvars ('theta', iou, 1, theta, 1., 0.)
      call getvars ('phi', iou, 1, phi, 1., 0.)
      call closefile (iou)

      inquire (file='../kmt.map', exist=exists)
      
      if (exists) then

!=======================================================================
!       read kmt map (if map exists assume it is corrected)
!=======================================================================

        print*, 'Warning: reading kmt from kmt.map'
      
        open (10, file='../kmt.map', status='UNKNOWN')
        do j=jmt,1,-1
          read (10,'(1000a)') (line(i),i=1,imt)
          do i=1,imt
            n = ichar(line(i))
            if (n .ge. 65 .and. n .le. 90) then
              n = n - 65 + 1
            elseif (n .ge. 97 .and. n .le. 122) then
              n = n - 97 + 27
            elseif (n .ge. 49 .and. n .le. 57) then
              n = n -49 + 53
            else
              n = 0
            endif
            var(i,j) = float(n)
          enddo
        enddo
        close (10)

        call fill_bays (var, imt, jmt)

      else

!=======================================================================
!       read elevation data
!=======================================================================

        call openfile (trim(path)//'elev.1-deg.nc', iou)
        call getvara ('lon', iou, id, (/1/), (/id/), xd, 1., 0.)
        call getvara ('lat', iou, jd, (/1/), (/jd/), yd, 1., 0.)
!       flip latitude
        yd(1:jd) = yd(jd:1:-1)
        call getvara ('data', iou, id*jd, (/1,1,1/), (/id,jd,1/)
     &,   data(:,:), 1., 0.)
!       flip data in latitude
        data(:,1:jd) = data(:,jd:1:-1)
        call closefile (iou)
        
!=======================================================================
!       rotate and interpolate data
!=======================================================================

        call rot_intrp_sclr (data, xd, yd, id, jd, var
     &,   xt, yt, imt, jmt, phi, theta, psi, -1.e20, 0)
     

!=======================================================================
!       smooth data
!=======================================================================

!        call smooth (var, imt, jmt, 0)

!=======================================================================
!       set to kmt to closest model level
!=======================================================================

        var(:,:) = -var(:,:)
        do j=1,jmt
          do i=1,imt
            dz = zw(km)
            n = 0
            do k=1,km
             if (abs(var(i,j) - zw(k)) .lt. dz) then
               dz = abs(var(i,j) - zw(k))
                n = k
              endif
            enddo
            if (var(i,j) .gt. 0) then
              var(i,j) = n
            else
              var(i,j) = 0
            endif          
          enddo
        enddo

        call fill_bays (var, imt, jmt)

!=======================================================================
!       write kmt map (if map does not exist)
!=======================================================================

        open (10, file='../kmt.map', status='UNKNOWN')
        do j=jmt,1,-1
          do i=1,imt
            n = int(var(i,j))
            if (n .ge. 1 .and. n .le. 26) then
              n = n + 65 - 1
            elseif (n .ge. 27 .and. n .le. 52) then
              n = n + 97 - 27
            elseif (n .ge. 53 .and. n .le. 61) then
              n = n + 49 - 53
            else
              n = ichar(".")
            endif
            line(i) = char(n)
          enddo
          write (10,'(1000a)') (line(i),i=1,imt)
        enddo
        close (10)

      endif

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      var(1,:) = var(imt-1,:)
      var(imt,:) = var(2,:)

!=======================================================================
!     write netcdf kmt data
!=======================================================================

      call opennew ("../kmt.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('kmt', iou, 2, (/id_xt,id_yt/), 0., 1000.
     &, ' ', 'I', 'ocean depth grid level', ' ' ,'1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('kmt', iou, imt*jmt, (/1,1/)
     &,            (/imt,jmt/), var, 1., 0.)
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

      print*, 'Warning: filling isolated bays'

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
