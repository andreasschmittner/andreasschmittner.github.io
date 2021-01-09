      program rivers_nc

!=======================================================================
!     creates river data file rivers.nc
!=======================================================================

      implicit none

      character(1), allocatable :: line(:)
      integer, allocatable :: maskv(:,:)
      real, allocatable :: var(:,:,:), xt(:), yt(:)

      integer id, jd
      parameter (id=102, jd=102)
      real data(id,jd,3), xd(id), yd(jd)

      integer i, imt, ios, iou, j, jmt, k, m, n, nb, ndp
      integer id_xt, id_yt

      real psi, theta, phi, wt

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
!      path = trim(path)//'/rivers/'
!     use local rivers.nc file
      path = '.'

!=======================================================================
!     read grid data
!=======================================================================

      call openfile ("../grid.nc", iou)
      call getdimlen ('xt', iou, imt)
      call getdimlen ('yt', iou, jmt)
      allocate ( xt(imt) )
      allocate ( line(imt) )
      allocate ( yt(jmt) )
      allocate ( var(imt,jmt,3) )
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
     &, var(:,:,1), 1., 0.)
      maskv = var(:,:,1)

!=======================================================================
!     interpolate from another rivers.nc file if it exists
!=======================================================================

      var(:,:,:) = 0.
      where (maskv(:,:) .eq. 0) var(:,:,1) = 1.
      
      inquire (file=trim(path)//'/rivers.nc', exist=exists)
      if (exists) then
        call openfile (trim(path)//'/rivers.nc', iou)
        call getvara ('xt', iou, id, 1, id, xd, 1., 0.)
        call getvara ('yt', iou, jd, 1, jd, yd, 1., 0.)
        call getvara ('rivers', iou, id*jd, (/1,1/), (/id,jd/)
     &,   data(:,:,1), 1., 0.)
        call getvara ('discharge', iou, id*jd, (/1,1/), (/id,jd/)
     &,   data(:,:,2), 1., 0.)
        call getvara ('weights', iou, id*jd, (/1,1/), (/id,jd/)
     &,   data(:,:,3), 1., 0.)
        call closefile (iou)
      
        do n=1,3
          call rot_intrp_sclr (data(:,:,n), xd, yd, id, jd, var(:,:,n)
     &,     xt, yt, imt, jmt, phi, theta, psi, -1.e20, 1)
        enddo
      endif

!=======================================================================
!     set cyclic boundary condition
!=======================================================================

      data(1,:,:) = data(imt-1,:,:)
      data(imt,:,:) = data(2,:,:)

!=======================================================================
!     read existing rivers.map file
!=======================================================================

      inquire (file='../rivers.map', exist=exists) 
      if (exists) then

!=======================================================================
!       read rivers map (if map exists assume it is corrected)
!=======================================================================

        print*, 'Warning: reading rivers from rivers.map'
      
        open (10, file='../rivers.map', status='UNKNOWN')
        read (10,*)
        do j=jmt-1,2,-1
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
            var(i,j,1) = float(n)
          enddo
        enddo
        read (10,*)
        ios = 0
        do while (ios .eq. 0) 
          read (10,'(a1,i3,F7.3)',iostat=ios) line(1), ndp, wt
          if (ios .eq. 0) then
            n = ichar(line(1))
            if (n .ge. 65 .and. n .le. 90) then
              n = n - 65 + 1
            elseif (n .ge. 97 .and. n .le. 122) then
              n = n - 97 + 27
            elseif (n .ge. 49 .and. n .le. 57) then
              n = n -49 + 53
            else
              n = 0
            endif
            do k=1,ndp
              read (10,'(2i4)') i, j
              var(i,j,2) = n
              var(i,j,3) = wt/float(ndp)
            enddo
          endif
        enddo
        close (10)
      
      else

!=======================================================================
!       write rivers map (if map does not exist)
!=======================================================================

        open (10, file='../rivers.tmp', status='UNKNOWN')
        write (10,*) ' '
        do j=jmt-1,2,-1
          do i=1,imt
            n = int(var(i,j,1))
            if (n .ge. 1 .and. n .le. 26) then
              n = n + 65 - 1
            elseif (n .ge. 27 .and. n .le. 52) then
              n = n + 97 - 27
            elseif (n .ge. 53 .and. n .le. 61) then
              n = n + 49 - 53
            else
              n = 0
            endif
            line(i) = char(n)
            if (n .eq. 0) line(i) = "#"
            if (maskv(i,j) .gt. 0 .and. n .ne. 0) line(i) = "#"
            if (maskv(i,j) .eq. 0 .and. n .eq. 0) line(i) = "."
             m = int(var(i,j,2))        
            if (maskv(i,j) .eq. 0 .and. m .ne. 0) line(i) = "!"           
            if (maskv(i,j) .ne. 0 .and. m .ne. 0) line(i) = "?"           
          enddo
          write (10,'(1000a)') (line(i),i=1,imt)
        enddo
        write (10,*) ' '
        
        print*, "Warning: check rivers.tmp, especially discharge points"
        print*, "         # = ocean point"
        print*, "         ! = discharge to an ocean point"
        print*, "         ? = discharge to a land point"
        print*, "rivers.tmp may be edited and moved to rivers.map but" 
        print*, "discharge points will have to be added by hand."
        
      endif

!=======================================================================
!     write netcdf data
!=======================================================================

      call opennew ("../rivers.nc", iou)
      call redef (iou)
      call defdim ('xt', iou, imt, id_xt)
      call defdim ('yt', iou, jmt, id_yt)
      call defvar ('xt', iou, 1, (/id_xt/), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
      call defvar ('yt', iou, 1, (/id_yt/), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
      call defvar ('rivers', iou, 2, (/id_xt,id_yt/), 0., 1000.
     &, ' ', 'I', 'basin river number', ' ' ,'1')
      call defvar ('discharge', iou, 2, (/id_xt,id_yt/), 0., 1000.
     &, ' ', 'I', 'discharge river number', ' ' ,'1')
      call defvar ('weights', iou, 2, (/id_xt,id_yt/), 0., 1000.
     &, ' ', 'F', 'river discharge weights', ' ' ,'1')
      call enddef (iou)
      call putvara ('xt', iou, imt, (/1/), (/imt/), xt, 1., 0.)
      call putvara ('yt', iou, jmt, (/1/), (/jmt/), yt, 1., 0.)
      call putvara ('rivers', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &,  var(:,:,1), 1., 0.)
      call putvara ('discharge', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &,  var(:,:,2), 1., 0.)
      call putvara ('weights', iou, imt*jmt, (/1,1/), (/imt,jmt/)
     &,  var(:,:,3), 1., 0.)
      call closefile (iou)

      end
