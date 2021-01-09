      subroutine mom_tbt_out (fname, ids, ide, jds, jde, imt, jmt, km
     &,                       nt, xt, yt, zt, xu, yu, zw, dxt, dyt, dzt
     &,                       dxu, dyu, dzw, ntrec, timunit, expnam
     &,                       avgper, relyr, stamp, mapt, tbt, stf
     &,                       tlat, tlon, ulat, ulon
     &,                       kmt, mskhr, tm)
!=======================================================================
!     output routine for ocean tracer balances

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, is, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!     imt, jmt ...       = global array dimensions
!     xt, yt ...         = global axes
!     ntrec              = number of time record (zero if not defined)
!     timunit            = time unit
!     expnam             = experiment name
!     avgper             = length of averaging period
!     relyr              = time in years
!     stamp              = time stamp
!     pass               = atmospheric transmission coefficient
!     tbt ...            = data to be written
!     tm                 = tracer grid mask

!   outputs:
!     ntrec              = number of time record in file

!   local variables
!     igs, ige, jgs, jge = global write domain start and end indicies
!     ig, jg             = global write domain size
!     ils, ile, jls, jle = local domain start and end indicies
!     it                 = t grid axis definitions (x,y,t default)
!     iu                 = u grid axis definitions (x,y,t default)
!     is                 = start for write on each axis (x,y,t default)
!     ic                 = count for write on each axis (x,y,t default)
!     id_...             = id's for axis (used for it, iu or defvar)
!     iou                = io unit (ncid)
!     ln                 = length of data to be written

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_tbt

      implicit none

      character(*) :: fname, timunit, expnam, stamp
      character(120) :: name, title, new_file_name
      character(3) :: a3
      character(10) :: mapt(nt)

      integer i, iou, j, ln, n, ntrec, imt, jmt, km, nt
      integer ids, ide, jds, jde, igs, ige, ig, jgs, jge, jg, kgs, kge
      integer kg, ils, ile, jls, jle, kls, kle, is(10), ic(10), it(10)
      integer iu(10), id_time, id_xt, id_xu, id_yt, id_yu, id_zt, id_zw
      integer id_xt_e, id_xu_e, id_yt_e, id_yu_e, id_zt_e,  id_zw_e
      integer kmt(ids:ide,jds:jde), mskhr(ids:ide,jds:jde)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      real xt(imt), xu(imt), yt(jmt), yu(jmt), zt(km), zw(km), avgper
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt), dzt(km), dzw(km)
      real tbt(ids:ide,jds:jde,km,nt,11), stf(ids:ide,jds:jde,nt)
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real ulat(ids:ide,jds:jde), ulon(ids:ide,jds:jde)
      real tm(ids:ide,jds:jde,km)
      real rkmt(ids:ide,jds:jde), rmskhr(ids:ide,jds:jde)
      real relyr, tmp, xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1)
      real yu_e(jmt+1), zt_e(km+1), zw_e(km+1), undef
      real c0, c1, c100, c400, c1e6, C2K, p1, p001, p035, cal2J
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpijk(:,:,:), tmpijkm(:,:,:)
      real, allocatable :: tmpi(:), tmpj(:), tmpk(:)
      real, allocatable :: tmpie(:), tmpje(:), tmpke(:)

      logical defined

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c400 = 400.
      c1e6 = 1.e6
      C2K = 273.15
      p1 = 0.1
      p001 = 0.001
      p035 = 0.035
      cal2J = 2.389e-05
      undef = 0.

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      defined = .true.
      name = new_file_name (fname)
      if (ntrec .eq. 0) defined = .false.
      call opennext (name, relyr, ntrec, iou)

!-----------------------------------------------------------------------
!     global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
# if defined cyclic
      igs = 2
      ige = imt-1
# else
      igs = 1
      ige = imt
# endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1
      kgs = 1
      kge = km
      kg  = kge-kgs+1

      if (.not. defined .or. ntrec .eq. 0) then

        if (ntrec .ne. 1) then
          call closefile (iou)
          call opennew (name, iou)
          ntrec = 1
        endif

!-----------------------------------------------------------------------
!       start definitions
!-----------------------------------------------------------------------
        call redef (iou)

!-----------------------------------------------------------------------
!       set global attributes
!-----------------------------------------------------------------------
        if (avgper .gt. 1.e-6) then
          title = 'Time Average Tracer Term Balances'
        else
          title = 'Snapshot Tracer Term Balances'
        endif
        call putglobal (iou, name, title, expnam, timunit)

!-----------------------------------------------------------------------
!       define dimensions
!-----------------------------------------------------------------------
        call defdim ('time', iou, 0, id_time)
        call defdim ('xt', iou, ig, id_xt)
        call defdim ('yt', iou, jg, id_yt)
        call defdim ('zt', iou, kg, id_zt)
        call defdim ('xu', iou, ig, id_xu)
        call defdim ('yu', iou, jg, id_yu)
        call defdim ('zw', iou, kg, id_zw)
        call defdim ('xt_edges', iou, ig+1, id_xt_e)
        call defdim ('yt_edges', iou, jg+1, id_yt_e)
        call defdim ('zt_edges', iou, kg+1, id_zt_e)
        call defdim ('xu_edges', iou, ig+1, id_xu_e)
        call defdim ('yu_edges', iou, jg+1, id_yu_e)
        call defdim ('zw_edges', iou, kg+1, id_zw_e)

!-----------------------------------------------------------------------
!       define 1d data (t)
!-----------------------------------------------------------------------
        it(1) = id_time
        call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &,   'time since initial condition', 'time', trim(timunit))
        call defvar ('year', iou, 1, it, c0, c0, ' ', 'F'
     &,   'year', ' ',' ')
        call defvar ('month', iou, 1, it, c0, c0, ' ', 'F'
     &,   'month', ' ',' ')
        call defvar ('day', iou, 1, it, c0, c0, ' ', 'F'
     &,   'day', ' ',' ')
        call defvar ('hour', iou, 1, it, c0, c0, ' ', 'F'
     &,   'hour', ' ',' ')
        call defvar ('minute', iou, 1, it, c0, c0, ' ', 'F'
     &,   'minute', ' ',' ')
        call defvar ('second', iou, 1, it, c0, c0, ' ', 'F'
     &,   'second', ' ',' ')
        if (avgper .gt. 1.e-6) then
          call defvar ('period', iou, 1, it, c0, c0, ' ', 'F'
     &,     'averaging period', ' ',' ')
        endif

!-----------------------------------------------------------------------
!       define 1d data (x, y or z)
!-----------------------------------------------------------------------
        it(1) = id_xt
        call defvar ('xt', iou, 1, it, c0, c0, 'X', 'D'
     &,   'longitude of the t grid', 'grid_longitude', 'degrees_east')
        it(1) = id_yt
        call defvar ('yt', iou, 1, it, c0, c0, 'Y', 'D'
     &,   'latitude of the t grid', 'grid_latitude', 'degrees_north')
        it(1) = id_zt
        call defvar ('zt', iou, 1, it, c0, c0, 'Z', 'D'
     &,   'depth of the t grid', 'depth', 'm')
        it(1) = id_xu
        call defvar ('xu', iou, 1, it, c0, c0, 'X', 'D'
     &,   'longitude of the u grid', 'grid_longitude', 'degrees_east')
        it(1) = id_yu
        call defvar ('yu', iou, 1, it, c0, c0, 'Y', 'D'
     &,   'latitude of the u grid', 'grid_latitude', 'degrees_north')
        it(1) = id_zw
        call defvar ('zw', iou, 1, it, c0, c0, 'Z', 'D'
     &,   'depth of the w grid', 'depth', 'm')
        it(1) = id_xt_e
        call defvar ('xt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'longitude of t grid edges', ' ', 'degrees')
        it(1) = id_yt_e
        call defvar ('yt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'latitude of t grid edges', ' ', 'degrees')
        it(1) = id_zt_e
        call defvar ('zt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'depth of t grid edges', ' ', 'm')
        it(1) = id_xu_e
        call defvar ('xu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'longitude of u grid edges', ' ', 'degrees')
        it(1) = id_yu_e
        call defvar ('yu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'latitude of u grid edges', ' ', 'degrees')
        it(1) = id_zw_e
        call defvar ('zw_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'depth of w grid edges', ' ', 'm')

!-----------------------------------------------------------------------
!       define 2d data (x,y)
!-----------------------------------------------------------------------
        it(1) = id_xt
        iu(1) = id_xu
        it(2) = id_yt
        iu(2) = id_yu
        call defvar ('kmt', iou, 2, it, c0, c1e6, ' ', 'I'
     &,   'ocean grid depth level', 'model_level_number' ,'1')
        call defvar ('mskhr', iou, 2, it, c0, c1e6, ' ', 'I'
     &,   'horizontal region mask', ' ' ,'1')
        call defvar ('tlat', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &,   'tracer grid latitude', 'latitude', 'degrees')
        call defvar ('tlon', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &,   'tracer grid longitude', 'longitude', 'degrees')
        call defvar ('ulat', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &,   'velocity grid latitude', 'latitude', 'degrees')
        call defvar ('ulon', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &,   'velocity grid longitude', 'longitude', 'degrees')

!-----------------------------------------------------------------------
!       define 3d data (x,y,t)
!-----------------------------------------------------------------------
        it(1) = id_xt
        iu(1) = id_xu
        it(2) = id_yt
        iu(2) = id_yu
        it(3) = id_time
        iu(3) = id_time
        do n=1,nt
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
          call defvar ('flux_'//trim(a3), iou ,3, it, -c1e6, c1e6, ' '
     &,     'F', 'tracer flux '//trim(a3), 'tracer_flux_'//trim(a3)
     &,     ' ')
        enddo

!-----------------------------------------------------------------------
!       define time dependent 4d data (x,y,z,t)
!-----------------------------------------------------------------------
        it(1) = id_xt
        iu(1) = id_xu
        it(2) = id_yt
        iu(2) = id_yu
        it(3) = id_zt
        iu(3) = id_zt
        it(4) = id_time
        iu(4) = id_time
        do n=1,nt
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
          call defvar ('tracer_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,     ' ', 'F', 'tracer '//trim(a3)
     &,     'tracer_'//trim(a3), ' ')
          call defvar ('d_'//trim(a3)//'dt', iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', 'd_'//trim(a3)//'/dt',' ', ' ')
          call defvar ('-xadv_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F','-zonal advection', ' ', ' ')
          call defvar ('-yadv_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', '-meridional advection', ' ', ' ')
          call defvar ('-zadv_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', '-vertical advection', ' ', ' ')
          call defvar ('xdiff_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', 'zonal diffusion', ' ', ' ')
          call defvar ('ydiff_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', 'meridional diffusion', ' ', ' ')
          call defvar ('zdiff_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', 'vertical diffusion', ' ', ' ')
# if defined source_term
          call defvar ('source_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', 'change due to sources', ' ', ' ')
# endif
          call defvar ('convect_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', 'change due to convection', ' ', ' ')
# if defined fourfil || defined firfil
          call defvar ('filter_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,     ' ', 'F', 'change due to filtering', ' ', ' ')
# endif
        enddo

!-----------------------------------------------------------------------
!       end definitions
!-----------------------------------------------------------------------
        call enddef (iou)
        if (ntrec .eq. 0) ntrec = 1

      endif

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)
      kls = max(1,kgs)
      kle = min(km,kge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpijm(ils:ile,jls:jle) )
      allocate ( tmpijk(ils:ile,jls:jle,kls:kle) )
      allocate ( tmpijkm(ils:ile,jls:jle,kls:kle) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      call putvars ('time', iou, ntrec, relyr, c1, c0)
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = nyear
      call putvars ('year', iou, ntrec, tmp, c1, c0)
      tmp = nmonth
      call putvars ('month', iou, ntrec, tmp, c1, c0)
      tmp = nday
      call putvars ('day', iou, ntrec, tmp, c1, c0)
      tmp = nhour
      call putvars ('hour', iou, ntrec, tmp, c1, c0)
      tmp = nmin
      call putvars ('minute', iou, ntrec, tmp, c1, c0)
      tmp = nsec
      call putvars ('second', iou, ntrec, tmp, c1, c0)
      if (avgper .gt. 1.e-6) then
        call putvars ('period', iou, ntrec, avgper, c1, c0)
      endif

      if (ntrec .eq. 1 .or. .not. defined) then

!-----------------------------------------------------------------------
!       write 1d data (x, y or z)
!-----------------------------------------------------------------------
        allocate ( tmpi(igs:ige) )
        allocate ( tmpj(jgs:jge) )
        allocate ( tmpk(kgs:kge) )
        allocate ( tmpie(igs:ige+1) )
        allocate ( tmpje(jgs:jge+1) )
        allocate ( tmpke(kgs:kge+1) )

        ib(1) = 1
        ic(1) = ig
        tmpi(igs:ige) = xt(igs:ige)
        call putvara ('xt', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = xu(igs:ige)
        call putvara ('xu', iou, ig, ib, ic, tmpi, c1, c0)

        ic(1) = jg
        tmpj(jgs:jge) = yt(jgs:jge)
        call putvara ('yt', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = yu(jgs:jge)
        call putvara ('yu', iou, jg, ib, ic, tmpj, c1, c0)

        ic(1) = kg
        tmpk(kgs:kge) = zt(kgs:kge)
        call putvara ('zt', iou, kg, ib, ic, tmpk, c100, c0)
        tmpk(kgs:kge) = zw(kgs:kge)
        call putvara ('zw', iou, kg, ib, ic, tmpk, c100, c0)

        ic(1) = ig + 1
        call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xt_e(igs:ige+1)
        call putvara ('xt_edges', iou, ig+1, ib, ic, tmpie, c1, c0)
        call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xu_e(igs:ige+1)
        call putvara ('xu_edges', iou, ig+1, ib, ic, tmpie, c1, c0)

        ic(1) = jg + 1
        call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
        call putvara ('yt_edges', iou, jg+1, ib, ic, tmpje, c1, c0)
        call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
        call putvara ('yu_edges', iou, jg+1, ib, ic, tmpje, c1, c0)

        ic(1) = kg + 1
        call edge_maker (1, zt_e, zt, dzt, zw, dzw, km)
        tmpke(kgs:kge+1) = zt_e(kgs:kge+1)
        call putvara ('zt_edges', iou, kg+1, ib, ic, tmpke, c100, c0)
        call edge_maker (2, zw_e, zt, dzt, zw, dzw, km)
        tmpke(kgs:kge+1) = zw_e(kgs:kge+1)
        call putvara ('zw_edges', iou, kg+1, ib, ic, tmpke, c100, c0)

        deallocate ( tmpi )
        deallocate ( tmpj )
        deallocate ( tmpk )
        deallocate ( tmpie )
        deallocate ( tmpje )
        deallocate ( tmpke )

!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
        tmpij(ils:ile,jls:jle) = kmt(ils:ile,jls:jle)
        call putvara ('kmt', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = mskhr(ils:ile,jls:jle)
        call putvara ('mskhr', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlat(ils:ile,jls:jle)
        call putvara ('tlat', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlon(ils:ile,jls:jle)
        call putvara ('tlon', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulat(ils:ile,jls:jle)
        call putvara ('ulat', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulon(ils:ile,jls:jle)
        call putvara ('ulon', iou, ln, ib, ic, tmpij, c1, c0)

      endif

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      do n=1,nt
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpij(ils:ile,jls:jle) = stf(ils:ile,jls:jle,n)
        call putvaramsk('flux_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   tmpijm, c1, c0)
      enddo

!-----------------------------------------------------------------------
!     write 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      is(1) = ils-igs+1
      ic(1) = ile-ils+1
      is(2) = jls-jgs+1
      ic(2) = jle-jls+1
      is(3) = kls-kgs+1
      ic(3) = kle-kls+1
      is(4) = ntrec
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      do n=1,nt
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,1)
        call putvaramsk('d_'//trim(a3)//'dt', iou, ln, is, ic, tmpijk
     &,   tmpijkm)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,2)
        call putvaramsk('-xadv_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,3)
        call putvaramsk('-yadv_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,4)
        call putvaramsk('-zadv_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,5)
        call putvaramsk('xdiff_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,6)
        call putvaramsk('ydiff_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,7)
        call putvaramsk('zdiff_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
# if defined source_term
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,8)
        call putvaramsk('source_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
# endif
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,9)
        call putvaramsk('convect_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
# if defined fourfil || defined firfil
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,10)
        call putvaramsk('filter_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
# endif
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,11)
        call putvaramsk('tracer_'//trim(a3), iou, ln, is, ic, tmpijk
     &,   tmpijkm)
      enddo

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)

      deallocate (tmpij)
      deallocate (tmpijm)
      deallocate (tmpijk)
      deallocate (tmpijkm)

#endif

      return
      end
