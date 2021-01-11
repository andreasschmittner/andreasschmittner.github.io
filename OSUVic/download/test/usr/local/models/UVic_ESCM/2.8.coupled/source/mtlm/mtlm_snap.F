      subroutine mtlm_snap_def (fname, imt, jmt, NPFT, NTYPE, xt, yt
     &,                         timunit, expnam, runstamp)
!=======================================================================
!     definition routine for land time averages or snapshots

!   inputs:
!     fname         = file name
!     imt, jmt ...  = global array dimensions
!     xt, yt ...    = global axes
!     timunit       = time unit
!     expnam        = experiment name
!     runstamp      = run stamp

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_mtlm

      implicit none

      character(*) :: fname, timunit, expnam, runstamp

      integer iou, j, imt, jmt
      integer igs, ige, ig, jgs, jge, jg
      integer it(10)
      integer iu(10), id_time, id_xt, id_xu, id_yt, id_yu, id_pft
      integer id_type, id_xt_e, id_xu_e, id_yt_e, id_yu_e, id_pft_e
      integer id_type_e

      integer NPFT, NTYPE

      real xt(imt), yt(jmt)

      real c0, c1, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global attributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)
      call putatttext (iou, 'global', 'time_unit', timunit)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      call defdim ('xt', iou, ig, id_xt)
      call defdim ('yt', iou, jg, id_yt)
      call defdim ('xu', iou, ig, id_xu)
      call defdim ('yu', iou, jg, id_yu)
      call defdim ('pft', iou, NPFT, id_pft)
      call defdim ('type', iou, NTYPE, id_type)
      call defdim ('xt_edges', iou, ig+1, id_xt_e)
      call defdim ('yt_edges', iou, jg+1, id_yt_e)
      call defdim ('xu_edges', iou, ig+1, id_xu_e)
      call defdim ('yu_edges', iou, jg+1, id_yu_e)
      call defdim ('pft_edges', iou, NPFT+1, id_pft_e)
      call defdim ('type_edges', iou, NTYPE+1, id_type_e)

!-----------------------------------------------------------------------
!     define 1d data (t)
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
      call defvar ('period', iou, 1, it, c0, c0, ' ', 'F'
     &,   'averaging period', ' ','days')

!-----------------------------------------------------------------------
!     define 1d data (x, y or n)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('xt', iou, 1, it, c0, c0, 'X', 'D'
     &,   'longitude of the t grid', 'grid_longitude', 'degrees_east')
      it(1) = id_yt
      call defvar ('yt', iou, 1, it, c0, c0, 'Y', 'D'
     &,   'latitude of the t grid', 'grid_latitude', 'degrees_north')
      it(1) = id_xu
      call defvar ('xu', iou, 1, it, c0, c0, 'X', 'D'
     &,   'longitude of the u grid', 'grid_longitude', 'degrees_east')
      it(1) = id_yu
      call defvar ('yu', iou, 1, it, c0, c0, 'Y', 'D'
     &,   'latitude of the u grid', 'grid_latitude', 'degrees_north')
      it(1) = id_pft
      call defvar ('pft', iou, 1, it, c0, c0, ' ', 'F'
     &,   'plant functional type', ' ', ' ')
      it(1) = id_type
      call defvar ('type', iou, 1, it, c0, c0, ' ', 'F'
     &,   'land type', ' ', ' ')
      it(1) = id_xt_e
      call defvar ('xt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'longitude of t grid edges', ' ', 'degrees')
      it(1) = id_yt_e
      call defvar ('yt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'latitude of t grid edges', ' ', 'degrees')
      it(1) = id_xu_e
      call defvar ('xu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'longitude of u grid edges', ' ', 'degrees')
      it(1) = id_yu_e
      call defvar ('yu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &,   'latitude of u grid edges', ' ', 'degrees')
      it(1) = id_pft_e
      call defvar ('pft_edges', iou, 1, it, c0, c0, ' ', 'F'
     &,   'plant functional type', ' ', ' ')
      it(1) = id_type_e
      call defvar ('type_edges', iou, 1, it, c0, c0, ' '
     &,   'F', 'land type', ' ', ' ')

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      call defvar ('tlat', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &,   'tracer grid latitude', 'latitude', 'degrees')
       call defvar ('tlon', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &,   'tracer grid longitude', 'longitude', 'degrees')
       call defvar ('ulat', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &,   'velocity grid latitude', 'latitude', 'degrees')
      call defvar ('ulon', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &,   'velocity grid longitude', 'longitude', 'degrees')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_time
      call defvar ('TS1', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'soil temperature', 'soil_temperature', 'K')
# if !defined uvic_embm
      call defvar ('TSTAR_GB', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'surface temperature', 'surface_temperature', 'K')
      call defvar ('ALBLAND', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'surface albedo', 'surface_albedo', '1')
      call defvar ('ET', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'evapotranspiration', 'evapotranspiration', 'kg m-2 s-1')
      call defvar ('M', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'soil moisture', 'soil_moisture_content', 'kg m-2')
# endif
      call defvar ('CS', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'soil carbon', 'soil_carbon_content', 'kg m-2 ')
      call defvar ('RESP_S', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'soil respiration', 'litter_carbon_flux'
     &,   'kg m-2 s-1')

!-----------------------------------------------------------------------
!     define 4d data (x,y,pft,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_pft
      it(4) = id_time
      call defvar ('GPP', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'gross primary productivity'
     &,   'gross_primary_productivity_of_carbon', 'kg m-2 s-1')
      call defvar ('NPP', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'net primary productivity'
     &,   'net_primary_productivity_of_carbon', 'kg m-2 s-1')
      call defvar ('HT', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'canopy height', 'canopy_height', 'm')
      call defvar ('LAI', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'leaf area index ', 'leaf_area_index', '1')
      call defvar ('C_VEG', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'vegetation carbon', 'vegetation_carbon_content', 'kg m-2 ')

!-----------------------------------------------------------------------
!     define 4d data (x,y,type,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_type
      it(4) = id_time
       call defvar ('FRAC', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &, 'areal coverage', 'vegetation_area_fraction', '1')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)

      return
      end

      subroutine mtlm_snap_out (fname, ids, ide, jds, jde, imt, jmt
     &,                         POINTS, NPFT, NTYPE, xt, yt, xu, yu
     &,                         dxt, dyt, dxu, dyu, avgper, relyr
     &,                         stamp, land_map, TS1,  CS, RESP_S, FRAC
     &,                         GPP, NPP, HT, LAI, C_VEG
#if !defined uvic_embm
     &,                         TSTAR_GB, ALBLAND, ET, M
#endif
     &,                         tlat, tlon, ulat, ulon, ntrec)
!=======================================================================
!     output routine for land time averages or snapshots

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname        = file name
!     ids, ide ... = start and end index for data domain
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     dxt, dyt ... = grid widths
!     avgper       = length of averaging period
!     relyr        = time in years
!     stamp        = time stamp
!     land_map     = land map
!     TS1, ...     = data to be written

!   outputs:
!     ntrec        = number of time record in file

!     based on code by: M. Eby
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, j, ln, n, ntrec, imt, jmt, ids, ide, jds, jde, igs
      integer ige, ig, jgs, jge, jg, ils, ile, jls, jle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec
      integer POINTS, NPFT, NTYPE, land_map(imt,jmt)

      real xt(imt), xu(imt), yt(jmt), yu(jmt), dxt(imt), dxu(imt)
      real dyt(jmt), dyu(jmt), avgper

      real relyr, tmp, xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1)
      real yu_e(jmt+1), tmpmask(imt,jmt), data(imt,jmt)
      real pft(NPFT), type(NTYPE), pft_e(NPFT+1), type_e(NTYPE+1)
      real c0, c1

# if !defined uvic_embm
      real TSTAR_GB(POINTS), ALBLAND(POINTS), ET(POINTS), M(POINTS)
# endif
      real TS1(POINTS), CS(POINTS), RESP_S(POINTS), FRAC(POINTS,NTYPE)
      real GPP(POINTS,NPFT), NPP(POINTS,NPFT), HT(POINTS,NPFT)
      real LAI(POINTS,NPFT), C_VEG(POINTS,NPFT)
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real ulat(ids:ide,jds:jde), ulon(ids:ide,jds:jde)
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      c0 = 0.
      c1 = 1.

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, relyr, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpijm(ils:ile,jls:jle) )

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

      if (ntrec .eq. 1) then

!-----------------------------------------------------------------------
!       write 1d data (x, y or z)
!-----------------------------------------------------------------------
        allocate ( tmpi(igs:ige) )
        allocate ( tmpj(jgs:jge) )
        allocate ( tmpie(igs:ige+1) )
        allocate ( tmpje(jgs:jge+1) )

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

        deallocate ( tmpi )
        deallocate ( tmpj )
        deallocate ( tmpie )
        deallocate ( tmpje )

        do n=1, NPFT
          pft(n) = n
          pft_e(n) = pft(n) - 0.5
        enddo
        pft_e(NPFT+1) = pft(NPFT) + 0.5
        ic(1) = NPFT
        call putvara ('pft', iou, NPFT, ib, ic, pft, c1, c0)
        ic(1) = NPFT + 1
        call putvara ('pft_edges', iou, NPFT+1, ib, ic, pft_e, c1, c0)
        do n=1, NTYPE
          type(n) = float(n)
          type_e(n) = type(n) - 0.5
        enddo
        type_e(NTYPE+1) = type(NTYPE) + 0.5
        ic(1) = NTYPE
        call putvara ('type', iou, NTYPE, ib, ic, type, c1, c0)
        ic(1) = NTYPE + 1
        call putvara ('type_edges', iou, NTYPE+1, ib, ic, type_e
     &,   c1, c0)

!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
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
      tmpmask(:,:) = 0.
      where (land_map(:,:) .gt. 0) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      call unloadland (POINTS, TS1, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('TS1', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
# if !defined uvic_embm
      call unloadland (POINTS, TSTAR_GB, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('TSTAR_GB', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      call unloadland (POINTS, ALBLAND, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('ALBLAND', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      call unloadland (POINTS, ET, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('ET', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
      call unloadland (POINTS, M, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('M', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
# endif
      call unloadland (POINTS, CS, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('CS', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
      call unloadland (POINTS, RESP_S, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('RESP_S', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)

!-----------------------------------------------------------------------
!     write 4d data (x,y,pft,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ic(3) = 1
      ib(4) = ntrec
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      do n=1,npft
        ib(3) = n
        call unloadland (POINTS, GPP(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('GPP', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
        call unloadland (POINTS, NPP(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('NPP', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
        call unloadland (POINTS, HT(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('HT', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
        call unloadland (POINTS, LAI(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('LAI', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
        call unloadland (POINTS, C_VEG(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('C_VEG', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
      enddo

!-----------------------------------------------------------------------
!     write 4d data (x,y,type,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ic(3) = 1
      ib(4) = ntrec
      ic(4) = 1
      do n=1,ntype
        ib(3) = n
        call unloadland (POINTS, FRAC(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('FRAC', iou, ln, ib, ic, tmpij, tmpijm, c1, c0)
      enddo

!-----------------------------------------------------------------------
!     close the file
!-----------------------------------------------------------------------
      call closefile (iou)
#endif

      return
      end
