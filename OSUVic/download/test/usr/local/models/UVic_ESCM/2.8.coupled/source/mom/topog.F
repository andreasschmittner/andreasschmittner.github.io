      subroutine topog (kmt, map, xt, yt, zt, xu, yu, zw, imt2,jmt2,km2)

!-----------------------------------------------------------------------
!     construct the topographic mask "kmt" which determines the geometry
!     and topography by defining the number of vertical levels on model
!     "t" points.

!     A base kmt field is constructed using one of the following:

!      (A) Construct a kmt field using one of the following:
!        -Dscripps_kmt : interpolate Scripps 1 degree topography to
!                        "kmt" on MOM grid
!        -Dread_my_kmt : read in your own "kmt" from a file "kmt.dta"
!        -Didealized_kmt : generate an "idealized" not very realistic
!                          topography without reading the database
!        -Drectangular_box : sets "kmt" to a flat bottom rectangular box

!       Options:
!        -Dcyclic        : sets cyclic conditons in zonal direction
!                          (otherwise solid_walls is assumed)
!        -Dflat_bottom   : sets "kmt" over ocean to maximum levels (km)

!     Problems may exist in the base kmt field. They are fixed by
!     iterating over the following sequence of operations.

!      (B) Fill isolated cells (optional)

!       Option:
!        -Dfill_isolated_cells : These are places where "T" cells are
!                                 deeper than nearest neighbors. Also
!                                 fills in "T" cells with land at all
!                                 four corners.

!      (C) Add USER CHANGES to kmt field (optional)

!      (D) Limit minimum levels based on parameter "kmt_min"
!          "kmt_min" is the minimum number of vertical levels
!           permitted at any ocean "T" cell.  kmt_min must be at >= 2.

!           Options:
!           -Dfill_shallow: makes all "t" cells for which
!                           kmt(i,jrow) < kmt_min into land.
!           -Ddeepen_shallow: deepens such cell to depth kmt_min.
!           -Dround_shallow: changes depth to the nearer of 0 and
!                            zw(kmt_min)

!      (E) Fix perimeter violations (cells in the ocean perimeter of
!                                    two different land masses)
!          Options:
!           -Dfill_perimeter_violations : build a land bridge between
!                                         these land masses (default)
!           -Dwiden_perimeter_violations: widen straits between two
!                                         different land masses to a
!                                         minimum of two "t" cells

!-----------------------------------------------------------------------
!     input:

!     xt  = longitudes of model grid "t" points (degrees)
!     yt  = latitudes of model grid "t" points (degrees)
!     zt  = depths of box "u/v" points on model grid (cm)
!     xu  = longitudes of model grid "u" points (degrees)
!     yu  = latitudes of model grid "u" points (degrees)
!     zw  = depths of box bottoms on model grid "t" boxes (cm)
!     imt = number of model longitudes
!     jmt = number of model latitudes
!     km  = maximum number of levels

!     output:

!     kmt = number of vertical levels on "t" points
!           0 is a land point.
!           n implies n boxes between the surface and bottom
!             where (n >= kmt_min) & (n <= km)
!     map = array indicating land masses and their perimeters
!           0 is a mid ocean point
!           n > 0 is a land mass
!           n < 0 is the ocean perimeter of land mass |n|

!-----------------------------------------------------------------------

#include "stdunits.h"
#include "size.h"
#include "isleperim.h"
#if defined uvic_read_kmt_map
      character(1) :: line(imt+10)
#else
      real tmpij(imt,jmt)
#endif
      integer ib(10), ic(10)

!     kmt_min is the minimum permitted number of levels at an ocean
!     "T" cell.  kmt_min must be at least 2.

      parameter (kmt_min = 2)

      dimension kmt(imt,jmt), kmu(imt,jmt), xt(imt), yt(jmt), zt(km)
      dimension xu(imt), yu(jmt), zw(km)

      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)

      character(40) :: cktext
      character(120) :: iotext
      character(60) :: expnam
      character(32) :: stamp

      character(120) :: fname, new_file_name
      logical exists

!     the variable auto_kmt_changes tells whether any changes have been
!     made to the kmt field as a result of define options.
!     nchanges counts how many changes in each group.

      auto_kmt_changes = .false.
      nchanges = 0

!-----------------------------------------------------------------------
!     check that grid sizes in argument list and file "size.h" agree
!-----------------------------------------------------------------------

      call size_check (imt2, jmt2, km2, 'topog', 'stop')

      ncase = 0

      write (stdout,'(//,t37,a,/)')
     & 'T O P O G R A P H Y    G E N E R A T I O N'
      write (stdout,'(a,i4,/)') 'kmt_min = ', kmt_min

#if defined scripps_kmt

!-----------------------------------------------------------------------
!     construct model "kmt" based on Scripps 1 deg topography data.
!-----------------------------------------------------------------------

      ncase = ncase + 1
      write (stdout,'(/a/)')
     & ' Interpolating Scripps 1 degree topography to the MOM grid:'

      call scripp (kmt, xt, yt, zt, xu, yu, zw, imt, jmt, km)
#endif
#if defined etopo_kmt

!-----------------------------------------------------------------------
!     construct model "kmt" based on NGDC ETOPO5 1/12 deg topography
!-----------------------------------------------------------------------

      ncase = ncase + 1
      write (stdout,'(/a/)')
     & ' Interpolating ETOPO 1/12 degree topography to the MOM grid:'

      call etopo (kmt, xt, yt, zt, xu, yu, zw, imt, jmt, km)

!-----------------------------------------------------------------------
!     show the kmt field
!-----------------------------------------------------------------------

      write (stdout,'(/,20x,a/)') ' The original "kmt" field follows:'
      call iplot (kmt, imt, imt, jmt)

#endif
#if defined read_my_kmt || uvic_read_kmt_map

      ncase = ncase + 1
# if defined uvic_read_kmt_map

      print*, jmt, imt
      fname = new_file_name ("kmt.map")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'topog in topog.f'
      endif
      call getunit (iou, 'kmt.map'
     &,             'formatted sequential rewind ieee')
      do j=jmt,1,-1
        read (iou,'(1000a)') (line(i),i=1,imt)
        do i=1,imt
          kmt(i,j) = ichar(line(i))
          if (kmt(i,j) .eq. 46) n = 46
          if (kmt(i,j) .ge. 65 .and. kmt(i,j) .le. 90) n = 64
          if (kmt(i,j) .ge.97) n = 70
          kmt(i,j) = kmt(i,j) - n
        enddo
      enddo
      call relunit (iou)
# else
      c0 = 0.
      c1 = 1.
      fname = new_file_name ("kmt.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'topog in topog.f'
      endif
      ib(:) = 1
      ic(:) = imt
      ic(2) = jmt
      call openfile (fname, iou)
      call getvara ('kmt', iou, imt*jmt, ib, ic, tmpij, c1, c0)
      call closefile (iou)
      kmt(:,:) = tmpij(:,:)
# endif

      nerror = 0
      do i=1,imt
        do j=1,jmt
          if (kmt(i,j) .lt. 0 .or. kmt(i,j) .gt. km) then
            nerror = nerror + 1
          endif
        enddo
      enddo
      if (nerror .gt. 0) then
        write (stdout,'(a,i4,a,a,i4,/)')
     &    '==>ERROR: There are ', nerror, ' depth values in the file '
     &,   '"kmt.dta" that are outside the range 0 <= kmt(i,j) <= km = '
     &,   km
        stop '==>topog'
      endif
      write (stdout,'(/a/)') '==> The kmt field has been imported.'
#endif
#if defined idealized_kmt

!-----------------------------------------------------------------------
!     construct a highly "idealized" world ... piece by piece

!     note: the purpose of this geometry/topography is to automatically
!           map into arbitrary resolution as grid dimensions "imt" and
!           "jmt" are changed, thereby facilitating the implementation
!           and verification of the model on various computer platforms
!           without referencing the topographic data base.  Although it
!           somewhat resembles the real world, it is NOT realistic.
!-----------------------------------------------------------------------

      ncase = ncase + 1
      write (stdout,'(/a/a/)')
     &' Constructing an idealized test case "kmt" which only crudely'
     &,' approximates world geometry. Topography is unrealistic.'

      call idealized (kmt, xt, yt, zw, imt, jmt, km)
#endif
#if defined rectangular_box
      ncase = ncase + 1
      write (stdout,'(/a,i4,a/)')
     & 'Constructing a flat bottomed rectangular box with ',km,' levels'

      do jrow=2,jmt-1
        do i=2,imt-1
          kmt(i,jrow) = km
        enddo
      enddo
#endif

!-----------------------------------------------------------------------
!     check that only one method of generating topography was chosen
!-----------------------------------------------------------------------

      if (ncase .eq. 0) then
        write (stdout,'(/a/a/a/)')
     &  '=>Error: One of the following options must be enabled:'
     &, '      idealized_kmt, scripps_kmt, read_my_kmt, rectangular_box'
     &, '      etopo_kmt'
        stop '=>topog'
      elseif (ncase .gt. 1) then
        write (stdout,'(/a/a/)')
     &  '=>Error: Only one of the following options may be enabled:'
     &, '      idealized_kmt, scripps_kmt, read_my_kmt, rectangular_box'
        stop '=>topog'
      endif

#if defined flat_bottom
      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .ne. 0) then
            kmt(i,jrow) = km
          endif
        enddo
      enddo
#endif

!-----------------------------------------------------------------------
!     apply boundary conditions
!-----------------------------------------------------------------------

      call kmtbc (kmt, imt, jmt)
      call area_volume (kmt, xu, yu, zw)
      cksum1 = checksumi (kmt, imt, jmt, 'Original "kmt" checksum =')

!-----------------------------------------------------------------------
!     iterate on changes to "kmt" until all are consistant
!-----------------------------------------------------------------------

      npass = 0
      n_del_kmt = 0
1000  continue
      npass = npass + 1

#if defined fill_isolated_cells

!-----------------------------------------------------------------------
!     fill isolated cells or potholes at all levels in kmt
!     filled kmt array is just the maximum of the surrounding kmu levels
!     symmetry conditions are automatic, if set on kmu
!-----------------------------------------------------------------------

      do i=1,imt-1
        do j=1,jmt-1
          kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
        enddo
      enddo
      n = 0
      write (*,'(/a)') 'Searching for isolated T cells...'
      do i=2,imt-1
        do j=2,jmt-1
          k = max (kmu(i,j), kmu(i-1,j), kmu(i,j-1), kmu(i-1,j-1))
          if (kmt(i,j) .ne. k) n = n + 1
          kmt(i,j) = k
        enddo
      enddo
      if (n .gt. 0) then
        write (*,*) '-> Found ',n,' and filled them in.'
      else
        write (*,*) '-> None Found.'
      endif

      call kmtbc (kmt, imt, jmt)
#endif

!-----------------------------------------------------------------------
!     USER INPUT
!     other changes to "kmt" may be made here.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     ALL USER MODIFICATIONS MUST PRECEDE THESE FINAL CHECKS
!-----------------------------------------------------------------------

      call kmtbc (kmt, imt, jmt)
#if defined obc
# if defined obc_south
      do jrow=1,jmt
        do i=1,imt
          if (jrow .eq. 2) then
            write(6,*)' Changing Topography at southern OB. '
            kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i,jrow+1)),
     &                  min(kmt(i,jrow),kmt(i,jrow+2)))
            kmt(i,jrow-1)=kmt(i,jrow)
            kmt(i,jrow+1)=kmt(i,jrow)
            kmt(i,jrow+2)=kmt(i,jrow+1)
          endif
        enddo
      enddo
# endif
# if defined obc_north
      do jrow=1,jmt
        do i=1,imt
          if (jrow .eq. jmt-1) then
            write(6,*)' Changing Topography at northern OB. '
            kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i,jrow-1)),
     &                  min(kmt(i,jrow),kmt(i,jrow-2)))
            kmt(i,jrow+1)=kmt(i,jrow)
            kmt(i,jrow-1)=kmt(i,jrow)
            kmt(i,jrow-2)=kmt(i,jrow-1)
          endif
        enddo
      enddo
# endif
# if defined obc_east
      do jrow=1,jmt
        do i=1,imt
          if (i .eq. 2) then
            write(6,*)' Changing Topography at eastern OB. '
            kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i+1,jrow)),
     &                  min(kmt(i,jrow),kmt(i+2,jrow)))
            kmt(i-1,jrow)=kmt(i,jrow)
            kmt(i+1,jrow)=kmt(i,jrow)
            kmt(i+2,jrow)=kmt(i+1,jrow)
          endif
        enddo
      enddo
# endif
# if defined obc_west
      do jrow=1,jmt
        do i=1,imt
          if (i .eq. imt-1) then
            write(6,*)' Changing Topography at western OB. '
            kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i-1,jrow)),
     &                  min(kmt(i,jrow),kmt(i-2,jrow)))
            kmt(i+1,jrow)=kmt(i,jrow)
            kmt(i-1,jrow)=kmt(i,jrow)
            kmt(i-2,jrow)=kmt(i-1,jrow)
          endif
        enddo
      enddo
# endif
#endif

!-----------------------------------------------------------------------
!     limit the minimum number of levels. kmt_min should be >= 2
!-----------------------------------------------------------------------

      call min_depth (kmt, kmt_min, xu, yu, zw)
      call kmtbc (kmt, imt, jmt)

!-----------------------------------------------------------------------
!     test for "island perimeter violations"
!     channels between different land masses must be at least 2 ocean
!     cells wide to prevent conflicting values of the stream function
!     on these island perimeter points.
!-----------------------------------------------------------------------

      write (stdout,'(/a/)')
     & 'Searching for and correcting island PERIMETER VIOLATIONS...'
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
      call kmtbc (kmt, imt, jmt)
      call area_volume (kmt, xu, yu, zw)

      write (cktext,'(a,i2,a)') 'Pass # ',npass
     &,                         ' complete.  kmt checksum ='
      cksum2 = checksumi (kmt, imt, jmt, cktext)
      if (cksum1 .ne. cksum2) then
        cksum1 = cksum2
        if (npass .gt. 15) then
          print *,'=>Error: changes to "kmt" may be incompatable'
          print *,'         arbitrarily stopping after 15 passes'
          stop
        endif
        go to 1000
      endif
      print '(/a/)', 'All changes have now been made to kmt'

!-----------------------------------------------------------------------
!     end of iteration on "kmt" changes
!-----------------------------------------------------------------------

#  if !defined skip_island_map
      print '(/a/)', 'Use -Dskip_island_map eliminate the following map'
      call showmap (map, imt, jmt, linewidth)
#  else
      print '(/a/)', 'Remove -Dskip_island_map to show perimeter map'
#  endif

!-----------------------------------------------------------------------
!     show the kmt field
!-----------------------------------------------------------------------

#  if !defined skip_kmt_map
      print '(/a/)', 'Use -Dskip_kmt_map to eliminate the following map'
      write (stdout,'(/,20x,a/)') ' The "kmt" field follows:'
      call iplot (kmt, imt, imt, jmt)
#  else
      print '(/a/)', 'Remove -Dskip_kmt_map to show the kmt map'
#  endif
      call area_volume (kmt, xu, yu, zw)

!---------------------------------------------------------------------
!     compute a topography checksum
!---------------------------------------------------------------------

      cksum1 = checksumi (kmt, imt, jmt, 'Final "kmt" checksum =')

      return
      end

      function checksumi (kmt, im, jm, text)
      dimension kmt(im,jm)
      character(*) :: text
      cksum = 0.0
      do jrow=1,jm
        do i=1,im
          cksum = cksum + i*jrow*kmt(i,jrow)
        enddo
      enddo
      print *, text, cksum
      checksumi = cksum
      return
      end

      subroutine area_volume(kmt, xu, yu, zw)

!-----------------------------------------------------------------------
!     compute surface area and volume of ocean ("t" cells and "u" cells)
!     (note that areas are defined at each level)
!-----------------------------------------------------------------------

#include "size.h"
#include "scalar.h"
      dimension kmt(imt, jmt), xu(imt), yu(jmt), zw(km)
      dimension tcella(km)
      dimension cst(jmt)
      dimension dxt(imt), dyt(jmt)

      do k=1,km
        tcella(k) = 0
      enddo
      ocnp   = 0
      tcellv = 0
      pi = 4.0 * atan(1.0)
      degrad = pi / 180.0
      do i=2,imt-1
        dxt(i) = radius * degrad * (xu(i) - xu(i-1))
      enddo
      do jrow = 2,jmt-1
        cst(jrow) = cos(0.5*(yu(jrow)+yu(jrow-1))*degrad)
        dyt(jrow) = radius * degrad * (yu(jrow) - yu(jrow-1))
      enddo

      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) +
     &           cst(jrow)*dxt(i)*dyt(jrow)
            enddo
            tcellv = tcellv + cst(jrow)*dxt(i)*dyt(jrow)*zw(kmt(i,jrow))
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
        enddo
      enddo

      print '(/)'
      print '(a,f8.0)',   'number of ocean T cells =', ocnp
      print '(a,e15.6,a)','surface area (T cells)  =',tcella(1),' cm**2'
      print '(a,e15.6,a)', 'ocean volume (T cells)  =', tcellv,' cm**3'
      print '(/)'

      return
      end

      subroutine kmtbc (kmt, imt, jmt)

!-----------------------------------------------------------------------
!     set lateral boundary conditions on kmt
!     set cyclic, solid wall and symmetry conditions on "t" grid
!-----------------------------------------------------------------------

      dimension kmt(imt,jmt)

#if defined cyclic

!     set cyclic conditions on eastern and western boundary

      do jrow=1,jmt
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
      enddo
#else

!     set solid wall conditions on eastern and western boundary
!     (if "cyclic" is not enabled then option "solid_walls" is assumed)

      do jrow=1,jmt
# if defined obc_west
        kmt(1,jrow)   = kmt(2,jrow)
# else
        kmt(1,jrow)   = 0
# endif
# if defined obc_east
        kmt(imt,jrow) = kmt(imt-1,jrow)
# else
        kmt(imt,jrow) = 0
# endif
      enddo
#endif
#if defined symmetry

!     set symmetry conditions at "u" row jmt-1

      do i=1,imt
        kmt(i,jmt) = kmt(i,jmt-1)
      enddo
#else
# if defined obc_north

!     extrapolate the domain to the northernmost row

      do i=1,imt
        kmt(i,jmt) = kmt(i,jmt-1)
      enddo
# else

!     close the domain at the northernmost row

      do i=1,imt
        kmt(i,jmt) = 0
      enddo
# endif
#endif
# if defined obc_south

!     extrapolate the domain to the southernmost row

      do i=1,imt
        kmt(i,1)   = kmt(i,2)
      enddo
# else

!     close the domain at the southernmost row

      do i=1,imt
        kmt(i,1)   = 0
      enddo
# endif
      return
      end

#if defined scripps_kmt
      subroutine scripp (kmt, xt, yt, zt, xu, yu, zw
     &,                  imt2, jmt2, km2)

!-----------------------------------------------------------------------
!     construct "kmt" from scripps 1 deg topography.
!-----------------------------------------------------------------------

# include "stdunits.h"
      parameter (iw=362, jw=182)
      parameter (lenw=20*iw)
# include "size.h"

      dimension kmt(imt2,jmt2), xt(imt2), yt(jmt2), zt(km2)
      dimension xu(imt2), yu(jmt2), zw(km2)
      dimension a(imt,jmt)

!     file 'scripps.top' has real *4 ieee format data

      real*4 lons4, lats4, cx4, cy4, d4
      dimension d(iw,jw), cx(iw), cy(jw), dt(iw), cxt(iw), work(lenw)
      dimension d4(iw,jw), cx4(iw), cy4(jw)

!     check that grid sizes in argument list and file "size.h" agree

      call size_check (imt2, jmt2, km2, 'scripps', 'stop')

!-----------------------------------------------------------------------
!     read  Scripps 1 deg topography data
!-----------------------------------------------------------------------

      call getunit (io_scripps, 'scripps.top'
     &,             'unformatted sequential rewind ieee')
      read (io_scripps)
      read (io_scripps) lons4, lats4, cx4, cy4, d4
      call relunit (io_scripps)
      lons = lons4
      lats = lats4
      do i=1,iw
        cx(i) = cx4(i)
      enddo
      do j=1,jw
        cy(j) = cy4(j)
      enddo
      do i=1,iw
        do j=1,jw
          d(i,j) = d4(i,j)
        enddo
      enddo
      write (stdout,'(/,a,/)')
     & ' => Preparing model "kmt" from Scripps 1 deg Topography'

!-----------------------------------------------------------------------
!     translate Scripps grid longitudes to eliminate non-monotonic
!     coordinates (359.5, 0.5, 1.5) across prime meridian from being in
!     the interior of the model grid. (only for limited domain grids
!     that contain the prime meridian)
!-----------------------------------------------------------------------

      call tranlon (d, iw, iw-2, jw, dt, cx, xt(2), imt-2, cxt)

#if !defined skip_translation_details
      write (stdout,'(/a/)')
     & 'Use -Dskip_translation_details to not show the following'
      write(stdout,'(///,50x, a30,/)') 'Scripps 1 deg topography grid'
      write (stdout,'(//,1x,a30,/)') '            latitudes:'
      write (stdout,'(1x,10f10.5)') cy
      write (stdout,'(//,1x,a30,/)') '  original longitudes:'
      write (stdout,'(1x,10f10.5)') cx
      write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
      write (stdout,'(1x,10f10.5)') cxt
#else
      write (stdout,'(/a/)')
     & 'Remove -Dskip_translation_details to show translation details'
#endif

!-----------------------------------------------------------------------
!     interpolate Scripps topography to MOM "t" grid points
!     intrp = 1 => model grid is coarser than 1 deg Scripps data
!     intrp = 2 => model grid is finer than 1 deg Scripps data
!-----------------------------------------------------------------------

!     USER INPUT: set "intrp" appropriately if this test fails

      avgdx = (xu(imt-1) - xu(1))/(imt-2)
      avgdy = (yu(jmt-1) - yu(1))/(jmt-2)
      if (avgdx .gt. 1.0 .and. avgdy .gt. 1.0) then
        intrp = 1
      else
        intrp = 2
      endif
!     end USER INPUT
      if (intrp .eq. 1) then
        print '(/a/a/)'
     &,' "intrp=1" assumes MOM resolution is coarser than Scripps 1 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ftc (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      elseif (intrp .eq. 2) then
        print '(/a/a/)'
     &,' "intrp=2" assumes MOM resolution is finer than Scripps 1 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ctf (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      endif

      do jrow=1,jmt
# if defined cyclic
        a(1,jrow)   = a(imt-1,jrow)
        a(imt,jrow) = a(2,jrow)
# else
        a(1,jrow)   = a(2,jrow)
        a(imt,jrow) = a(imt-1,jrow)
# endif
      enddo

!-----------------------------------------------------------------------
!     convert depths from meters to cm and change sign.
!-----------------------------------------------------------------------

      do jrow=1,jmt
        do i=1,imt
          if (a(i,jrow) .lt. 0.0) then
            a(i,jrow) = -a(i,jrow)*100.0
          else
            a(i,jrow) = 0.0
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     set topographic modifications on model grid here
!-----------------------------------------------------------------------

!     set artificial walls at northern & southern boundaries

      do i=1,imt
        a(i,1 )  = 0.0
        a(i,jmt) = 0.0
      enddo

!-----------------------------------------------------------------------
!     discretize topography "kmt" to nearest "zw" levels
!-----------------------------------------------------------------------

      do jrow=1,jmt
        do i=2,imt-1
          kmt(i,jrow) = 0
          if (a(i,jrow) .ne. 0.0) then
            kmt(i,jrow) = indp (a(i,jrow), zw, km)
          endif
        enddo
# if defined cyclic
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
# else
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
# endif
      enddo

      return
      end
#endif

#if defined etopo_kmt
      subroutine etopo (kmt, xt, yt, zt, xu, yu, zw
     &,                  imt2, jmt2, km2)

!-----------------------------------------------------------------------
!     construct "kmt" from NGDC ETOPO5 1/12 deg topography.
!     This dataset is not available from GFDL.
!     The ETOPO5 dataset can be purchased from the
!     Marine Geology and Geophysics Division of the
!     National Geophysical Data Center.

!     based on code by: Bill Hurlin
!-----------------------------------------------------------------------

#include "stdunits.h"
      parameter (iw=4322, jw=2160, iwm2=iw-2)
      parameter (lenw=20*iw)
#include "size.h"

      dimension kmt(imt2,jmt2), xt(imt2), yt(jmt2), zt(km2)
      dimension xu(imt2), yu(jmt2), zw(km2)
      dimension a(imt,jmt)

!     file 'etopo' has real*4 ieee format data

      dimension lon(iwm2)
      dimension d(iw,jw), cx(iw), cy(jw), dt(iw), cxt(iw), work(lenw)

!     check that grid sizes in argument list and file "size.h" agree

      call size_check (imt2, jmt2, km2, 'etopo', 'stop')

!-----------------------------------------------------------------------
!     read  etopo 1/12 deg topography
!-----------------------------------------------------------------------

      call getunit (io_etopo, '/net/rcp/ETOPO5.NGDCunformat_ieee'
     &,             'unformatted sequential rewind ieee')

      twelvdg = 1.0/12.0
      cx(1)=-twelvdg
      do i = 2,iw
        cx(i) = cx(i-1) + twelvdg
      enddo

      print *,' Reading Etopo 1/12 deg data...'
      rlat = 90.0
      do j=1,jw
        jj = jw - j + 1
        read (io_etopo) lon
        cy(jj) = rlat
        rlat = rlat - twelvdg
        do i=1,iwm2
          d(i+1,jj) = lon(i)
        enddo
        d(1,jj)  = d(iw-1,jj)
        d(iw,jj) = d(2,jj)
!#define debug_etopo
# if defined debug_etopo
        write(*,105) jj, cy(jj), (d(ii,jj),ii=1,2)
     &,              (d(ii,jj),ii=iw-1,iw)
 105    format(1x,'jj=',i5,' cy(jj)=',f8.2
     &,         ' dpt(1..2)=',2f10.4,' dpt(iw-1..iw)=',2f10.4)
# endif
      enddo
      call relunit (io_etopo)
      write (stdout,'(/,a,/)')
     & ' => Preparing model "kmt" from etopo 1/12 deg Topography'

!-----------------------------------------------------------------------
!     translate Etopo grid longitudes to eliminate non-monotonic
!     coordinates (360-1/12, 0, 1/12) across prime meridian from being
!     in the interior of the model grid. (only for limited domain grids
!     that contain the prime meridian)
!-----------------------------------------------------------------------

      call tranlon (d, iw, iw-2, jw, dt, cx, xt(2), imt-2, cxt)

      write(stdout,'(///,50x, a30,/)') 'Etopo 1/12 deg topography grid'
      write (stdout,'(//,1x,a30,/)') '            latitudes:'
      write (stdout,'(1x,10f10.5)') cy
      write (stdout,'(//,1x,a30,/)') '  original longitudes:'
      write (stdout,'(1x,10f10.5)') cx
      write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
      write (stdout,'(1x,10f10.5)') cxt

!-----------------------------------------------------------------------
!     interpolate Scripps topography to MOM "t" grid points
!     intrp = 1 => model grid is coarser than 1 deg Scripps data
!     intrp = 2 => model grid is finer than 1 deg Scripps data
!-----------------------------------------------------------------------

!     USER INPUT: set "intrp" appropriately if this test fails

      avgdx = (xu(imt-1) - xu(1))/(imt-2)
      avgdy = (yu(jmt-1) - yu(1))/(jmt-2)
      if (avgdx .ge. twelvdg .and. avgdy .ge. twelvdg) then
        intrp = 1
      else
        intrp = 2
      endif
      intrp=1
!     end USER INPUT
      if (intrp .eq. 1) then
        print '(/a/a/)'
     &,'"intrp=1" assumes MOM resolution is coarser than Etopo 1/12 deg'
     &,'if not ... then reset it in the USER INPUT section of topog.F'
        call ftc (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      elseif (intrp .eq. 2) then
        print '(/a/a/)'
     &,'"intrp=2" assumes MOM resolution is finer than Etopo 1/12 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ctf (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      endif

      do jrow=1,jmt
#if defined cyclic
        a(1,jrow)   = a(imt-1,jrow)
        a(imt,jrow) = a(2,jrow)
#else
        a(1,jrow)   = a(2,jrow)
        a(imt,jrow) = a(imt-1,jrow)
#endif
      enddo

!-----------------------------------------------------------------------
!     convert depths from meters to cm and change sign.
!-----------------------------------------------------------------------

      do jrow=1,jmt
        do i=1,imt
          if (a(i,jrow) .lt. 0.0) then
            a(i,jrow) = -a(i,jrow)*100.0
          else
            a(i,jrow) = 0.0
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     set topographic modifications on model grid here
!-----------------------------------------------------------------------

!     set artificial walls at northern & southern boundaries

      do i=1,imt
        a(i,1 )  = 0.0
        a(i,jmt) = 0.0
      enddo

!-----------------------------------------------------------------------
!     discretize topography "kmt" to nearest "zw" levels
!-----------------------------------------------------------------------

      do jrow=1,jmt
        do i=2,imt-1
          kmt(i,jrow) = 0
          if (a(i,jrow) .ne. 0.0) then
            kmt(i,jrow) = indp (a(i,jrow), zw, km)
          endif
        enddo
#if defined cyclic
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
#else
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
#endif
      enddo

      return
      end

#endif

#if defined idealized_kmt
      subroutine idealized (kmt, xt, yt, zw, imt, jmt, km)

!-----------------------------------------------------------------------
!     construct a highly "idealized" world ... piece by piece

!     note: the purpose of this geometry/topography is to automatically
!           map into arbitrary resolution as grid dimensions "imt" and
!           "jmt" are changed, thereby facilitating the implementation
!           and verification of the model on various computer platforms
!           without referencing the topographic data base.  Although it
!           somewhat resembles the real world, it is NOT realistic.

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      dimension kmt(imt,jmt), xt(imt), yt(jmt), zw(km)

      do jrow=2,jmt-1
        do i=2,imt-1
          kmt(i,jrow) = km
        enddo
      enddo

      do jrow=1,jmt
        kmt(imt,jrow) = kmt(2,jrow)
        kmt(1,jrow)   = kmt(imt-1,jrow)
      enddo
      do i=1,imt
        kmt(i,jmt) = 0
        kmt(i,1)   = 0
      enddo

!     antarctica

      call setkmt (kmt, xt, yt, imt, jmt, -90.0, 0.0, 360.0, -80.0, 0.0
     &,            360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 360.0-25.0, 360.0
     &,            -70.0, 360.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 0.0, 360.0, -70.0, 0.0
     &,            170.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 360.0-135.0, 360.0-60.0
     &,            -68.0, 360.0-75.0, 360.0-60.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -70.0, 0.0, 155.0, -67.0, 50.0
     &,            145.0, 0)

!     australia

      call setkmt (kmt, xt, yt, imt, jmt, -35.0, 116.0, 120.0, -31.0
     &,            114.0, 130.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -38.0, 140.0, 151.0, -31.0
     &,            130.0, 151.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -31.0, 115.0, 153.0, -20.0
     &,            113.0, 149.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 113.0, 149.0, -11.0
     &,            131.0, 143.0, 0)

!     south america

      call setkmt (kmt, xt, yt, imt, jmt, -50.0, 360.0-74.0, 360.0-68.0
     &,            -40.0, 360.0-73.0, 360.0-62.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -40.0, 360.0-73.0, 360.0-62.0
     &,            -20.0, 360.0-70.0, 360.0-40.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 360.0-70.0, 360.0-40.0
     &,            -16.0, 360.0-81.0, 360.0-35.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -16.0, 360.0-81.0, 360.0-35.0
     &,            0.0, 360.0-80.0, 360.0-50.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 0.0, 360.0-80.0, 360.0-50.0
     &,            11.0, 360.0-75.0, 360.0-60.0, 0)

!     central america

      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 360.0-78.0, 360.0-75.0
     &,            20.0, 360.0-105.0, 360.0-97.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 20.0, 360.0-105.0, 360.0-97.0
     &,            30.0, 360.0-115.0, 360.0-94.0, 0)

!     north america

      call setkmt (kmt, xt, yt, imt, jmt, 25.0, 360.0-82.0, 360.0-80.0
     &,            30.0, 360.0-85.0, 360.0-81.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 30.0, 360.0-115.0, 360.0-80.0
     &,            40.0, 360.0-124.0, 360.0-74.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 360.0-124.0, 360.0-74.0
     &,            50.0, 360.0-124.0, 360.0-57.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 50.0, 360.0-124.0, 360.0-57.0
     &,            60.0, 360.0-140.0, 360.0-64.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 360.0-165.0, 360.0-64.0
     &,            65.0, 360.0-140.0, 360.0-64.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 65.0, 360.0-140.0, 360.0-64.0
     &,            70.0, 360.0-162.0, 360.0-72.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 360.0-162.0, 360.0-140.0
     &,            72.0, 360.0-157.0, 360.0-157.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 360.0-130.0, 360.0-70.0
     &,            75.0, 360.0-120.0, 360.0-80.0, 0)

!     greenland

      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 360.0-45.0, 360.0-45.0
     &,            75.0, 360.0-58.0, 360.0-19.0, 0)

!     africa

      call setkmt (kmt, xt, yt, imt, jmt, -35.0, 19.0, 28.0, 6.0, 8.0
     &,            50.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 0.0, 50.0, 18.0, 0.0
     &,            56.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 18.0, 0.0, 56.0, 26.0, 0.0
     &,            59.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 360.0-10.0, 360.0, 18.0
     &,            360.0-18.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 18.0, 360.0-18.0, 360.0, 26.0
     &,            360.0-15.0, 360.0, 0)

!     northern africa &  europe & asia

      call setkmt (kmt, xt, yt, imt, jmt, 26.0, 360.0-15.0, 360.0, 40.0
     &,            360.0-7.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 360.0-7.0, 360.0, 50.0
     &,            360.0, 360.0, 0)

      call setkmt (kmt, xt, yt, imt, jmt, 8.0, 77.0, 78.0, 26.0, 65.0
     &,            90.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 4.0, 99.0, 100.0, 26.0, 90.0
     &,            115.0, 0)

      call setkmt (kmt, xt, yt, imt, jmt, 26.0, 0.0, 126.0, 40.0, 0.0
     &,            122.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 0.0, 130.0, 50.0, 0.0
     &,            140.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 50.0, 0.0, 140.0, 60.0, 8.0
     &,            140.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 8.0, 163.0, 65.0, 13.0
     &,            180.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 65.0, 13.0, 188.0, 70.0, 20.0
     &,            180.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 70.0, 180.0, 75.0, 90.0
     &,            100.0, 0)

!     add an "idealized" undulating topography

      bot = zw(km)
      pi  = 4.0*atan(1.0)
      do jrow=2,jmt
        do i=2,imt-1
          if (kmt(i,jrow) .ne. 0) then
            arg = bot*(1-0.4*abs(cos((jrow*pi)/jmt)*sin((i*2*pi)/imt)))
            kmt(i,jrow) = indp (arg, zw, km)
          endif
        enddo
      enddo

!     add "idealized" ridges

      level = indp (0.666*zw(km), zw, km)

      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 360.0-20.0, 360.0-10.0
     &,            30.0, 360.0-45.0, 360.0-35.0, level)
      call setkmt (kmt, xt, yt, imt, jmt, 30.0, 360.0-45.0, 360.0-35.0
     &,            60.0, 360.0-20.0,  360.0-30.0, level)
      call setkmt (kmt, xt, yt, imt, jmt, -60.0,360.0-100.0, 360.0-130.0
     &,            40.0, 360.0-160.0, 180.0, level)

      level = indp (0.5*zw(km), zw, km)

      call setkmt (kmt, xt, yt, imt, jmt, -50.0, 360.0-120.0
     &,            360.0-120.0, 30.0, 190.0, 190.0, level)
      return
      end
#endif

      subroutine setkmt (kmt, xt, yt, imt, jmt, alat1, slon1, elon1
     &,                  alat2, slon2, elon2, num)

!-----------------------------------------------------------------------
!     set the topography mask "kmt(i,j)" = "num" within the area of
!     the parallelogram bounded by vertices:
!     (alat1,slon1), (alat1,elon1), (alat2,slon2), & (alat2,elon2)

!     inputs:

!     xt = longitudes of "t" points in degrees
!     yt = latitudes of "t" points in degrees
!     imt = number of model longitudes
!     jmt = number of model latitudes
!     alat1 = southern latitude of parallelogram (degrees)
!     slon1 = starting longitude of southern edge of parallelogram (deg)
!     elon1 = ending longitude of southern edge of parallelogram (deg)
!     alat2 = northern latitude of parallelogram (degrees)
!     slon2 = starting longitude of notthern edge of parallelogram (deg)
!     elon2 = ending longitude of notthern edge of parallelogram (deg)
!     num   = number of vertical levels

!     outputs:

!     kmt   = mask of vertical levels on model "t" points

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      dimension kmt(imt,jmt), xt(imt), yt(jmt)

!     convert the four vertices into model indices
!     (js,is1), (js,ie1), (je,is2), (je,ie2)

      j1 = indp (alat1, yt, jmt)
      j2 = indp (alat2, yt, jmt)
      js = min (j1,j2)
      je = max (j1,j2)

      i1  = indp (slon1, xt, imt)
      i2  = indp (elon1, xt, imt)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)

      i1  = indp (slon2, xt, imt)
      i2  = indp (elon2, xt, imt)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)

      is = is1
      ie = ie1

!     fill in the area bounded by (js,is1), (js,ie1), (je,is2), (je,ie2)
!     the nudging of 1.e-5 is to insure that the test case resolution
!     generates the same topography and geometry on various computers.

      c1 = 1.0
      if (js .eq. je) then
        rdj = c1
      else
        rdj = c1/(je-js)
      endif
      do jrow=js,je
        do i=is,ie
          kmt(i,jrow) = num
        enddo
        is = nint(rdj*((jrow-js)*is2 + (je-jrow)*is1) + 1.0e-5)
        ie = nint(rdj*((jrow-js)*ie2 + (je-jrow)*ie1) + 1.0e-5)
      enddo
      return
      end

      subroutine min_depth (kmt, kmt_min, xu, yu, zw)

!     limit the minimum number of levels. kmt_min should be >= 2

#include "size.h"
#include "stdunits.h"

      dimension kmt(imt,jmt)
      dimension xu(imt), yu(jmt), zw(km)
      write (stdout,'(/a/)')
     & 'Searching for and correcting minimum level violations'
      n = 0
      do i=2,imt-1
        do jrow=jmt-1,2,-1
          if (kmt(i,jrow) .ne. 0 .and. kmt(i,jrow) .lt. kmt_min) then
            n = n + 1
            kmt_shallow = kmt(i,jrow)
#if defined round_shallow || !(defined deepen_shallow || defined fill_shallow)
            if (zw(kmt(i,jrow)) .lt. 0.5*zw(kmt_min)) then
              kmt(i,jrow) = 0
            else
              kmt(i,jrow) = kmt_min
            endif
#endif
#if defined fill_shallow
            kmt(i,jrow) = 0
#endif
#if defined deepen_shallow
            kmt(i,jrow) = kmt_min
#endif
#if defined debug_topog
            write (stdout,'(a,i4,a,i4,a,i4,a,i4)')
     &      'a) changed shallow kmt(',i,',',jrow,') from ', kmt_shallow
     &,     ' to ', kmt(i,jrow)
#endif
          endif
        enddo
      enddo
      if (n .gt. 0) then
        write (stdout,'(a,i5,a/)')  '->Modified', n,' shallow cells'
      else
        write (stdout,'(a/)') '->No modifications needed'
      endif
      return
      end
