      subroutine gcoord (maxlen, imt, jmt, km, dxtdeg, dytdeg, dxudeg
     &,                  dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)

!=======================================================================

!                      G R I D    C O O R D I N A T E S

!     Construct grid point coordinates and resolution

!     input:

!     maxlen = maximum number of grid cells in latitude, longitude,
!              and depth

!              set grid specifications in USER INPUT section.

!     output:

!     imt    = number of longitudes
!     jmt    = number of latitudes
!     km     = number of depths
!     dxtdeg = width of "t" grid cells (degrees)
!     dytdeg = height of "t" grid cells (degrees)
!     dxudeg = width of "u" grid cells (degrees)
!     dyudeg = height of "u" grid cells (degrees)
!     dzt    = thickness of "t" grid cells (cm)
!     dzw    = thickness of "w" grid cells (cm)
!     xt     = longitude at centres of "t" grid cells (degrees)
!     xu     = longitude at centres of "u" grid cells (degrees)
!     yt     = latitude at centres of "t" grid cells (degrees)
!     yu     = latitude at centres of "u" grid cells (degrees)
!     zt     = depth at centres of "t" grid cells (centimetres)
!     zw     = depth at centres of "u" grid cells (centimetres)

!     based on code by: R. C. Pacanowski
!=======================================================================

#include "stdunits.h"
      parameter (maxbounds=11)
      dimension xt(imt), yt(jmt), xu(imt), yu(jmt)
      dimension zw(km), zt(km)

      dimension dxtdeg(imt), dytdeg(jmt), dzt(km)
      dimension dxudeg(imt), dyudeg(jmt), dzw(0:km)

      dimension x_lon(maxbounds),   dx_lon(maxbounds)
      dimension y_lat(maxbounds),   dy_lat(maxbounds)
      dimension z_depth(maxbounds), dz_depth(maxbounds)

      integer ib(10), ic(10)
      character(120) :: fname, new_file_name
      logical exists

!-----------------------------------------------------------------------
!     USER INPUT ===>  read this section and specify the grid.
!-----------------------------------------------------------------------

!     How to define a grid resolution:

!     The ocean model domain is composed of one or more regions in
!     latitude, longitude, and depth. Latitude and longitude are in
!     degrees east of Greenwich and depth is in centimetres. Each
!     region is defined by its bounds and resolution at those bounds.
!     Within each region, resolution may be constant or smoothly
!     varying but there must be an integral number of grid cells
!     contained within the region`s bounds.

!     If the resolution at both bounds of a region is the same, then
!     resolution within the region is constant. If the bounding
!     resolutions differ, resolution varies smoothly across the region
!     to make the transition. The functional for the variation is
!     taken to be a cosine. Regions sharing a common boundary have the
!     same resolution at that boundary AND the first derivative is
!     zero there to minimize effects of numerics. In the vertical, the
!     last region allows a stretching factor to provide a more drastic
!     fall off of resolution to the bottom if desired.

!     Example:

!     Define a domain composed of two regions of longitude: The first
!     is bounded by longitudes 120E and 180E with constant 1.0 degree
!     resolution. The second extends east of 180E to 300E where the
!     resolution is desired to be 3 degrees (at 300E). To do this,
!     specify the following:

!      parameter (nxlons=3)
!      data (x_lon(i), i=1,nxlons)  /120.0, 180.0, 300.0/
!      data (dx_lon(i),i=1,nxlons)  /  1.0,   1.0,   3.0/

!      where...

!     "nxlons"  = number of bounding longitudes to define regions {3}
!     "x_lon"   = bounding longitudes {120.0E, 180.0E, 300.0E}
!     "dx_lon"  = resolution centred at "x_lon" {1.0, 1.0, 3.0}

!     Note: Region bounds are not defined exactly as in MOM 1. Now they
!           are defined on "u" grid points (which are at the edges of
!           "T" cells) as described below.

!     The western edge of the 1st "T" cell will be at x_lon(1) and
!     the eastern edge of the last "T" cell within the first region
!     will be at x_lon(2). Since the resolution at both bounds is
!     1.0 degree, there will be a constant 1 degree resolution within
!     this first region.

!     The western edge of the 1st "T" cell within the second region is
!     at x_lon(2) and the eastern edge of the last "T" cell is at
!     x_lon(3). Since resolutions at the bounds of region 2 differ,
!     resolution will gradually taper from 1.0 degree at the western
!     bound to 3.0 degrees at the eastern bound. Note that the
!     resolution of the 1st "T" cell adjacent to x_lon(2) will
!     be slightly greater than 1.0 degree (it would be 1.0 degree if
!     the grid cell were centred exactly at x_lon(2) but it is not).

!     The number of "T" cells within each region will be equal to the
!     width of the region divided by the average resolution within the
!     region. The average resolution within the region is simply the
!     average of the bounding resolutions. In region #1, it will be
!     (1.0 + 1.0)/2 = 1.0 degree and in region #2 it will be
!     (1.0 + 3.0)/2 = 2.0 degree. Since region #1 is 60 degree wide,
!     the number of "T" cells will be 60.0/1.0 = 60 and for region #2
!     there will be 120.0/2.0 = 60 cells also.

!     How are "u" cells defined?

!     Since MOM uses an Arakawa "B" grid, the "u" cells are offset from
!     the "T" cells. The "u" cell grid points are at the vertices of the
!     "T" cells and visa versa. The "u" cell with coordinates (i,j) is
!     on the northeast vertex of the "T" cell with coordinates (i,j). In
!     the vertical, "u" cells and "T" cells are at the same level.
!     However, advective vertical velocity is defined at the base of
!     "T" and "u" cells. Thus "w" cells can be thought of as vertical
!     advective velocity cells with the grid being staggered in the
!     vertical. The relationship between "T" (or "u") cells and "w"
!     cells is that the "w" cell with index (k) is at the base of the
!     "T" (or "u") cell with index k.

!     What happens when the grid cells are stretched?

!     Grid stretching (variable resolution) is not done exactly as in
!     MOM 1. In MOM 1, the "u" cell resolution was the average of the
!     resolutions of the adjacent "T" cells. For instance: "dxt" was
!     specified first... then dxu(i) = 0.5*(dxt(i+1) + dxt(i)).
!     In a stretched grid, the "T" grid points were always in the middle
!     of "T" cells and the "u" grid points were off centre in "u" cells.

!     In MOM 2, it`s the other way around: the "T" cell resolution is
!     the average of the resolutions of the adjacent "u" cells. So,
!     "dxu" is specified first... then dxt(i) = 0.5*(dxu(i-1) + dxu(i))
!     Also, "u" grid points are always in the middle of "u" cells and
!     "T" grid points are off centre in "T" cells.

!     This allows for more accurate advection velocities to be defined
!     on the faces of the cells when the resolution varies. The
!     implication is that advection will be more accurate for tracers
!     in the horizontal and both tracers and momentum in the vertical.
!     Advection will not be more accurate for momentum in the horizontal
!     because "u" grid  points will always be in the centre of the "u"
!     cells.

!     What if I want the resolution from MOM 1?

!     It is possible to closely approximate the horizontal resolution
!     from MOM 1 by enabling the "-Dcentered_t" option. Whether this is
!     done or not, the "u" grid points are always on the edges of "T"
!     cells and the "T" grid points are always on the edges of "u"
!     cells. In the vertical, MOM 1 used a Gaussian for stretching the
!     resolution. MOM 2 uses the same cosine functional for varying
!     resolution in the vertical as in the horizontal. Additionally,
!     a stretching factor is used to give the functionality of enhanced
!     stretching toward the bottom of the ocean domain if desired.

!     Another alternative is to read in whatever grid discretization you
!     like by enabling the -D"read_my_grid" option and being aware of
!     the potential consequences. The crucial point is that if you have
!     enough resolution to adequately resolve the scales of motion (but
!     who does?), it does not much matter how the resolution is
!     distributed. However, keep in mind the above remarks about
!     advection.

!-----------------------------------------------------------------------
!     The following specifies the global test grid resolution:
!-----------------------------------------------------------------------

!     constant 4.0 degree resolution in longitude

      parameter (nxlons=2)
      data (x_lon(i), i=1,nxlons) / 0.0, 360.0/
      data (dx_lon(i),i=1,nxlons) / 4.0,   4.0/

!     constant 3.0 degree resolution in latitude

      parameter (nylats=2)
      data (y_lat(j), j=1,nylats) / -87.0,  90.0/
      data (dy_lat(j),j=1,nylats) /  3.0,    3.0/

!     two vertical regions: constant 25m in the upper 100m, and variable
!     resolution in the lower region. Note that by setting the stretch
!     factor "stretch_z" > 1.0, the last region in the vertical can
!     have the grid cells stretched further. To see the effect of
!     stretching on vertical cell thickness and number of vertical
!     cells, try using "stretch_z"=1.1 and gradually increase it.

      parameter (nzdepths=3)
      data (z_depth(k), k=1,nzdepths) /  0.0e2,   100.0e2,  5600.0e2/
      data (dz_depth(k),k=1,nzdepths) / 25.0e2,    25.0e2,   975.0e2/
      data stretch_z /1.0/

!     Some other examples:

!      4.0 degrees at the poles down to 2.0 deg at the equator

!      parameter (nylats=3)
!      data (y_lat(j), j=1,nylats) / -90.0,  0.0,  90.0/
!      data (dy_lat(j),j=1,nylats) /   4.0,  2.0,   4.0/

!      set 1 deg at -30 deg lat down to 1/3 deg at -10 deg lat
!      constant 1/3 deg resolution from -10 deg to +10 deg lat and
!      1/3 deg at +10 deg lat to 1.0 deg at +30 deg lat

!      parameter (nylats=4)
!      data (y_lat(j), j=1,nylats) / -30.0,  -10.0,    10.0,  30.0/
!      data (dy_lat(j),j=1,nylats) /   1.0,  0.3333,  0.3333,  1.0/

# if defined obctest
!-----------------------------------------------------------------------
!     The following specifies the box for the open boundary test:
!     see Stevens (1990), Chap. 5: A Passive Open Boundary Condition
!-----------------------------------------------------------------------

!     constant 1.0 degree resolution

      parameter (nxlons=2)
      data (x_lon(i), i=1,nxlons) / -0.5,  18.5/
      data (dx_lon(i),i=1,nxlons) /  1.0,   1.0/

!     constant 0.5 degree resolution

      parameter (nylats=2)
# if defined obc_south && defined obc_north
      data (y_lat(j), j=1,nylats)  / 56.75,  67.75/
# else
      data (y_lat(j), j=1,nylats)  / 52.75,  72.25/
# endif
      data (dy_lat(j),j=1,nylats) /   0.5,    0.5 /

      parameter (nzdepths=2)
      data (z_depth(k), k=1,nzdepths) /   0.0e2,  2700.0e2/
      data (dz_depth(k),k=1,nzdepths) /  50.0e2,   850.0e2/
      data stretch_z /1.0/

# endif
# if defined obctest2
!-----------------------------------------------------------------------
!     The following specifies the box for the second open boundary test:
!     see Stevens (1990), 6. An Active Open Boundary Condition
!-----------------------------------------------------------------------

!     constant 2.0 degree resolution

      parameter (nxlons=2)
      data (x_lon(i), i=1,nxlons) / -1.0,  37.0/
      data (dx_lon(i),i=1,nxlons) /  2.0,   2.0/

!     constant 2.0 degree resolution

      parameter (nylats=2)
# if defined obc_south
      data (y_lat(j), j=1,nylats)  / 36.00,  54.00/
# else
      data (y_lat(j), j=1,nylats)  / 16.00,  54.00/
# endif
      data (dy_lat(j),j=1,nylats) /   2.0,    2.0 /

      parameter (nzdepths=2)
      data (z_depth(k), k=1,nzdepths) /   0.0e2,  3750.0e2/
      data (dz_depth(k),k=1,nzdepths) /  50.0e2,  1200.0e2/
      data stretch_z /1.0/

# endif
!-----------------------------------------------------------------------
!     ==> end of USER INPUT
!-----------------------------------------------------------------------

!     set some constants

      p5 = 0.5

      ncase = 0
#if defined read_my_grid
      ncase = ncase + 1
      c0 = 0.
      c1 = 1.
      c100 = 100.
      fname = new_file_name ("grid.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'gcoord in grids.f'
      endif
      dxtdeg(:) = c0
      dytdeg(:) = c0
      dxudeg(:) = c0
      dyudeg(:) = c0
      dzt(:) = c0
      dzw(:) = c0
      xt(:) = c0
      yt(:) = c0
      xu(:) = c0
      yu(:) = c0
      zt(:) = c0
      zw(:) = c0
      ib(:) = 1
      call openfile (fname, iou)
      ic(:) = imt
      call getvara ('dxt', iou, imt, ib, ic, dxtdeg, c1, c0)
      call getvara ('dxu', iou, imt, ib, ic, dxudeg, c1, c0)
      call getvara ('xt', iou, imt, ib, ic, xt, c1, c0)
      call getvara ('xu', iou, imt, ib, ic, xu, c1, c0)
      ic(1) = jmt
      call getvara ('dyt', iou, jmt, ib, ic, dytdeg, c1, c0)
      call getvara ('dyu', iou, jmt, ib, ic, dyudeg, c1, c0)
      call getvara ('yt', iou, jmt, ib, ic, yt, c1, c0)
      call getvara ('yu', iou, jmt, ib, ic, yu, c1, c0)
      ic(1) = km
      call getvara ('dzt', iou, km, ib, ic, dzt, c100, c0)
      call getvara ('dzw', iou, km, ib, ic, dzw(1:km), c100, c0)
      call getvara ('zt', iou, km, ib, ic, zt, c100, c0)
      call getvara ('zw', iou, km, ib, ic, zw, c100, c0)
      call closefile (iou)

      if (zt(km) .gt. 8000.e2) then
!       assume depths are in cm
        print*, "==> Warning: converting depths to cm in gcoord"
        zt(1:km) = zt(1:km)/100.
        dzt(1:km) = dzt(1:km)/100.
        zw(1:km) = zw(1:km)/100.
        dzw(1:km) = dzw(1:km)/100.
      endif
      dzw(0)  = zt(1)
      dzw(km) = zw(km) - zt(km)
#endif
#if defined generate_a_grid
      ncase = ncase + 1

      write (stdout,'(//,36x,a,/)') 'G R I D   G E N E R A T I O N'
      write (stdout,'(/,43x,a,/)') 'Grid resolution:'

      if (nxlons .gt. maxbounds .or. nylats .gt. maxbounds .or.
     &    nzdepths .gt. maxbounds) then
        write (stdout,'(/a/)')
     &  ' Increase parameter "maxbounds" to contain all regions.'
        stop
      endif

!-----------------------------------------------------------------------
!     Calculate resolution in longitude. Add one boundary cell at the
!     start and end of the domain so that calculations are meaningful
!     for grid cells i=2,imt-1
!-----------------------------------------------------------------------

      nbpts = 2
      write (stdout,'(/a,i1,a)')
     & ' Generating the longitudinal resolution ( ', nbpts
     &, ' extra boundary points will be used).'

!     convert negative longitudes (degrees west of Greenwich) to
!     positive longitudes

      offset = 0.0
      do n=1,nxlons
        if (x_lon(n) .lt. 0.0) then
          offset = 360.0
        endif
      enddo
      if (offset .gt. 0.0) then
          write (stdout,'(/,a/,a/,a)')
     &    '=>Warning: adding 360.0 degrees to all longitudes to convert'
     &,   '           negative values (degrees west of Greenwich) to'
     &,   '           positive values (degrees east of Greenwich)'
      endif
      do n=1,nxlons
        x_lon(n) = x_lon(n) + offset
      enddo

!     if a region contains the Greenwich meridian, compensate by adding
!     360.0 degrees to all remaining regions

      do n=2,nxlons
        if (x_lon(n-1) .gt. x_lon(n)) then
          x_lon(n) = x_lon(n) + 360.0
          write (stdout,'(/,a,i1,a,a/)')
     &    ' Warning: adding 360.0 degrees to x_lon(',n,')'
     &,   ' to insure region boundaries increase monotonically'
        endif
      enddo

!     if width of domain exceeds 360.0 deg... limit width to 360.0

      if (x_lon(nxlons) - x_lon(1) .gt. 360.0) then
        write (stdout,'(/a,a,g14.7/)')
     &  '=>Warning: Domain width exceeds 360 deg. Restricting last'
     &, ' x_lon to ',x_lon(1)+360.0
          do m=1,nxlons
            write (stdout,'(i3,f10.5)') m, x_lon(m)
          enddo
      endif
      if (x_lon(nxlons) - x_lon(1) .eq. 360.0) then
        if (dx_lon(nxlons) .ne. dx_lon(1)) then
          write (stdout,'(/a,a)')
     &    '=>Error: dx_lon(1) must equal dx_lon(last) when domain'
     &,   ' width = 360.0 degrees'
        endif
      endif

      stretch_x = 1.0
      call gcell (maxlen, nxlons, x_lon,  dx_lon, nbpts
     &,           imt, dxtdeg, dxudeg, stretch_x)

!-----------------------------------------------------------------------
!     Build the longitudinal grid points for a "B" grid
!     (account for an extra boundary point at the start)
!-----------------------------------------------------------------------

      xt(1) = x_lon(1) - p5*dx_lon(1)
      xu(1) = x_lon(1)

      do i=2,imt
        xu(i) = xu(i-1) +  dxtdeg(i)
        xt(i) = xt(i-1) +  dxudeg(i-1)
      enddo

!-----------------------------------------------------------------------
!     Calculate resolution in latitude. Add one boundary cell at the
!     start and end of the domain so that calculations are meaningful
!     for grid cells jrow=2,jmt-1
!-----------------------------------------------------------------------

      nbpts = 2
      write (stdout,'(/a,i1,a)')
     & ' Generating the latitudinal resolution ( ', nbpts
     &, ' extra boundary points will be used).'

!     if width of domain exceeds 180.0 deg... limit width to 180.0

      if (y_lat(nylats) - y_lat(1) .gt. 180.0) then
        write (stdout,'(/a,a,g14.7/)')
     &  '=>Warning: Latitudinal domain width exceeds 180 deg.'
     &, ' Restricting last y_lat to ',y_lat(1)+180.0
      endif
      do n=2,nylats
        if (y_lat(n-1) .gt. y_lat(n)) then
          write (stdout,'(/,a,/a/)')
     &    ' =>Error: latitude boundaries must increase monotonically'
     &,   '          check the specifications in the USER INPUT section'
          do m=1,nylats
            write (stdout,'(i3,f10.5)') m, y_lat(m)
          enddo
          stop '==>grids'
        endif
      enddo
      stretch_y = 1.0
      call gcell (maxlen, nylats, y_lat, dy_lat, nbpts
     &,           jmt, dytdeg, dyudeg, stretch_y)

!-----------------------------------------------------------------------
!     Build the latitudinal grid points on a "B" grid
!     (account for an extra boundary point at the start)
!-----------------------------------------------------------------------

      yt(1) = y_lat(1) - p5*dy_lat(1)
      yu(1) = y_lat(1)
      do jrow=2,jmt
        yu(jrow) = yu(jrow-1) +  dytdeg(jrow)
        yt(jrow) = yt(jrow-1) +  dyudeg(jrow-1)
      enddo

!-----------------------------------------------------------------------
!     Calculate resolution in depth. No boundary cells added here so
!     calculations are meaningful for k=1,km. Allow the bottom region
!     to be stretched further if desired.
!-----------------------------------------------------------------------

      nbpts = 0
      write (stdout,'(/a,i1,a)')
     & ' Generating the vertical resolution ( ', nbpts
     &, ' extra boundary points will be used).'
      do n=2,nzdepths
        if (z_depth(n-1) .gt. z_depth(n)) then
          write (stdout,'(/,a,/a/)')
     &    ' =>Error: depth boundaries must increase monotonically'
     &,   '          check the specifications in the USER INPUT section'
          stop '==>grids'
        endif
      enddo
      call gcell (maxlen, nzdepths, z_depth, dz_depth, nbpts
     &,           km, dzt, dzw(1), stretch_z)

!-----------------------------------------------------------------------
!     Build the vertical grid points on a "B" grid. The "T" and "u"
!     cells are staggered in the horizontal but at the same level in
!     the vertical. However, the "w" cells here refer to the vertical
!     advection velocities at the bottoms of the "u" and "T" cells.
!     (no extra boundary point at the start)
!-----------------------------------------------------------------------

      zt(1) = z_depth(1) + p5*dz_depth(1)
      zw(1) = z_depth(1) + dzt(1)
      do k=2,km
        zw(k) = zw(k-1) + dzt(k)
        zt(k) = zt(k-1) + dzw(k-1)
      enddo

!     set "w" cell thickness at surface and bottom to ocean part of cell

      dzw(0)  = zt(1)
      dzw(km) = zw(km) - zt(km)

!-----------------------------------------------------------------------
!     Print grid "t" cell resolution in longitude, latitude and depth
!     (also for "w" cells in depth)
!-----------------------------------------------------------------------

      write (stdout,9101)
     &'Vertical resolution of cells "dzw(k)" k=0,', km, 'cm'
      write (stdout,9002) (dzw(k),k=0,km)
      write (stdout,9101)
     &'Vertical resolution of cells "dzt(k)" k=1,', km, 'cm'
      write (stdout,9002) (dzt(k),k=1,km)
      write (stdout,9101)
     &'Longitudinal resolution of cells "dxtdeg(i)" i=1,', imt, 'deg'
      write (stdout,9001) (dxtdeg(i),i=1,imt)
      write (stdout,9101)
     &'Latitudinal resolution of cells "dytdeg(j)" jrow=1,', jmt, 'deg'
      write (stdout,9001) (dytdeg(jrow),jrow=1,jmt)
#endif
      if (ncase .eq. 0) then
        write (stdout,'(/a/a/)')
     & '=>Error: One of the following options must be enabled:'
     &,'         generate_a_grid    read_my_grid'
        stop '=>grid'
      elseif (ncase .gt. 1) then
        write (stdout,'(/a/a/)')
     & '=>Error: Only one of the following options may be enabled:'
     &,'         generate_a_grid    read_my_grid'
        stop '=>grid'
      endif

!-----------------------------------------------------------------------
!     Check if the "t" grid resolution is an average of the
!     "u" cell resolution. This insures more accurate advection of
!     tracers within a stretched grid.
!-----------------------------------------------------------------------

      num  = 0
      tolr = 1.e-5
      write (stdout,'(/)')
      do i=2,imt-1
        dxubar = p5*(dxudeg(i) + dxudeg(i-1))
        if (abs(dxubar-dxtdeg(i)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta x at i=',i
     &,   ' is not an average of adjacent "u" cell delta x`s'
        endif
      enddo

      do jrow=2,jmt-1
        dyubar = p5*(dyudeg(jrow) + dyudeg(jrow-1))
        if (abs(dyubar-dytdeg(jrow)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta y at jrow=',jrow
     &,   ' is not an average of adjacent "u" cell delta y`s'
        endif
      enddo

      tolr = 1.e0
      do k=2,km-1
        dzwbar = p5*(dzw(k) + dzw(k-1))
        if (abs(dzwbar-dzt(k)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta z at k=',k
     &,   ' is not an average of adjacent "w" cell delta z`s'
        endif
      enddo

      if (num .ne. 0) then
        write (stdout,'(/a/a/a/a//a,a/)')
     &  '==>Warning, At the above locations, advection of tracers is'
     &, 'not as accurate as it could be. If you are reading in your own'
     &, 'grid or constructing a grid as in MOM 1, we assume you want to'
     &, 'define the grid this way and we let you proceed from here...'
     &, 'Please read ALL the information in the USER INPUT section to '
     &, 'understand what this means'
      endif

!-----------------------------------------------------------------------
!     Print all grid coordinates
!-----------------------------------------------------------------------

      write (stdout
     &,'(//,40x,a,//,a,g14.7,a,/a/,a,g14.7,a/a,/,a,g14.7,a)')
     &  ' Grid Point Coordinate details: '
     &, ' The western edge of the 2nd "t" cell is at longitude:'
     &,   xu(1),' (deg)',' (the 1st "t" cell is a boundary cell)'
     &, ' The southern edge of the 2nd "t" cell is at latitude:'
     &,   yu(1),' (deg)',' (the 1st "t" cell is a boundary cell)'
     &,' The top edge of the 1st "t" cell is at z =',z_depth(1),' (cm)'
      write (stdout,'(/,a,g14.7,a/a/,a,g14.7,a/a/,a,g14.7,a/)')
     &  ' The western edge of the 1st "u" cell is at longitude:', xt(1)
     &, ' (deg)',' (the 1st "u" cell is a boundary point)'
     &, ' The southern edge of the 1st "u" cell is at latitude:', yt(1)
     &, ' (deg)',' (the 1st "u" cell is a boundary point)'
     &, ' The top edge of the 1st "w" cell is at z =',zt(1),' (cm)'
      write (stdout,9103) km
      write (stdout,9002) (zt(k),k=1,km)
      write (stdout,9104) km
      write (stdout,9002) (zw(k),k=1,km)
      write (stdout,9105) jmt
      write (stdout,9001) (yt(jrow),jrow=1,jmt)
      write (stdout,9106) jmt
      write (stdout,9001) (yu(jrow),jrow=1,jmt)
      write (stdout,9107) imt
      write (stdout,9001) (xt(i),i=1,imt)
      write (stdout,9108) imt
      write (stdout,9001) (xu(i),i=1,imt)

#if defined symmetry

!-----------------------------------------------------------------------
!     insure that yu(jmt-1) = 0.0 (equator) when using symmetry
!-----------------------------------------------------------------------

      if (yu(jmt-1) .ne. 0.0) then
        write (stdout,*) '=> Error: yu(jmt-1) must = 0.0 for symmetry'
        stop '=>gcoord'
      endif
#endif

!---------------------------------------------------------------------
!     compute a grid checksum
!---------------------------------------------------------------------

      cksum = 0.0
      cksum = cksum + checksum (xt, imt, 1)
      cksum = cksum + checksum (yt, jmt, 1)
      cksum = cksum + checksum (zt, km, 1)
      cksum = cksum + checksum (xu, imt, 1)
      cksum = cksum + checksum (yu, jmt, 1)
      cksum = cksum + checksum (zw, km, 1)
      cksum = cksum + checksum (dxtdeg, imt, 1)
      cksum = cksum + checksum (dytdeg, jmt, 1)
      cksum = cksum + checksum (dxudeg, imt, 1)
      cksum = cksum + checksum (dyudeg, jmt, 1)
      cksum = cksum + checksum (dzt, km, 1)
      cksum = cksum + checksum (dzw, km+1, 1)
      write (stdout,'(/)')
      write (stdout,*) 'Grid checksum = ',cksum
      write (stdout,'(/)')
      return
9001  format (1x,10f10.4)
9002  format (1x,10f10.2)
#if defined generate_a_grid
9101  format (/,  a,i4,' in units of ',a,' as follows:')
#endif
9103  format (/,' Depth to "t" & "u" grid points (cm): zt(k) k=1,',i3)
9104  format (/,' Depth to "w" grid points (cm): zw(k) k=1,',i3)
9105  format (/,' Latitude of "t" points (deg): yt(j) j=1,',i4)
9106  format (/,' Latitude of "u" points (deg): yu(j) j=1,',i4)
9107  format (/,' Longitude of "t" points (deg): xt(i) i=1,',i4)
9108  format (/,' Longitude of "u" points (deg): xu(i) i=1,',i4)
      end

      subroutine gcell (maxlen, n_bounds, bounds, d_bounds, nbpts
     &,                 num, deltat, deltau, stretch)

!=======================================================================

!              G R I D   C E L L   C O N S T R U C T I O N

!     A domain is composed of one or more regions:
!     Build "num" "t"  cells with resolution "deltat(n) n=1,num"
!     within the domain composed of regions bounded by "bounds".
!     Also construct "num" "u"  cells of resolution "deltau(n) n=1,num"
!     with the relation between "t" and "u" cells given by:
#if defined centered_t
!     deltau(n) = 0.5*(deltat(n+1) + deltat(n))
#else
!     deltat(n) = 0.5*(deltau(n-1) + deltau(n))
#endif
!     Resolution may be constant or smoothly varying within each
!     region AND there must be an integral number of grid cells within
!     each region. The domain is the sum of all regions.

!     inputs:

!     maxlen   = maximum length of "deltat" and "deltau"
!     n_bounds = number of bounds needed to define the regions
!     bounds   = latitude, longitude, or depth at each bound
!     d_bounds = delta (resolution) at each of the "bounds"
!     nbpts    = number of extra boundary cells to add to the domain.
!                (usually one at the beginning and end)
!     stretch  = stretching factor for last region (should only be used
!                in the vertical to provide increased stretching of grid
!                points. "stretch" = 1.0 gives no increased stretching.
!                "stretch" = 1.2 gives increased stretching...etc

!     outputs:

!     num    = total number of grid cells within the domain
!     deltau = resolution of "u" grid cells: n=1,num
!     deltat = resolution of "t" grid cells: n=1,num

!     based on code by: R. C. Pacanowski
!=======================================================================

#include "stdunits.h"
      dimension deltat(maxlen), deltau(maxlen)
      dimension d_bounds(n_bounds), bounds(n_bounds)

!     Set some constants

      p5 = 0.5
      pi = 4.0*atan(1.0)

!     Do all regions, one at a time, to construct the domain

      num  = 1
      do l=1,n_bounds-1

        write (stdout,'(2x,a,i2,a,g14.7,a,g14.7,a,g14.7,a,g14.7,a)')
     & ' region # ',l,'  going from ',bounds(l),' (res=',d_bounds(l)
     &,') to ',  bounds(l+1),' (res=',d_bounds(l+1),')'

!       avg_res = average resolution of "t" cells within region
!       chg_res = change in resolution across the region
!       wid     = width of region
!       tol     = tolerance for fitting "t" cells within region width

!       provide for stretching last region if needed

        if (l .eq. n_bounds-1) then
          avg_res = p5*(d_bounds(l) + stretch*d_bounds(l+1))
          chg_res = (stretch*d_bounds(l+1) - d_bounds(l))
        else
          avg_res = p5*(d_bounds(l) + d_bounds(l+1))
          chg_res = (d_bounds(l+1) - d_bounds(l))
        endif

        tol = 1.e-5
        wid = abs(bounds(l+1) - bounds(l))
        an  = wid/avg_res
        m   = nint(an)

#if defined centered_t

!       Calculate resolution of "T" cells: "deltat". Note that
!       "T" grid points are centered in these cells (as in MOM 1)
!       n = number of "t" cells fitting within the region boundaries

        sum = 0.0
        n   = 0
        do i = 1,100000
          del = avg_res - p5*chg_res*cos((pi/m)*(i-0.5))
          if (sum + del .le. wid*(1.0 + tol)) then
            sum = sum + del
            if (num+i-1 .gt. maxlen) then
              write (stdout,*) "=>Error: maxlen exceeded in gcell. "
     &,                        " ...increase size of maxlen"
              stop
            endif
            deltat(num+i-1) = del
            n = n + 1
          else
            go to 100
          endif
        enddo
#else

!       Calculate resolution of "u" cells: "deltau"
!       "u" grid points will be centered in these cells
!       n = number of "t" cells fitting within the region boundaries
!       note: "sum" initially discounts half of the "u" cells widths
!       at the boundaries

        sum = 0.5*d_bounds(l) - 0.5*d_bounds(l+1)
        n   = 0
        do i = 1,100000
          del = avg_res - p5*chg_res*cos((pi/m)*i)
          if (sum + del .le. wid*(1.0 + tol)) then
            sum = sum + del
            if (num+i-1 .gt. maxlen) then
              write (stdout,*) "=>Error: maxlen exceeded in gcell. "
     &,                        " ...increase size of maxlen"
              stop
            endif
            deltau(num+i-1) = del
            n = n + 1
          else
            go to 100
          endif
        enddo
#endif
100     continue
        if (l .eq. n_bounds-1 .and. stretch .ne. 1.0) then
          write (stdout,'(a,i3,a,f5.2)')
     & '    constructed ',n,' cells with a stretch factor of ', stretch
          write (stdout,'(/2(a,g14.7),/2(a,g14.7),/a,a/)')
     &    'Note: you specified the ocean bottom at ',bounds(l+1)
     &,   ' cm with a bottom cell thickness of ',d_bounds(l+1)
     &,   '      The stretch factor puts the bottom at ',bounds(l)+sum
#if defined centered_t
     &,   ' cm with a bottom cell thickness of ',deltat(num+n-1)
#else
     &,   ' cm with a bottom cell thickness of '
     &,   0.5*(deltau(num+n-1) + deltau(num+n-2))
#endif
     &,   '      Adjust "stretch_z" in subroutine "gcoord" to get'
     &,   ' closer to the desired specifications if needed.'

        else
          write (stdout,'(a,g14.7,a)')
     &   '    constructed ',an,' grid cells for this region'
          if (abs(an-n) .gt. 0.01) then
            write (stdout, '(/,a,i2,/,a,g14.7/,a,g14.7,a//a/a)')
     & '==>Error: non integral number of cells in region #',l
     &,'          average resolution within region =',avg_res
     &,'          this implies ',an,' grid cells'
     &,'          Change grid specifications within USER INPUT section'
     &,'          Here is some help...'
            d_new = (2.0*wid)/(n-1) - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n-1,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            d_new = (2.0*wid)/n - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            d_new = (2.0*wid)/(n+1) - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n+1,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            stop '=>gcell'
          endif
        endif
        num = num + n
      enddo

!     adjust "num" to reflect the total number of cells contained in
!     all regions

      num = num - 1

      do i=1,num
#if defined centered_t

!       build resolution for "u" cells: "deltau". Note that
!       variable resolution (stretched grid) implies "u" points are
!       off centre as in MOM 1

        if (i .eq. num) then
          deltau(i) = d_bounds(n_bounds)
        else
          deltau(i) = p5*(deltat(i+1) + deltat(i))
        endif
#else

!       build resolution for "T" cells: "deltat". Note that
!       variable resolution (stretched grid) implies "T" points are
!       off centre

        if (i .eq. 1) then
          deltat(i) = p5*(d_bounds(1) + deltau(i))
        else
          deltat(i) = p5*(deltau(i) + deltau(i-1))
        endif
#endif
      enddo

!     add boundary points if needed

      if (nbpts .ne. 0) then
        do i=num,1,-1
          deltat(i+1) = deltat(i)
          deltau(i+1) = deltau(i)
        enddo
        deltat(1)     = deltat(2)
        deltau(1)     = d_bounds(1)
        deltat(num+2) = deltat(num+1)
        deltau(num+2) = deltau(num+1)
        num           = num + 2
      endif
      return
      end

#if !defined driver_only
      subroutine grids

!=======================================================================
!     set up a staggered "B" grid for MOM and compute grid related
!     variables

!     based on code by: R. C. Pacanowski
!=======================================================================

#include "param.h"
#include "coord.h"
#include "grdvar.h"
#if !defined implicitvmix || defined isopycmix || defined isoneutralmix
# include "accel.h"
#endif
#include "scalar.h"
#if defined uvic_mom
# include "hmixc.h"
# include "vmixc.h"
#endif

      integer ib(10), ic(10)
      character(120) :: fname, new_file_name

!-----------------------------------------------------------------------
!     set some constants
!-----------------------------------------------------------------------

      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      degtcm = radius/radian

!-----------------------------------------------------------------------
!     calculate coordinates for "t" and "u" grid cells.
!-----------------------------------------------------------------------

      maxlen = max(imt,jmt,km)
      imt2 = imt
      jmt2 = jmt
      km2 = km
      call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)

!-----------------------------------------------------------------------
!     verify that the number of grid points match the number set in
!     the parameter statement in "size.h".
!-----------------------------------------------------------------------

      call size_check (imt2, jmt2, km2, 'sub grids', 'stop')

!-----------------------------------------------------------------------
!     convert grid resolution to cm
!-----------------------------------------------------------------------

      do jrow=1,jmt
        dyt(jrow) = dytdeg(jrow)*degtcm
        dyu(jrow) = dyudeg(jrow)*degtcm
      enddo

      do i=1,imt
        dxt(i) = dxtdeg(i)*degtcm
        dxu(i) = dxudeg(i)*degtcm
      enddo
#if defined cyclic
      dxt(1)   = dxt(imt-1)
      dxt(imt) = dxt(2)
      dxu(1)   = dxu(imt-1)
      dxu(imt) = dxu(2)
#endif

!-----------------------------------------------------------------------
!     compute all quantities derived from the grid spacings
!-----------------------------------------------------------------------

      do k=1,km
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
      enddo

      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)

      do k=1,km
        dzwr(k-1)    = c1/dzw(k-1)
        dzw2r(k-1)   = p5/dzw(k-1)
      enddo

      do k=1,km
#if defined implicitvmix || defined isopycmix || defined isoneutralmix
        dztur(k) = c1/(dzw(k-1)*dzt(k))
        dztlr(k) = c1/(dzw(k)*dzt(k))
#endif
#if defined tcvmix
        dzwur(k) = c1/(dzt(k)*dzw(k))
        if (k .lt. km) dzwlr(k) = c1/(dzt(k+1)*dzw(k))
#endif
        dztr(k)  = c1/dzt(k)
      enddo

      tiny = 1.e-20
      do jrow=1,jmt
        dytr(jrow)  = c1/dyt(jrow)
        dyt2r(jrow) = p5/dyt(jrow)
        dyt4r(jrow) = p25/dyt(jrow)
        dyur(jrow)  = c1/dyu(jrow)
        dyu2r(jrow) = p5/dyu(jrow)
        dyu4r(jrow) = p25/dyu(jrow)
        phi(jrow)   = yu(jrow)/radian
        phit(jrow)  = yt(jrow)/radian
        cst(jrow)   = cos(phit(jrow))
        csu(jrow)   = cos(phi (jrow))
        sine(jrow)  = sin(phi(jrow))
        if (cst(jrow)  .eq. 0.0) then
          print '(/a,e14.7,a,i4,/a)'
     &  ,' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
          cst(jrow)  = tiny
        endif
        if (csu(jrow)  .eq. 0.0) then
          print '(/a,e14.7,a,i4,/a)'
     &  ,' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
          csu(jrow)  = tiny
        endif
        cstr(jrow)     = c1/cst(jrow)
        csur(jrow)     = c1/csu(jrow)
        tng(jrow)      = sine(jrow)/csu(jrow)
        cstdytr(jrow)  = c1/(cst(jrow)*dyt(jrow))
        cstdyt2r(jrow) = cstdytr(jrow)*p5
        csudyur(jrow)  = c1/(csu(jrow)*dyu(jrow))
        csudyu2r(jrow) = p5/(csu(jrow)*dyu(jrow))
        cst_dytr(jrow) = cst(jrow)/dyt(jrow)
        csu_dyur(jrow) = csu(jrow)/dyu(jrow)

      enddo

      do i=1,imt
        dxtr(i)  = c1/dxt(i)
        dxt2r(i) = p5/dxt(i)
        dxt4r(i) = p25/dxt(i)
        dxur(i)  = c1/dxu(i)
        dxu2r(i) = p5/dxu(i)
        dxu4r(i) = p25/dxu(i)
      enddo

      do i=2,imtm1
        dxmetr(i) = c1/(dxt(i) + dxt(i+1))
      enddo

      do i=1,imt
        duw(i) = (xu(i) - xt(i))*degtcm
      enddo
      do i=1,imtm1
        due(i) = (xt(i+1) - xu(i))*degtcm
      enddo
# if defined cyclic
      due(imt) = due(2)
# else
      due(imt) = due(imtm1)
# endif

      do jrow=1,jmt
        dus(jrow) = (yu(jrow) - yt(jrow))*degtcm
      enddo

      do jrow=1,jmtm1
        dun(jrow) = (yt(jrow+1) - yu(jrow))*degtcm
      enddo
      dun(jmt) = dun(jmtm1)

# if defined uvic_mom
#  if !defined implicitvmix || defined isopycmix || defined isoneutralmix

!     for convection code, compute values needed to include effects
!     of tracer timestep acceleration on effective layer thicknesses

      do k=1,km
        dztxcl(k) = dzt(k)/dtxcel(k)
      enddo

      do k=1,kmm1
        dzwxcl(k) = c1/(dztxcl(k)+dztxcl(k+1))
      enddo
      dzwxcl(km) = c0

#  endif
# endif

# if defined quicker

!     coefficients for 3rd order advective scheme

      do i=2,imt-1
        if (i .lt. imt-1) then
          ip2 = i+2
        else
          ip2 = 3
        endif
        quick_x(i,1) = c2*dxt(i+1)/(dxt(i+1)+dxt(i))
        quick_x(i,2) = c2*dxt(i  )/(dxt(i+1)+dxt(i))

        curv_xp(i,1) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i-1)+2.0*dxt(i)+dxt(i+1))*(dxt(i)+dxt(i+1)))
        curv_xp(i,2) =-c2*dxt(i)*dxt(i+1)/((dxt(i)+dxt(i+1))
     &                                  *(dxt(i-1)+dxt(i)))
        curv_xp(i,3) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i-1)+2.0*dxt(i)+dxt(i+1))*(dxt(i-1)+dxt(i)))

        curv_xn(i,1) = c2*dxt(i)*dxt(i+1)/
     &             ((dxt(i)+2.0*dxt(i+1)+dxt(ip2))*(dxt(i+1)+dxt(ip2)))
        curv_xn(i,2) =-c2*dxt(i)*dxt(i+1)/((dxt(i+1)+dxt(ip2))
     &                                  *(dxt(i)+dxt(i+1)))
        curv_xn(i,3) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i)+2.0*dxt(i+1)+dxt(ip2))*(dxt(i)+dxt(i+1)))
      enddo
      call setbcx (quick_x, imt, 2)
      call setbcx (curv_xp, imt, 3)
      call setbcx (curv_xn, imt, 3)

      do j=1,jmt
        jp2 = min(j+2,jmt)
        jp1 = min(j+1,jmt)
        jm1 = max(j-1,1)
        quick_y(j,1) = c2*dyt(jp1)/(dyt(jp1)+dyt(j))
        quick_y(j,2) = c2*dyt(j  )/(dyt(jp1)+dyt(j))

        curv_yp(j,1) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(jm1)+2.0*dyt(j)+dyt(jp1))*(dyt(j)+dyt(jp1)))
        curv_yp(j,2) =-c2*dyt(j)*dyt(jp1)/((dyt(j)+dyt(jp1))
     &                                  *(dyt(jm1)+dyt(j)))
        curv_yp(j,3) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(jm1)+2.0*dyt(j)+dyt(jp1))*(dyt(jm1)+dyt(j)))

        curv_yn(j,1) = c2*dyt(j)*dyt(jp1)/
     &             ((dyt(j)+2.0*dyt(jp1)+dyt(jp2))*(dyt(jp1)+dyt(jp2)))
        curv_yn(j,2) =-c2*dyt(j)*dyt(jp1)/((dyt(jp1)+dyt(jp2))
     &                                  *(dyt(j)+dyt(jp1)))
        curv_yn(j,3) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(j)+2.0*dyt(jp1)+dyt(jp2))*(dyt(j)+dyt(jp1)))
      enddo

      do k=1,km
        kp2 = min(k+2,km)
        kp1 = min(k+1,km)
        km1 = max(k-1,1)
        quick_z(k,1) = c2*dzt(kp1)/(dzt(kp1)+dzt(k))
        quick_z(k,2) = c2*dzt(k  )/(dzt(kp1)+dzt(k))

        curv_zp(k,1) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(km1)+2.0*dzt(k)+dzt(kp1))*(dzt(k)+dzt(kp1)))
        curv_zp(k,2) =-c2*dzt(k)*dzt(kp1)/((dzt(k)+dzt(kp1))
     &                                  *(dzt(km1)+dzt(k)))
        curv_zp(k,3) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(km1)+2.0*dzt(k)+dzt(kp1))*(dzt(km1)+dzt(k)))

        curv_zn(k,1) = c2*dzt(k)*dzt(kp1)/
     &             ((dzt(k)+2.0*dzt(kp1)+dzt(kp2))*(dzt(kp1)+dzt(kp2)))
        curv_zn(k,2) =-c2*dzt(k)*dzt(kp1)/((dzt(kp1)+dzt(kp2))
     &                                  *(dzt(k)+dzt(kp1)))
        curv_zn(k,3) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(k)+2.0*dzt(kp1)+dzt(kp2))*(dzt(k)+dzt(kp1)))
      enddo

# endif
      fname = new_file_name ("grid.nc")
      ib(:) = 1
      ic(:) = imt
      ic(2) = jmt
      call openfile (fname, iou)
      call getvara ('lat_t', iou, imt*jmt, ib, ic, tlat, c1, c0)
      call getvara ('lon_t', iou, imt*jmt, ib, ic, tlon, c1, c0)
      call getvara ('lat_u', iou, imt*jmt, ib, ic, ulat, c1, c0)
      call getvara ('lon_u', iou, imt*jmt, ib, ic, ulon, c1, c0)
      call closefile (iou)

      return
      end
#endif
