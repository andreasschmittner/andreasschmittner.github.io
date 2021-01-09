!======================= include file "size.h" =========================

!-----------------------------------------------------------------------
!     USER INPUT:
!-----------------------------------------------------------------------

!     imt    = number of grid points in the longitudinal direction
!              (calculated points are from 2 through imt-1. End points
!               are boundaries)
!     jmt    = number of grid points (latitude rows) in the latitudinal
!              direction (calculated points are from 2 through jmt-1.
!              End points are boundaries)
!     km     = number of grid points in the vertical direction
!              (calculated points are from 1 through km)
!     nt     = number of tracers (temperature, salinity, ...)
!     nsrc   = number of tracer with sources
!     kpzd   = depth for limited npzd model
!     jmz    = size for "unrotated" zonal averages
!     jmzm1  = jmz minus one
!     mnisle = maximum number of islands (unconnected land masses)
!     maxipp = maximum number of all island perimeter points
!-----------------------------------------------------------------------
# include "derived_options.h"

      integer imt, jmt, km, nt, nsrc, kpzd, nat, jmz, jmzm1, mnisle
      integer maxipp, nprocessors, jmw, jsmw, jemw, jextra

      parameter (imt=  102, jmt=  102, km= 19)
      parameter (nt=2
#if defined uvic_carbon
     $             +1
# if defined uvic_carbon_14
     $             +1
# endif
#endif
#if defined uvic_cfc11
     $             +1
#endif
#if defined uvic_cfc12
     $             +1
#endif
#if defined uvic_alk
     $             +1
#endif
#if defined uvic_o2
     $             +1
#endif
#if defined uvic_npzd
     $             +4
# if defined uvic_nitrogen
     $             +2
# endif
#endif
     $               )
      parameter (nsrc=0
#if defined uvic_carbon
     $               +1
# if defined uvic_carbon_14
     $               +1
# endif
#endif
#if defined uvic_alk
     $               +1
#endif
#if defined uvic_o2
     $               +1
#endif
#if defined uvic_npzd
     $               +4
# if defined uvic_nitrogen
     $               +2
# endif
#endif
     $                 )
      parameter (kpzd=km)

      parameter (nat=2
#if defined uvic_carbon && defined uvic_carbon_co2_2d
     $              +1
#endif
     $, jmz=jmt, jmzm1=jmz-1)
      parameter (mnisle=50, maxipp=5000)

#if defined coarse_grained_parallelism && defined ramdrive
      parameter (nprocessors=4)
#endif

#if defined obctest
# if defined obc_south || defined obc_north
      parameter (imt= 21, jmt= 24, km= 6 )
# else
      parameter (imt= 21, jmt= 41, km= 6 )
# endif
#endif
#if defined obctest2
# if defined obc_south
      parameter (imt= 21, jmt= 11, km= 6 )
# else
      parameter (imt= 21, jmt= 21, km= 6 )
# endif
#endif

#if defined obc_south || defined obc_north || defined obc_west || defined obc_east && !defined obc
# define obc
#endif
#if !defined coarse_grained_parallelism
# if !defined uvic_min_window
      parameter (jmw=jmt)
# else

!     for UNI-TASKING: "jmw" is set to the minimum for each option class
!     "jmw" may be increased up to "jmt"

#  if defined fourth_order_window
#   if defined pressure_gradient_average
#    if defined biharmonic  || defined fourth_order_tracer_advection || defined fct || defined quicker
      parameter (jmw=5)
#    else
      parameter (jmw=4)
#    endif
#   else
      parameter (jmw=4)
#   endif
#  else
      parameter (jmw=3)
#  endif
# endif
#else

!     for MICROTASKING: window size for coarse grained parallelism only
!     (these settings should not be changed)

# if defined fourth_order_window
#  if defined pressure_gradient_average
#   if defined biharmonic  || defined fourth_order_tracer_advection || defined fct || defined quicker
      parameter (jmw=7)
#   else
      parameter (jmw=5)
#   endif
#  else
      parameter (jmw=5)
#  endif
# else
      parameter (jmw=3)
# endif
#endif

!-----------------------------------------------------------------------
!     set first and last calculated row within the MW. other rows
!     are used as buffers
!-----------------------------------------------------------------------

!     jsmw   = 1st calculated row within the MW
!     jemw   = last calculated row within the MW

      parameter (jsmw=2, jemw=jmw-1)

!     jextra = extra buffer rows needed for coarse_grained_parallelism.
!              "jextra" rows are added to the top and bottom of the MW
!              so the MW size increases by 2*jextra

#if defined coarse_grained_parallelism
      parameter (jextra = 0
# if defined fourth_order_window
#  if defined pressure_gradient_average
     &                    + 1
#  endif
#  if defined biharmonic || defined fourth_order_tracer_advection || defined fct || defined quicker
     &                    + 1
#  endif
# endif
     &                    )
#else
          parameter (jextra = 0)
#endif
#if !defined fct && !defined fourth_order_tracer_advection && !defined quicker && !defined second_order_tracer_advection
# define second_order_tracer_advection
#endif
#if defined uvic_mtlm
! POINTS = Maximum number of points in grid.
! STEPSM = Maximum number of timesteps in a day.
! klmax = maximum ocean depth levels over which the land model can exist

      integer POINTS, STEPSM, klmax
      parameter (POINTS=4600, STEPSM=24, klmax=0)

! NNVG  = Number of non-vegetation surface types.
! NPFT  = Number of plant functional types.
! NTYPE = Number of surface types.
! SOIL  = Index of the surface type 'Soil'
! Land surface types :
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!     6 - Soil

      integer NNVG, NPFT, NTYPE, SOIL
      parameter (NNVG=4, NPFT=5, NTYPE=6, SOIL=6)
#endif
