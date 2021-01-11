#include "derived_options.h"
      integer ncat           ! ice categories (not counting open water)
      integer ntilay         ! equal to the sum of the elements of nilay
#if defined uvic_ice_cpts1
      parameter(ncat=1)
      parameter(ntilay=4)
#elif defined uvic_ice_cpts3
      parameter(ncat=3)
      parameter(ntilay=14)
#elif defined uvic_ice_cpts5
      parameter(ncat=5)
      parameter(ntilay=26)
#elif defined uvic_ice_cpts10
      parameter(ncat=10)
      parameter(ntilay=56)
#endif
      integer itme           ! number of time levels kept
      parameter(itme = 2)

!     upper limit on the number of layers in a single category
!     fine if bigger than necessary
      integer nmax
      parameter (nmax=10)

!     these common blocks must be initialized at setup, see setembm.F

      real dtau  ! time step (2*dtatm=>leapfrog, dtatm=>forward)
      common/tstep/dtau

      integer idx            ! index for time step
      common/counters/idx

      real M_def(ncat,ncat)  !     area frac from cat i that goes into j
      real N_def(ncat,ncat)  !   volume frac from cat i that goes into j
      real HN_def(ncat,ncat) ! H*volume frac from cat i that goes into j
      common /defmat/ M_def,N_def,HN_def

      common /layinf/ nilay, layer1, layern, ncrel
      integer nilay(ncat)      ! the number of layers
      integer layer1(ncat)     ! position of the top layer
      integer layern(ncat)     ! position of the bottom layer
      integer ncrel(ncat,ncat) ! helps remap layers when ridging

      common /small/ asmall, hstar
      real hstar(0:ncat)  ! category limits (m)
      real asmall(0:ncat) ! minimum fract. area of allowed

      common /salt/ saltz(nmax+1,ncat),tmelz(nmax,ncat), salnew

      real saltz   ! salinity of each layer (ppm)
      real tmelz   ! melting temperature of each layer (C)
      real salnew  ! salinity of new ice

!     default units are mks so use these factors to
!     convert various parameters from mks to cgs
!     set them equal to 1. if no conversion is desired
      real deci  ! factor of 10.
      real centi ! factor of 100.
      real kilo  ! factor of 1000.

      parameter ( deci=10., centi=100.0, kilo=1000.)
!      parameter ( deci = 1.,  centi = 1.,    kilo = 1.)

!     parameter for small things
      real TINY
      parameter( TINY=1.0e-10 )

!     if open water area > GSTAR then no mechanical redist.
      real GSTAR                   ! portion of itd that may partic.
      parameter ( GSTAR = 0.15 )

!     max ridged ice thickness is 2*sqrt(cK*Hi), Hi is orig thickness
      real cK  ! used to determine max. ridged thickness (m)
      parameter ( cK = 1.e2*centi)

!     if using Rothrock75 energetics argument for pressure
      real cpe ! coef. for potential energy term(N/m**3)
      parameter ( cpe = 450./deci)

!     minimum allowable open water/ice fraction
!     best if small as possible, can cause round off error if too small
      real a0small
      real aismall
      parameter ( a0small=0.0001, aismall=0.0005 )

!     minimum snow thickness allowed
      real hsstar
      parameter ( hsstar=0.00001*centi )
