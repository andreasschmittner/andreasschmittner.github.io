!====================== include file "ptraj.h" =========================

!     inputs:

!     nptraj  = number of particle trajectories. particles are
!               assumed to be of zero mass, neutrally buoyant and move
!               with the local three dimensional flow during each time
!               step.

!               initially, particles are randomly distributed within
!               the volume given by:

!     ptslon  = starting longitude for initial particle distribution
!     ptelon  = ending longitude for initial particle distribution
!     ptslat  = starting latitude for initial particle distribution
!     ptelat  = ending latitude for initial particle distribution
!     ptsdpt  = starting depth for initial particle distribution
!     ptedpt  = ending depth for initial particle distribution

!     outputs:

!     pxyz    = particle coordinates. index (1,2,3) is for particle
!               (longitude, latitude, depth).

!     pijk    = the particle is bounded by the volume with verticies
!               given by the eight nearest surrounding model grid points
!               on the "xu","yu", and "z" grids. index (1,2,3) locates
!               the (longitude, latitude, depth) index of the deepest
!               northeast corner of this bounding volume.

!     em      = matrix of deformation rates for calculation of
!               lyapunov exponents

!     ptdone  = boolean for limiting multiple passes on trajectories
!               within one time step

!     initpt  = boolean for initiailizing particle positions
!               (t,f) = (initialize, do not initialize)

!     based on code by: R. C. Pacanowski

      parameter (nptraj = 200)
      integer pijk
      logical ptdone, initpt

      common /cptrji/ pijk(3,nptraj)
      common /cptraj/ pxyz(3,nptraj)
      common /cptraj/ ptslon, ptelon, ptslat, ptelat, ptsdpt, ptedpt
#if defined lyapunov
      common /cptraj/ em(2,2,nptraj)
#endif
      common /cptra2/ ptdone(nptraj), initpt
