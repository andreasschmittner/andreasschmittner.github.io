!====================== include file "termbal.h" =======================

!     term balances are instantaneous breakdowns of all terms in the
!     momentum & tracer equations. They are averaged over ocean volumes
!     defined by horizontal and vertical regional masks:

!     termbt   = term balance components for time rate of change of
!                tracers within a  volume. the total time rate of change
!                is broken down into components as follows:
!                the form is d( )/dt = terms (2) ... (10) where each
!                term has the units of "tracer units/sec" using
!                schematic terms for illustration.
!                (1)  = total time rate of change for the tracer
!                (2)  = change due to zonal nonlinear term: (UT)x
!                (3)  = change due to meridional nonlinear term: (VT)y
!                (4)  = change due to vertical nonlinear term: (WT)z
!                (5)  = change due to zonal diffusion: Ah*Txx
!                (6)  = change due to meridional diffusion: Ah*Tyy
!                (7)  = change due to vertical diffusion:  kappa_h*Tzz
!                (8)  = change due to source term
!                (9)  = change due to explicit convection
!                (10) = change due to filtering
!     the nonlinear terms can be broken into two parts: advection and a
!     continuity part: The physically meaningful part is advection.
!     eg: Zonal advection of tracer "A" is -U(A)x = A(U)x - (UA)x
!                (11) = zonal advection U(Ax)
!                (12) = meridional advection V(Ay)
!                (13) = vertical advection W(Az)
!                (14) = change of tracer variance (tracer**2 units)
!                (15) = average tracer within volume (tracer units)
!     terr     = error term = (1) - sum (2) ... (10)
!     asst     = average sea surface tracer for regional surface areas
!     stflx    = average surface tracer flux for regional surface areas
!                tracer (#1,#2) units = (cal/cm**2/sec, gm/cm**2/sec)

!     termbm   = term balance components for time rate of change of
!                momentum within a volume. the total time rate of change
!                is broken down into components as follows:
!                the form is d( )/dt = terms (2) ... (13) where each
!                term has the units of "cm/sec**2" and "Q" is the
!                momentum component {zonal or meridional} using
!                schematic terms for illustration.
!                (1)  = total time rate of change for the momentum
!                (2)  = change due to the pressure gradient: grad_p
!                       without the surface pressure gradients
!                       (i.e., for computing the internal modes)
!                (3)  = change due to zonal nonlinear term: (UQ)x
!                (4)  = change due to meridional nonlinear term: (VQ)y
!                (5)  = change due to vertical nonlinear term: (wQ)z
!                (6)  = change due to zonal viscosity: Am*Qxx
!                (7)  = change due to meridional viscosity: Am*Qyy
!                (8)  = change due to vertical viscosity: kappa_m*Qzz
!                (9)  = change due to metric diffusion terms
!                (10) = change due to coriolis terms: fQ
!                (11) = change due to source terms
!                (12) = change due to surface pressure gradient
!                       this is obtained after solving the external mode
!                       in the stream function technique. It is solved
!                       directly from the elliptic equation for the
!                       prognostic surface pressure technique
!                (13) = change due to metric advection
!     the nonlinear terms can be broken into two parts: advection and a
!     continuity part: The physically meaningful part is advection.
!     eg: Zonal advection of vel component "Q" is -U(Q)x = Q(U)x - (UQ)x
!                (14) = zonal advection U(Qx)
!                (15) = meridional advection V(Qy)
!                (16) = vertical advection W(Qz)
!                (17) = average velocity component
!     smflx    = average surface momentum flux for regional surf areas
!                in dynes/cm**2
!     avgw     = average vertical velocity (cm/sec)

!     ustf     = names & units for surface tracer fluxes

!     based on code by: R. C. Pacanowski

      parameter (ntterms=15, nuterms=17)

      character(15) :: ustf(nt,2)
      common /termb1/ ustf
      common /termb2/asst(nt,0:nhreg), avgw(numreg)
      common /termb2/termbt(0:km,ntterms,nt,0:numreg)
      common /termb2/termbm(0:km,nuterms,2,numreg)
      common /termb2/smflx(2,0:nhreg)
      common /termb2/stflx(nt,0:nhreg), terr(nt)
