!======================= include file "hmixc.h" ========================

!                    horizontal mixing coefficients

!     visc_cnu = viscosity coeff for northern face of "u" cell
!     visc_ceu = viscosity coeff for eastern face of "u" cell
!     diff_cnt = diffusion coeff for northern face of "T" cell
!     diff_cet = diffusion coeff for eastern face of "T" cell

!     am     = constant lateral viscosity coeff for momentum
!     ah     = constant lateral diffusion coeff for tracers
!     am3    = viscosity coeff for metric term on "u" cell
!     am4    = another viscosity coeff for metric term on "u" cell
!     ambi   = constant lateral biharmonic viscosity coeff for momentum
!     ahbi   = constant lateral biharmonic diffusion coeff for tracers

!     based on code by: R. C. Pacanowski
!=======================================================================

#if defined consthmix
      common /diffus0/ am, ambi, am3(jmt), am4(jmt,2)
      common /diffus0/ ah, ahbi
      common /diffus0/ visc_ceu, visc_cnu
      common /diffus0/ amc_north(jmt), amc_south(jmt)

# if defined bryan_lewis_horizontal

!     bryan_lewis mixing case

      common /diffus0/ Ahh(km)
      common /diffus0/ diff_cnt(km), diff_cet(km)
      common /diffus0/ ahc_north(jmt,km), ahc_south(jmt,km)
# else
      common /diffus0/ diff_cnt, diff_cet
      common /diffus0/ ahc_north(jmt), ahc_south(jmt)
# endif
#endif

#if defined smagnlmix

!     non-linear horizontal viscosity after Smagorinsky 1963,
!     as described in Rosati & Miyakoda (jpo,vol 18,#11,1988)
!     see Smagorinsky 1963, Mon Wea Rev, 91, 99-164.
!     Also see Deardorff 1973 J. Fluid Eng. Sep., 429-438.

!     strain = tension(1) and shearing(2) rates of strain
!     smag_metric  = metric term
!     diff_c_back = background diffusion coeff for "t" cell (cm**2/sec)

# if defined coarse_grained_parallelism
!DIR$ TASKCOMMON diffus
!CDIR TASKLOCAL (/diffus/)
!IBM* THREADLOCAL /diffus/
# endif
      common /diffus/ strain(imt,km,1:jemw,2)
      common /diffus/ am_lambda(imt,km,1:jemw), am_phi(imt,km,1:jemw)
      common /diffus/ smag_metric(imt,km,jsmw:jemw)
      common /diffus/ diff_c_back
      common /diffus/ visc_ceu(imt,km,jsmw:jemw)
      common /diffus/ visc_cnu(imt,km,1:jemw)
      common /diffus/ diff_cet(imt,km,jsmw:jemw)
      common /diffus/ diff_cnt(imt,km,1:jemw)
#endif

#if defined held_larichev

!     variables for held_larichev diffusion

!     hl_u     = Held/Larichev diffusion coefficient (cm**2/sec)
!                defined over "u" cells
!     hl_n     = avg_x(hl_u) defined at northern face of "T" cells
!     hl_e     = avg_y(hl_u) defined at eastern face of "T" cells
!     hl_b     = avg_x(avg_y(hl_u))
!     hl_depth = depth of integration (cm)
!     hl_back  = background mixing coeff (cm**2/sec)
!     hl_max   = max allowed mixing coeff (cm**2/sec)
!     droz     = vertical difference of rho

# if defined coarse_grained_parallelism
!DIR$ TASKCOMMON diffus1
!CDIR TASKLOCAL (/diffus1/)
!IBM* THREADLOCAL /diffus1/
# endif
      common /diffus1/ hl_depth, hl_back, hl_max
      common /diffus1/ hl_u(imt,1:jemw)
      common /diffus1/ hl_n(imt,1:jemw), hl_e(imt,jsmw:jemw)
      common /diffus1/ hl_b(imt,jsmw:jemw)
      common /diffus1/ droz(imt,km,jmw), rich_inv(imt,km,1:jemw)
#endif
