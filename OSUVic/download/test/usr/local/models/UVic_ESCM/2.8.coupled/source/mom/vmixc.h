!====================== include file "vmixc.h" =========================

!         vertical mixing coefficients and related variables

!     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
!     kappa_m = constant vertical viscosity coefficient (cm**2/sec)

!     visc_cbu  = viscosity coeff at bottom of "u" cell
!     diff_cbt  = diffusion coeff at bottom of "T" cell
!     visc_cbu_limit = largest allowable "visc_cbu"
!     diff_cbt_limit = largest allowable "diff_cbt"
!     aidif = coefficient for implicit time differencing for
!             vertical diffusion. aidif=1 gives the fully implicit
!             case. aidif=0 gives the fully explicit case
!             note: not used unless "implicitvmix" or "isopycmix"
!                   is enabled

!     based on code by: R. C. Pacanowski
!=======================================================================

      real kappa_h,  kappa_m
      common /vmixr0/ visc_cbu_limit, diff_cbt_limit, aidif
      common /vmixr0/ kappa_h, kappa_m
#if !defined ppvmix
# if defined coarse_grained_parallelism
!DIR$ TASKCOMMON vmixr1
!CDIR TASKLOCAL (/vmixr1/)
!IBM* THREADLOCAL /vmixr1/
# endif
      common /vmixr1/ visc_cbu(imt,km,jsmw:jemw)
      common /vmixr1/ diff_cbt(imt,km,jsmw:jemw)
#endif

#if defined bryan_lewis_vertical
      common /vmixr0/ Ahv(km)
#endif

#if defined ppvmix

!     variables for pacanowski-philander vertical diffusion

!     rhom1z = rho(k)-rho(k+1)
!     uzsq   = (u(k)-u(k+1))**2 + (v(k)-v(k+1))**2
!     visc_cbt  = viscosity coeff at bottom of "T" cell
!     fricmx = max vertical mixing coefficient
!     wndmix = min vertical mixing in level 1 to simulate wind mixing
!     diff_cbt_back = background "diff_cbt"
!     visc_cbu_back = background "visc_cbu"

# if defined coarse_grained_parallelism
!DIR$ TASKCOMMON vmixr2
!CDIR TASKLOCAL (/vmixr2/)
!IBM* THREADLOCAL /vmixr2/
# endif
      common /vmixr0/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
      common /vmixr2/ rhom1z(imt,km,jmw), uzsq(imt,km,jmw)
      common /vmixr2/ diff_cbt(imt,km,jsmw:jmw)
      common /vmixr2/ visc_cbt(imt,km,jsmw:jmw)
      common /vmixr2/ visc_cbu(imt,km,jsmw:jemw)
#endif
