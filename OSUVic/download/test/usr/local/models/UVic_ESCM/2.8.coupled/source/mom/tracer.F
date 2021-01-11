      subroutine tracer (joff, js, je, is, ie)

!=======================================================================
!     compute tracers at "tau+1" for rows js through je in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

!     based on code by: R. C. Pacanowski
!=======================================================================

#include "param.h"
      parameter (istrt=2, iend=imt-1)
#include "accel.h"
#if defined obc_north
# include "cobc.h"
#endif
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#if defined isopycmix || defined isoneutralmix
# include "isopyc.h"
#endif
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "timeavgs.h"
#include "tmngr.h"
#include "vmixc.h"
      dimension twodt(km)
#if defined llnl_plume
      dimension  work(imt,km,jsmw:jmw)
#endif
#include "fdift.h"
#if defined uvic_save_convection || defined uvic_carbon_14
# include "diaga.h"
#endif
#if defined uvic_ice
# include "ice.h"
#endif
#if defined uvic_npzd
# if defined uvic_embm
#  include "atm.h"
# endif
# include "npzd.h"

      integer index
      real snpzd(ntnpzd), tnpzd(ntnpzd)
      real rctheta, declin, gl, impo, expo, npp
      real remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real phin, dz, prca, dprca, nud, bct, tap, fo2, so2
# if defined uvic_nitrogen
      real npp_D, graz_D, morp_D, no3flag, deni, nfix
# endif
#endif
#if defined uvic_carbon_14
      real fy, fyz
#endif
#if defined uvic_npzd || defined uvic_carbon_14
      real src(imt,km,jsmw:jmw,nsrc)
      src(:,:,:,:) = 0.
#endif

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

#if defined fourth_order_window || defined isopycmix || defined isoneutralmix
      limit = min(je+1+joff,jmt) - joff
      do j=js,limit
#else
      do j=js,je
#endif
        jrow = j + joff
        do i=istrt-1,iend
          cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
          cstdxur(i,j)    = cstr(jrow)*dxur(i)
#if defined consthmix
          ah_cstdxur(i,j) = diff_cet*cstr(jrow)*dxur(i)
#endif
        enddo
      enddo
#if defined llnl_plume

!-----------------------------------------------------------------------
!     find depth of penetration for the llnl simple plume model
!-----------------------------------------------------------------------

!     set parameters
      zmax = 160.0e2 ! use constant depth penetration (if cont=0)
!      cont = 0.0    ! use rho contrast for penetration (if cont>0)
      cont = 0.4e-3  ! 0.4 kg/m3 = 0.4e-3g/cm3

      do j=js,je
        jrow = j + joff
        do i=is,ie
          subz(i,jrow) = zmax
        enddo
        if (cont .gt. c0) then
          call state_ref (t(1,1,1,1,tau), t(1,1,1,2,tau)
     &,                   work(1,1,jsmw), j, j, 1, imt, 1)
          do i=is,ie
            kmax = kmt(i,jrow)
            subz(i,jrow) = zw(1)
            do k=2,kmax
              drho = work(i,k,j) - work(i,1,j)
              drhom1 = work(i,k-1,j) - work(i,1,j)
              if (drho .ge. cont .and. drhom1 .lt. cont)
     &          subz(i,jrow) = zt(k-1) + (zt(k) - zt(k-1))
     &                      *(cont - drhom1)/(drho - drhom1)
            enddo
            if (drho .le. cont) subz(i,jrow) = zw(kmax)
          enddo
        endif

!-----------------------------------------------------------------------
!     calculate "3d stf" multiplier for the llnl simple plume model
!-----------------------------------------------------------------------

        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            subz(i,jrow) = max(min(subz(i,jrow),zw(kmt(i,jrow))),zw(1))
            work(i,1,j) = c1/subz(i,jrow)
            do k=2,km
              wt =  min(max(c0,(subz(i,jrow) - zw(k-1))),dzt(k))*dztr(k)
              work(i,k,j) = wt*work(i,1,j)
            enddo
          else
            do k=1,km
              work(i,k,j) = c0
            enddo
          endif
        enddo
      enddo
#endif
#if defined uvic_npzd

!-----------------------------------------------------------------------
!     calculation of biological interactions
!-----------------------------------------------------------------------

      declin = sin((mod(relyr,1.) - 0.22)*2.*pi)*0.4   ! declination

      do k=1,km
        twodt(k) = c2dtts*dtxcel(k)
        nbio(k) = twodt(k)/dtnpzd
        dtbio(k) = twodt(k)/nbio(k)
        rdtts(k) = 1./twodt(k)
        rnbio(k) = 1./nbio(k)
      enddo
      tap = 2.*alpha*par

      do j=js,je
        jrow = j + joff

        do i=is,ie

          if (kmt(i,jrow) .gt. 0) then

!           calculate day fraction and incoming solar
!           angle of incidence = lat - declin, refraction index = 1.33
            rctheta = max(-1.5, min(1.5, tlat(i,jrow)/radian - declin))
            rctheta = kw/sqrt(1. - (1. - cos(rctheta)**2.)/1.33**2.)
            dayfrac = min( 1., -tan(tlat(i,jrow)/radian)*tan(declin))
            dayfrac = max(1e-12, acos(max(-1., dayfrac))/pi)
# if defined uvic_embm
            swr = dnswr(i,jrow)*1e-3   ! convert to W/m^2
#  if defined uvic_ice
     &          *(1.0 + aice(i,jrow,2) ! attenuation by sea ice & snow
     &          *(exp(-ki*(hice(i,jrow,2) + hsno(i,jrow,2))) -1.))
#  endif
# else
            swr = 200.
# endif
            expo = 0.0
            impo = 0.0
            phin = 0.0 ! integrated phytoplankton
            prca = 0.0 ! integrated production of calcite

            kmax = min(kmt(i,jrow), kpzd)
            do k=1,kmax

!-----------------------------------------------------------------------
!             limit tracers to positive values
!-----------------------------------------------------------------------

              tnpzd(1) = max(t(i,k,j,ipo4,taum1), trcmin)
              tnpzd(2) = max(t(i,k,j,iphyt,taum1), trcmin)
              tnpzd(3) = max(t(i,k,j,izoop,taum1), trcmin)
              tnpzd(4) = max(t(i,k,j,idetr,taum1), trcmin)
# if defined uvic_nitrogen
              tnpzd(5) = max(t(i,k,j,ino3,taum1), trcmin)
              tnpzd(6) = max(t(i,k,j,idiaz,taum1), trcmin)
# endif

              swr = swr*exp(-kc*phin)
              phin = phin + tnpzd(2)*dzt(k)
              gl = tap*swr*exp(ztt(k)*rctheta)
              impo = expo*dztr(k)
              bct = bbio**(cbio*t(i,k,j,itemp,taum1))
# if defined uvic_o2
!             decrease remineralisation rate in oxygen minimum zone
              nud = nud0*(0.65+0.35*tanh(t(i,k,j,io2,taum1)*1000.-6.))
# else
              nud = nud0
# endif

!-----------------------------------------------------------------------
!             call the npzd model
!-----------------------------------------------------------------------

              call npzd_src (tnpzd, nbio(k), dtbio(k), gl, bct, impo
     &,                      dzt(k), dayfrac, wd(k), rkwz(k), nud
     &,                      snpzd, expo, graz, morp, morz
# if defined uvic_npzd_out
     &,                      npp, morpt, remi, excr
#  if defined uvic_nitrogen
     &,                      npp_D, graz_D, morp_D, nfix
#  endif
# endif
     &                      )

              snpzd(1:4) = snpzd(1:4)*rdtts(k)
# if defined uvic_nitrogen
              snpzd(5:6) = snpzd(5:6)*rdtts(k)
# endif
              expo = expo*rnbio(k)
# if defined uvic_npzd_out
              rexpo(i,k,j) = expo
              rgraz(i,k,j) = graz*rnbio(k)
              rmorp(i,k,j) = morp*rnbio(k)
              rmorz(i,k,j) = morz*rnbio(k)
              rnpp(i,k,j) = npp*rnbio(k)
              rmorpt(i,k,j) = morpt*rnbio(k)
              rremi(i,k,j) = remi*rnbio(k)
              rexcr(i,k,j) = excr*rnbio(k)
#  if defined uvic_nitrogen
              rnpp_D(i,k,j) = npp_D*rnbio(k)
              rgraz_D(i,k,j) = graz_D*rnbio(k)
              rmorp_D(i,k,j) = morp_D*rnbio(k)
              rnfix(i,k,j) = nfix*rnbio(k)
#  endif
# endif
!-----------------------------------------------------------------------
!             calculate detritus at the bottom
!-----------------------------------------------------------------------

              if (k .eq. kmt(i,jrow)) then
                snpzd(4) = snpzd(4) + expo
                expo = 0.
# if defined uvic_npzd_out
                rexpo(i,k,j) = 0.
# endif
              elseif (k .eq. kpzd .and. expo .ne. 0.) then

                expo = expo*dzt(kpzd)

                do m=kpzd+1,kmt(i,jrow)

                  expo = expo*dztr(m)
!                 if below the npzd model move all biology to export
                  src(i,m,j,isphyt) = -rdtts(m)*t(i,m,j,iphyt,taum1)
                  t(i,m,j,iphyt,taup1) = 0.
                  src(i,m,j,iszoop) = -rdtts(m)*t(i,m,j,izoop,taum1)
                  t(i,m,j,izoop,taup1) = 0.
                  src(i,m,j,isdetr) = -rdtts(m)*t(i,m,j,idetr,taum1)
                  t(i,m,j,idetr,taup1) = 0.
# if defined uvic_nitrogen
                  src(i,m,j,isdiaz) = -rdtts(m)*t(i,m,j,idiaz,taum1)
                  t(i,m,j,idiaz,taup1) = 0.
                  expo = expo - src(i,m,j,isphyt) - src(i,m,j,iszoop)
     &                 - src(i,m,j,isdetr) - src(i,m,j,isdiaz)
# else
                  expo = expo - src(i,m,j,isphyt) - src(i,m,j,iszoop)
     &                 - src(i,m,j,isdetr)
# endif
                  if (m .eq. kmt(i,jrow)) then
!                   remineralize all detritus if at the bottom
                    src(i,m,j,ispo4) = redptn*expo
# if defined uvic_npzd_out
                    rremi(i,m,j) = expo
                    rexpo(i,m,j) = 0.
# endif
                  else
!                   calculate the proportion of detritus remineralized
                    nud = nud0*bbio**(cbio*t(i,m,j,itemp,taum1))
#  if defined uvic_o2
!                   decrease remineralisation in oxygen minimum zone
     &                  *(0.65+0.35*tanh(t(i,m,j,io2,taum1)*1000.-6.))
#  endif
                    remi = expo*nud/(nud + wd(m))
                    src(i,m,j,ispo4) = redptn*remi
# if defined uvic_npzd_out
                    rremi(i,m,j) = remi
                    rexpo(i,m,j) = expo - remi
# endif
!                   calculate the detritus exported to the next level
                    expo = (expo - remi)*dzt(m)
                  endif
# if defined uvic_nitrogen
                  src(i,m,j,isno3) = src(i,m,j,ispo4)*redntp
# endif
# if defined uvic_carbon
                  src(i,m,j,isdic) = src(i,m,j,ispo4)*redctp
# endif
# if defined uvic_alk
                  src(i,m,j,isalk) = -src(i,m,j,ispo4)*redntp*1.e-3
# endif
                enddo
              endif

!-----------------------------------------------------------------------
!             set source terms
!-----------------------------------------------------------------------

              src(i,k,j,ispo4) = snpzd(1)
              src(i,k,j,isphyt) = snpzd(2)
              src(i,k,j,iszoop) = snpzd(3)
              src(i,k,j,isdetr) = snpzd(4)
# if defined uvic_nitrogen
              src(i,k,j,isno3) = snpzd(5)
              src(i,k,j,isdiaz) = snpzd(6)
# endif
# if !defined uvic_npzd_no_calcite
!             production of calcite
              dprca = (morp+morz+graz*(1.-gamma1))*capr*redctn*rnbio(k)
              prca = prca + dprca*dzt(k)
#  if defined uvic_carbon
              src(i,k,j,isdic) = (snpzd(1)*redctp - dprca)
#  endif
#  if defined uvic_alk
              src(i,k,j,isalk) = (-snpzd(1)*redntp*1.e-3 - 2.*dprca)
#  endif
# endif
# if defined time_averages && defined uvic_npzd_out

!-----------------------------------------------------------------------
!             accumulate time averages
!-----------------------------------------------------------------------

              if (timavgperts .and. .not. euler2) then
                ta_rnpp(i,k,jrow) = ta_rnpp(i,k,jrow) + rnpp(i,k,j)
                ta_rgraz(i,k,jrow) = ta_rgraz(i,k,jrow) + rgraz(i,k,j)
                ta_rmorp(i,k,jrow) = ta_rmorp(i,k,jrow) + rmorp(i,k,j)
                ta_rmorpt(i,k,jrow)= ta_rmorpt(i,k,jrow) + rmorpt(i,k,j)
                ta_rmorz(i,k,jrow) = ta_rmorz(i,k,jrow) + rmorz(i,k,j)
                ta_rexcr(i,k,jrow) = ta_rexcr(i,k,jrow) + rexcr(i,k,j)
#  if defined uvic_nitrogen
                ta_rnpp_D(i,k,jrow) = ta_rnpp_D(i,k,jrow)
     &                              + rnpp_D(i,k,j)
                ta_rgraz_D(i,k,jrow) = ta_rgraz_D(i,k,jrow)
     &                               + rgraz_D(i,k,j)
                ta_rmorp_D(i,k,jrow) = ta_rmorp_D(i,k,jrow)
     &                               + rmorp_D(i,k,j)
                ta_rnfix(i,k,jrow) = ta_rnfix(i,k,jrow) + rnfix(i,k,j)
#  endif
              endif
# endif
!             calculate total export to get total import for next layer
              expo = expo*dzt(k)

            enddo

# if defined uvic_o2
            kmax = kmt(i,jrow)
            do k=1,kmax
!             limit oxygen consumption below concentrations of
!             5umol/kg as recommended in OCMIP
              fo2 = 0.5*tanh(t(i,k,j,io2,taum1)*1000. - 5.)
!             sink of oxygen
              so2 = src(i,k,j,ispo4)*redotp
              src(i,k,j,iso2) = -so2*(0.5 + fo2)
#  if defined uvic_nitrogen
!             add denitrification as source term for NO3
              no3flag = 0.5+sign(0.5,t(i,k,j,ino3,taum1)-trcmin)
!             800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
              deni = 800.*no3flag*so2*(0.5 - fo2)
              src(i,k,j,isno3) = src(i,k,j,isno3) - deni
#   if defined uvic_npzd_out
              rdeni(i,k,jrow) = deni
#   endif
#  endif
            enddo
# endif
# if defined time_averages && defined uvic_npzd_out

!-----------------------------------------------------------------------
!           accumulate time averages for full depth variables
!-----------------------------------------------------------------------
            kmax = kmt(i,jrow)
            do k=1,kmax
              if (timavgperts .and. .not. euler2) then
                ta_rremi(i,k,jrow) = ta_rremi(i,k,jrow) + rremi(i,k,j)
                ta_rexpo(i,k,jrow) = ta_rexpo(i,k,jrow) + rexpo(i,k,j)
#  if defined uvic_nitrogen && defined uvic_o2
                ta_rdeni(i,k,jrow) = ta_rdeni(i,k,jrow) + rdeni(i,k,j)
#  endif
              endif
            enddo
# endif
# if !defined uvic_npzd_no_calcite

!-----------------------------------------------------------------------
!           remineralize calcite
!-----------------------------------------------------------------------
            kmax = kmt(i,jrow)
            do k=1,kmax-1
#  if defined uvic_carbon
              src(i,k,j,isdic) = src(i,k,j,isdic) + prca*rcak(k)
#  endif
#  if defined uvic_alk
              src(i,k,j,isalk) = src(i,k,j,isalk) + 2.*prca*rcak(k)
#  endif
            enddo
#  if defined uvic_carbon
            src(i,kmax,j,isdic) = src(i,kmax,j,isdic) + prca*rcab(kmax)
#  endif
#  if defined uvic_alk
            src(i,kmax,j,isalk) = src(i,kmax,j,isalk)+2.*prca*rcab(kmax)
#  endif
# endif
          endif

        enddo

      enddo
#endif
#if defined uvic_carbon && defined uvic_carbon_14

!-----------------------------------------------------------------------
!     set source for c14
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
# if defined uvic_npzd
              src(i,k,j,isc14) = src(i,k,j,isdic)*rstd
     &                         - 3.836e-12*t(i,k,j,ic14,taum1)
# else
              src(i,k,j,isc14) = - 3.836e-12*t(i,k,j,ic14,taum1)
# endif
            enddo
          endif
        enddo
      enddo
#endif

!-----------------------------------------------------------------------
!     solve for one tracer at a time
!     n = 1 => temperature
!     n = 2 => salinity
!     n > 2 => other tracers (if applicable)
!-----------------------------------------------------------------------

      do n=1,nt
#if defined linearized_advection

!-----------------------------------------------------------------------
!       calculate 2*linear advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux_lin (joff, js, je, is, ie, n)
#endif
#if defined second_order_tracer_advection

!-----------------------------------------------------------------------
!       calculate 2* 2nd order advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux_2nd (joff, js, je, is, ie, n)
#endif
#if defined fct

!-----------------------------------------------------------------------
!       calculate 2* FCT tracer flux
!-----------------------------------------------------------------------

        call adv_flux_fct(joff, js, je, is, ie, n)
#endif
#if defined fourth_order_tracer_advection

!-----------------------------------------------------------------------
!       calculate 2* 4th order advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux_4th (joff, js, je, is, ie, n)
#endif
#if defined quicker

!-----------------------------------------------------------------------
!       calculate 2* 3rd order advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux_quick (joff, js, je, is, ie, n)
#endif

!-----------------------------------------------------------------------
!       calculate diffusive flux across eastern and northern faces
!       of "T" cells due to various parameterizations for diffusion.
!-----------------------------------------------------------------------

#if defined consthmix && !defined biharmonic

!       diffusive flux on eastern face of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) =
#  if defined bryan_lewis_horizontal
     &                         diff_cet(k)*cstdxur(i,j)*
#  else
     &                         ah_cstdxur(i,j)*
#  endif
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
# if defined isopycmix || defined isoneutralmix

!       diffusive flux on northern face of "T" cells
!       (background for isopycnal mixing)

        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =
#   if defined bryan_lewis_horizontal
     &                         diff_cnt(k)*
#   else
     &                         diff_cnt*
#   endif
     &           csu_dyur(jrow)*(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
# else

!       diffusive flux on northern face of "T" cells is not calculated.
!       instead, it is incorporated into DIFF_Ty for performance issues.

# endif
#endif
#if defined consthmix && defined biharmonic

!-----------------------------------------------------------------------
!       diffusive flux across eastern face of "T" cell
!       diffusive flux across northern face of "T" cell
!-----------------------------------------------------------------------

        m = 2 + n
        do j=js,je
          jrow = j + joff
          ahbi_cstr = diff_cet*cstr(jrow)
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = ahbi_cstr*dxur(i)*
     &                         (del2(i+1,k,j,m)-del2(i,k,j,m))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          ahbi_csu_dyur = diff_cnt*csu(jrow)*dyur(jrow)
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = ahbi_csu_dyur*
     &                        (del2(i,k,j+1,m) - del2(i,k,j,m))
            enddo
          enddo
        enddo
#endif
#if defined smagnlmix

!       diffusive flux on eastern and northern faces of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = diff_cet(i,k,j)*cstdxur(i,j)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = diff_cnt(i,k,j)*csu_dyur(jrow)*
     &                     (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#endif

!-----------------------------------------------------------------------
!       calculate diffusive flux across bottom face of "T" cells
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km-1
            do i=istrt,iend
              diff_fb(i,k,j) = diff_cbt(i,k,j)*dzwr(k)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo

#if defined isopycmix || defined isoneutralmix

!-----------------------------------------------------------------------
!       compute isopycnal diffusive flux through east, north,
!       and bottom faces of T cells.
!-----------------------------------------------------------------------

        call isoflux (joff, js, je, is, ie, n)
#endif

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "T" cells for diffusion
!       and advection. for isopycnal diffusion, set adiabatic boundary
!       conditions.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow   = j + joff
          do i=istrt,iend
            kb              = kmt(i,jrow)
#if defined uvic_replacst
            diff_fb(i,0,j)  = c0
#else
            diff_fb(i,0,j)  = stf(i,j,n)
#endif
            diff_fb(i,kb,j) = btf(i,j,n)
            adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                        t(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbt(i,km,j)*t(i,km,j,n,tau)
          enddo
        enddo

#if defined source_term

!-----------------------------------------------------------------------
!       set source term for "T" cells
!-----------------------------------------------------------------------

        source(:,:,:) = c0

# if defined uvic_npzd || defined uvic_carbon_14
        if (itrc(n) .ne. 0) then
          do j=js,je
            do k=1,km
              do i=istrt,iend
                source(i,k,j) = src(i,k,j,itrc(n))
              enddo
            enddo
          enddo
        endif
# endif
# if defined sponges

!-----------------------------------------------------------------------
!       add newtonian damping term for sponge layers to source
!-----------------------------------------------------------------------

        call sponge1 (joff, js, je, istrt, iend, n, t(1,1,1,1,taum1)
     &,               source)
# endif
# if defined shortwave

!-----------------------------------------------------------------------
!       incorporate short wave penetration into source
!-----------------------------------------------------------------------

        if (n .eq. 1) then
          call swflux0 (joff, js, je, istrt, iend, source)
        endif
# endif
#endif

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer using statement functions to represent
!       each component of the calculation
!-----------------------------------------------------------------------

#if defined obc_west || defined obc_east || defined obc_north
# include "tracer_obc.inc"
#else

!       1st: solve using all components which are treated explicitly

        do j=js,je
          jrow   = j + joff
          do k=1,km
            twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &          DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j)
     &          - ADV_Tx(i,k,j) -  ADV_Ty(i,k,j,jrow,n) -  ADV_Tz(i,k,j)
# if defined isopycmix && defined gent_mcwilliams && !defined fct && !defined quicker
     &          - ADV_Txiso(i,k,j,n) - ADV_Tyiso(i,k,j,jrow,n)
     &          - ADV_Tziso(i,k,j)
# endif
# if defined source_term
     &          + source(i,k,j)
# endif
# if defined llnl_plume
     &          + work(i,k,j)*subflux(i,jrow,n)
# endif
     &          )*tmask(i,k,j)
            enddo
          enddo
        enddo
# if defined implicitvmix || defined isopycmix || defined isoneutralmix || defined redi_diffusion

!       2nd: add in portion of vertical diffusion handled implicitly

        call ivdift (joff, js, je, istrt, iend, n, twodt)
# endif
#endif

        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, js, je, istrt, iend, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

#if defined uvic_replacst

      gamma = zw(1)/twodt(1)
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          stf(i,j,1) = gamma*(sbc(i,jrow,isst) - t(i,1,j,1,taup1))
          stf(i,j,2) = gamma*(sbc(i,jrow,isss) - t(i,1,j,2,taup1))
          t(i,1,j,1,taup1) = sbc(i,jrow,isst)
          t(i,1,j,2,taup1) = sbc(i,jrow,isss)
        enddo
      enddo
#endif
#if defined uvic_convect_brine

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!     convect brine rejected under all ice categories
!-----------------------------------------------------------------------

      call convect_brine (joff, js, je, is, ie)
#else
# if !defined implicitvmix || defined isopycmix || defined isoneutralmix

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!-----------------------------------------------------------------------

#  if defined fullconvect
      call convct2 (t(1,1,1,1,taup1), joff, js, je, is, ie, kmt)
      do j=js,je
        do n=1,nt
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
#  else
      call convct (t(1,1,1,1,taup1), ncon, joff, js, je, istrt, iend
     &,            kmt)
#  endif
# endif
#endif
#if defined uvic_save_convection
      if (timavgperts .and. eots) then
        if (joff .eq. 0) nta_conv = nta_conv + 1
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            ta_totalk(i,jrow) = ta_totalk(i,jrow) + totalk(i,j)
            ta_vdepth(i,jrow) = ta_vdepth(i,jrow) + vdepth(i,j)
            ta_pe(i,jrow) = ta_pe(i,jrow) + pe(i,j)
          enddo
        enddo
      endif
#endif

!-----------------------------------------------------------------------
!     construct diagnostics after convection
!-----------------------------------------------------------------------

      idiag = 10
      call diagt2 (joff, js, je, istrt, iend, idiag)

#if defined fourfil || defined firfil

!-----------------------------------------------------------------------
!     filter tracers at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filt (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in tracer'
        stop '=>tracer'
      endif
#endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics after filtering (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, js, je, istrt, iend, idiag)

#if !defined uvic_replacst

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------

      call asbct (joff, js, je, istrt, iend, isst, itemp)
      call asbct (joff, js, je, istrt, iend, isss, isalt)
#endif
#if defined uvic_carbon
      call asbct (joff, js, je, istrt, iend, issdic, idic)
# if defined uvic_carbon_14
      call asbct (joff, js, je, istrt, iend, issc14, ic14)
# endif
#endif
#if defined uvic_alk
      call asbct (joff, js, je, istrt, iend, issalk, ialk)
#endif
#if defined uvic_o2
      call asbct (joff, js, je, istrt, iend, isso2, io2)
#endif
#if defined uvic_npzd
      call asbct (joff, js, je, istrt, iend, isspo4, ipo4)
# if defined uvic_npzd_vflux
      call asbct (joff, js, je, istrt, iend, issphyt, iphyt)
      call asbct (joff, js, je, istrt, iend, isszoop, izoop)
      call asbct (joff, js, je, istrt, iend, issdetr, idetr)
# endif
# if defined uvic_nitrogen
      call asbct (joff, js, je, istrt, iend, issno3, ino3)
#  if defined uvic_npzd_vflux
      call asbct (joff, js, je, istrt, iend, issdiaz, idiaz)
#  endif
# endif
#endif
#if defined uvic_cfc11
      call asbct (joff, js, je, istrt, iend, isscfc11, icfc11)
#endif
#if defined uvic_cfc12
      call asbct (joff, js, je, istrt, iend, isscfc12, icfc12)
#endif
#if defined uvic_carbon && defined uvic_carbon_14

!-----------------------------------------------------------------------
!     calculate diagnostic delta carbon 14
!-----------------------------------------------------------------------

      if (tsiperts .or. timavgperts .or. snapts) then
        rrstd = 1000./rstd
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              dc14(i,k,j) = (rrstd*t(i,k,j,ic14,taup1)
     &                      /(t(i,k,j,idic,taup1) + epsln) - 1000.)
     &                      *tmask(i,k,j)
            enddo
          enddo
        enddo
      endif
      if (tsiperts .and. eots) then
        if (js+joff .eq. 2) dc14bar = 0.
        do j=js,je
          jrow = j + joff
          fy = cst(jrow)*dyt(jrow)
          do k=1,km
            fyz = fy*dzt(k)
            do i=istrt,iend
              dc14bar = dc14bar + dc14(i,k,j)*dxt(i)*fyz*tmask(i,k,j)
            enddo
          enddo
        enddo
      endif
      if (timavgperts .and. .not. euler2) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              ta_dc14(i,k,jrow) = ta_dc14(i,k,jrow) + dc14(i,k,j)
            enddo
          enddo
        enddo
      endif
#endif

#if defined trace_indices && !defined obc_west && !defined obc_east
      write (stdout,'(2x,5(a,i4))')
     & "=> In tracer: js=",js," je=",je," joff=",joff
     &,", calculating latitude rows: ",js+joff," to ",je+joff
#endif

      return
      end

      subroutine diagt1 (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     construct diagnostics associated with tracer component "n"

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = (1,2) = (u,v) velocity component
!       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

#include "param.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#if defined tracer_averages
# include "ctavg.h"
#endif
#include "diag.h"
#include "diaga.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#if defined isopycmix || defined isoneutralmix
# include "isopyc.h"
#endif
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
#if defined meridional_tracer_budget
# include "ctmb.h"
      dimension stor(imt,km), div(imt,km), sorc(imt,km), dif(imt,km)
      dimension t1(imt), t2(imt), t3(imt), t4(imt)
#endif
#if defined time_step_monitor
      dimension temp1(imt,km), temp2(imt,km), temp3(imt,km)
#endif
      dimension twodt(km)
#include "fdift.h"

#if defined save_mixing_coeff

!-----------------------------------------------------------------------
!     diagnostic: estimate mixing coefficients on east, north, and
!                 bottom face of T cells from the flux

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      if (cmixts .and. n .eq. 1 .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
              dtdx = (t(i+1,k,j,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i+1,k,j)*tmask(i,k,j)
     &                *cstr(jrow)*dxur(i) + epsln
              ce(i,k,j,2) = diff_fe(i,k,j)/dtdx
# if !defined consthmix || defined biharmonic || defined isopycmix || defined isoneutralmix
              dtdy = (t(i,k,j+1,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i,k,j+1)*tmask(i,k,j)
     &                *dyur(jrow) + epsln
              cn(i,k,j,2) = csur(jrow)*diff_fn(i,k,j)/dtdy
# else
              cn(i,k,j,2) = csur(jrow)*ah
# endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              dtdz = (t(i,k,j,1,taum1)-t(i,k+1,j,1,taum1))
     &                *tmask(i,k,j)*tmask(i,k+1,j)
     &                *dzwr(k) + epsln
              cb(i,k,j,2) = diff_fb(i,k,j)/dtdz
# if defined isopycmix || defined isoneutralmix
     &                      + diff_fbiso(i,k,j)/dtdz
# endif
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,2) = 0.0
          enddo
        enddo

        do j=js,je
          call setbcx (ce(1,1,j,2), imt, km)
          call setbcx (cn(1,1,j,2), imt, km)
          call setbcx (cb(1,1,j,2), imt, km)
        enddo
      endif
#endif

# if defined save_convection

!-----------------------------------------------------------------------
!     diagnostic: initialize temperature before convection

!     based on code by: S. Griffies
!-----------------------------------------------------------------------

      if (exconvts .and. n .eq. 1 .and. eots) then
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv0(i,k,j) = t(i,k,j,1,taup1)
            enddo
          enddo
        enddo
      endif

#endif

#if defined time_step_monitor

!-----------------------------------------------------------------------
!     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
!                 globally
!     based on code by: R. C. Pacanowski
!               (based on diagnostic by M. Cox)
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              darea      = dzt(k)*dxt(i)*cosdyt*tmask(i,k,j)
              temp3(i,k) = t(i,k,j,n,tau)*darea
              temp1(i,k) = t(i,k,j,n,tau)**2*darea
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     darea*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo
      endif
#endif

#if defined tracer_averages

!-----------------------------------------------------------------------
!     diagnostic: accumulate tracers for averages under horizontal
!                 regions (use units of meters, rather than cm)

!     based on code by: K. Dixon
!-----------------------------------------------------------------------

      if (tavgts .and. eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*dzt(k)*tmask(i,k,j)*0.01
             enddo
            endif
          enddo
        enddo
      endif
#endif

#if defined tracer_yz

!----------------------------------------------------------------------
!     diagnostic: integrate tracer equations in longitude

!     based on code by: R. C. Pacanowski
!----------------------------------------------------------------------

      if (tyzts .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            rtwodt = c1/twodt(k)
            sumdx = c0
            do m=1,5
              tyz(jrow,k,n,m) = c0
            enddo
            do i=is,ie
              delx  = dxt(i)*tmask(i,k,j)
              sumdx = sumdx + delx
              tyz(jrow,k,n,1) = tyz(jrow,k,n,1)+ delx*t(i,k,j,n,tau)
              tyz(jrow,k,n,2) = tyz(jrow,k,n,2) + rtwodt*delx*
     &          (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              tyz(jrow,k,n,3) =  tyz(jrow,k,n,3)
     &          - delx*(ADV_Ty(i,k,j,jrow,n) + ADV_Tz(i,k,j))
              tyz(jrow,k,n,4) =  tyz(jrow,k,n,4)
     &          + delx*(DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j))
# if defined source_term
              tyz(jrow,k,n,5) = tyz(jrow,k,n,5) + delx*source(i,k,j)
# endif
            enddo
            sumdxr = c1/(sumdx + epsln)
            do m=1,5
              tyz(jrow,k,n,m) = tyz(jrow,k,n,m)*sumdxr
            enddo
          enddo
        enddo
        if (js+joff .eq. 2) then
          do m=1,5
            do k=1,km
              tyz(1,k,n,m)   = c0
              tyz(jmt,k,n,m) = c0
            enddo
          enddo
        endif
      endif
#endif

#if defined meridional_tracer_budget

!----------------------------------------------------------------------
!     diagnostic: integrate equations in depth, lon, and time for
!                 each basin and jrow. basin # 0 is used to catch
!                 values in land areas

!     based on code by: R. C. Pacanowski
!----------------------------------------------------------------------

      if (tmbperts .and. eots) then
        do j=js,je
          jrow = j + joff
          if (jrow .eq. 2 .and. n .eq. 1) numtmb = numtmb + 1
          cosdyt = cst(jrow)*dyt(jrow)
          do k=1,km
            do i=is,ie
              stor(i,k) = 0.0
              div(i,k)  = 0.0
              sorc(i,k) = 0.0
              dif(i,k)  = 0.0
            enddo
          enddo
          do k=1,km
            rtwodt = c1/twodt(k)
            do i=is,ie
              dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
              dxdydz          = dxdy*dzt(k)
              stor(i,k) = stor(i,k) + rtwodt*dxdydz*
     &                            (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              div(i,k)  = div(i,k)  - dxdydz*ADV_Ty(i,k,j,jrow,n)
# if defined source_term
              sorc(i,k) = sorc(i,k) + dxdydz*source(i,k,j)
# endif
              dif(i,k)  = dif(i,k)  + dxdydz*DIFF_Ty(i,k,j,jrow,n)
            enddo
          enddo

!         the accumulation is done this way for issues of speed

          do i=is,ie
            t1(i) = stor(i,1)
            t2(i) = div(i,1)
            t3(i) = sorc(i,1)
            t4(i) = dif(i,1)
          enddo

          do k=2,km
            do i=is,ie
              t1(i) = t1(i) + stor(i,k)
              t2(i) = t2(i)  + div(i,k)
              t3(i) = t3(i) + sorc(i,k)
              t4(i) = t4(i)  + dif(i,k)
            enddo
          enddo

          do i=is,ie
            m               = msktmb(i,jrow)
            tstor(jrow,n,m) = tstor(jrow,n,m) + t1(i)
            tdiv(jrow,n,m)  = tdiv(jrow,n,m)  + t2(i)
            tsorc(jrow,n,m) = tsorc(jrow,n,m) + t3(i)
            tdif(jrow,n,m)  = tdif(jrow,n,m)  + t4(i)
          enddo

          k = 1
          do i=is,ie
            m               = msktmb(i,jrow)
            dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
            tflux(jrow,n,m) = tflux(jrow,n,m) + dxdy*stf(i,j,n)
          enddo
          if (n .eq. 1) then
            do k=1,km
              do i=is,ie
                m               = msktmb(i,jrow)
                dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
                dxdydz          = dxdy*dzt(k)
                smdvol(jrow,m)  = smdvol(jrow,m)  + dxdydz
              enddo
            enddo
          endif
        enddo
      endif
#endif

#if defined gyre_components

!-----------------------------------------------------------------------
!     diagnostic: compute the northward transport components of
!                 each tracer

!     based on code by: R. C. Pacanowski
!               (based on diagnostic by K. Bryan)
!-----------------------------------------------------------------------

      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)
#endif

#if defined term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the tracer equations
!                 over specified regional volumes.

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)
#endif

#if defined xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt1 (joff, js, je, n)
#endif
#if defined uvic_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations
!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt1 (joff, js, je, n)
#endif

      return
      end

      subroutine diagt2 (joff, js, je, is, ie, idiag)

!-----------------------------------------------------------------------
!     construct d(tracer)/dt diagnostics

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       idiag = 1  => total tracer change
!       idiag = 10 => change of tracer due to filtering(also convection)
!-----------------------------------------------------------------------

#include "param.h"
#include "coord.h"
#include "diaga.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "timeavgs.h"

#if defined save_convection

!-----------------------------------------------------------------------
!     diagnostic: save tracer time change due to explicit convection
!     idiag = 10 signifies diagnostics immediately after convection

!     based on code by: S. Griffies
!-----------------------------------------------------------------------

      if (exconvts .and. idiag .eq. 10 .and. eots) then
        rdt = c1/c2dtts

!       excnv1 = epsln over land cells and d(convect)/dt over ocean

        do j=js,je
          do k=1,km
            do i=1,imt
              excnv1(i,k,j) = tmask(i,k,j)*
     &         (t(i,k,j,1,taup1)-excnv0(i,k,j))*rdt
     &                       + (c1-tmask(i,k,j))*epsln
            enddo
          enddo
        enddo

        reltim = relyr
        if (joff + js .eq. 2) then
          write (stdout,*) ' =>Writing explicit convection at ts=',itt
     &  , ' ',stamp
          call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

          period = 0.0
          iotext = 'read(iocv) reltim, period, imt, jmt, km, flag'
          write (iocv) stamp, iotext, expnam
          write (iocv) reltim, period, imt, jmt, km, epsln

          iotext = 'read(iocv) (xt(i),i=1,imt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, xt, imt)

          iotext = 'read(iocv) (yt(j),j=1,jmt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, yt, jmt)

          iotext = 'read(iocv) (zt(k),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, zt, km)

          call relunit (iocv)
        endif

        call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

        do j=js,je
          jrow = j+joff

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocv) (convect(i,k),i=1,imt),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, excnv1(1,1,j), imt*km)

        enddo
        call relunit(iocv)
      endif
#endif

#if defined term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate d/dt(tracer) over specified regional volumes
!                  after convection and filtering

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)
#endif

#if defined xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt2 (joff, js, je, idiag)
#endif
#if defined uvic_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering

!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt2 (joff, js, je, idiag)
#endif

      return
      end

      subroutine asbct (joff, js, je, is, ie, isbc, itr)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C. (surface boundary conditions)

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       isbc  = index for sbc
!       itr   = index for tracer

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

#include "param.h"
#include "csbc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"

!     initialize the Atmos S.B.C. at the start of each ocean segment
!     (do not alter values in land)

      if (isbc .le. 0 .or. itr .le. 0) return

      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) sbc(i,jrow,isbc) = c0
          enddo
        enddo
      endif

!     accumulate surface tracers for the Atmos S.B.C. every time step

      if (eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            sbc(i,jrow,isbc) = sbc(i,jrow,isbc)+t(i,1,j,itr,taup1)
          enddo
        enddo
      endif

!     average the surface tracers for the Atmos S.B.C. at the end of
!     each ocean segment. (do not alter values in land)

      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0)
     &        sbc(i,jrow,isbc) = rts*sbc(i,jrow,isbc)
          enddo
        enddo
      endif

      return
      end

#if defined implicitvmix || defined isopycmix || defined isoneutralmix || defined redi_diffusion
      subroutine ivdift (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     solve vertical diffusion of tracers implicitly

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

# include "param.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt(km)

!     store terms to compute implicit vertical mixing on
!     diagnostic time steps

# if defined xbts || defined uvic_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
# else
#  if defined term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#  endif
# endif

# if defined tcvmix
      call invtri (tp1(1,1,1,n), stf(1,1,n), btf(1,1,n), vdca(???)
     &,              twodt, kmt, tmask(1,1,1), is, ie, joff, js, je)
      if (imt .gt. 1) then
        print *,'Error: ivdif.F is not converted for ifvef tvcmi x'
        stop '=>ivdif'
      endif
# else
      call invtri (t(1,1,1,n,taup1), stf(1,1,n), btf(1,1,n)
     &, diff_cbt(1,1,jsmw), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je)
# endif

!     compute residual implicit vertical mixing

# if defined xbts || defined uvic_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
# else
#  if defined term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
#  endif
# endif

      return
      end
#endif

#if defined sponges
      subroutine sponge1 (joff, js, je, is, ie, n, tm1, tsource)

!=======================================================================
!     newtonian damping variables for sponge regions adjacent to
!     artificial southern and northern boundaries in limited domain
!     basins. data must be prepared using the "sponge" routines
!     included in the programs for working with the MOM dataset.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       n    = (1,2) = tracer component
!       tm1  = tracer at "tau-1"

!     output:
!       tsource = newtonian damping term

!     based on code by: R. C. Pacanowski
!=======================================================================

# include "param.h"
# include "iounit.h"
# include "calendar.h"
# if defined equatorial_thermocline
#  include "coord.h"
#  include "mw.h"
# endif
# include "sponge.h"
# include "switch.h"
# include "tmngr.h"

      dimension tsource(imt,km,jsmw:jemw)
      dimension tm1(imt,km,jmw,nt)

# if defined equatorial_thermocline

!       only appropriate for idealized equatorial studies

!       spng_width = width (deg) of sponges
!       spng_damp  = reciprocal of damping period (1/sec)

        spng_width = 3.0
        spng_damp  = 1.0/(5.0*daylen)

        do j=js,je
          jrow = joff + j
          if ((yt(jrow) - yu(1)) .le. spng_width
     &      .or. (yu(jmt-1) - yt(jrow)) .le. spng_width) then
            do k=1,km
              do i=is,ie
                tsource(i,k,j) = tsource(i,k,j) -
     &                  spng_damp*(tm1(i,k,j,n) - tbarz(k,n))
              enddo
            enddo
          endif
        enddo
# else

      if ((joff + js .eq. 2) .and. n .eq. 1 .and. is .eq. 2) then

!-----------------------------------------------------------------------
!       decide whether to read sponge data or not
!-----------------------------------------------------------------------

        begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
        method = 3
        call timeinterp (begtim, indxsp, tspng, spgdpm, 12, .true.
     &,           method, inextd, iprevd, wprev, readsp, inext, iprev)

!       read in the next data record from disk when needed

        if (readsp) then
          call getunit (ionew2, 'sponges', opt_sponge)
          read (ionew2, rec=inextd) stnext, spdpmn, im, kk, jm, spngs
     &,     spngn, spbuf1
          do in=1,4
            do k=1,km
              do i=1,imt
                spbuf(i,k,in,inext) = spbuf1(i,k,in)
              enddo
           enddo
          enddo
          write (stdout,'(/a,i3,a,i2,a,i2,a,a/)')
     &    '=> read sponge record =',inextd,' into buffer =', inext
     &,   ' method #',method,' at ', stamp
          call relunit (ionew2)
        endif
      endif
      if (n .le. 2) then

!-----------------------------------------------------------------------
!       construct newtonian damping term using sponge data
!-----------------------------------------------------------------------

        do j=js,je
          jrow = joff + j
          if (spngs(jrow) .ne. c0) then
            tnext = c1-wprev
            do k=1,km
              do i=is,ie
                data = tnext*spbuf(i,k,n,inext)
     &               + wprev*spbuf(i,k,n,iprev)
                tsource(i,k,j) = tsource(i,k,j) -
     &                            spngs(jrow)*(tm1(i,k,j,n) - data)
              enddo
            enddo
          endif

          if (spngn(jrow) .ne. c0) then
            tnext = c1-wprev
            do k=1,km
              do i=is,ie
                data = tnext*spbuf(i,k,n+2,inext)
     &               + wprev*spbuf(i,k,n+2,iprev)
                tsource(i,k,j) = tsource(i,k,j) -
     &                            spngn(jrow)*(tm1(i,k,j,n) - data)
              enddo
            enddo
          endif
        enddo
      endif
# endif

      return
      end
#endif

#if defined shortwave
      subroutine swflux0 (joff, js, je, is, ie, source)

!=======================================================================
!       incorporate short wave penetration via the "source"
!       term. note that the divergence of shortwave for the first
!       level "divpen(1)" is compensating for the effect of having
!       the shortwave component already included in the total
!       surface tracer flux "stf(i,j,1)"

!       incorporating the penetrative shortwave radiative flux into
!       the vertical diffusive flux term directly is not correct when
!       vertical diffusion is solved implicitly.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

!     output:
!       source = source term penetrative short wave

!     based on code by: R. C. Pacanowski
!=======================================================================

# include "param.h"
# include "csbc.h"
# include "cshort.h"

      dimension source(imt,km,jsmw:jemw)

      if (isw .ne. 0) then
        do j=js,je
          jrow   = j + joff
          do k=1,km
            do i=is,ie
              source(i,k,j) = source(i,k,j)
     &                       + sbc(i,jrow,ipsw)*divpen(k)
            enddo
          enddo
        enddo
      endif

      return
      end
#endif
