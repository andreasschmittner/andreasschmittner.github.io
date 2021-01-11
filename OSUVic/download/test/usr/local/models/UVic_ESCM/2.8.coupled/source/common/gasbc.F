      subroutine gasbc (is, ie, js, je)

!=======================================================================
!     calculate boundary conditions for the atmospheric model

!     based on code by: A. Fanning and M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "param.h"
# include "coord.h"
# include "csbc.h"
# include "ice.h"
# include "switch.h"
# include "tmngr.h"
# include "cembm.h"
# include "atm.h"
# include "insolation.h"
# include "calendar.h"
# include "grdvar.h"
# include "levind.h"
# if defined uvic_embm_adv_q || defined uvic_embm_adv_t
#  include "solve.h"
# endif
# if defined uvic_mtlm
#  include "mtlm.h"
# endif

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1, l

      real sss, sst, phlo, phhi, sit_in, pt_in, atmpres, xconv, ta_in
      real dic_in, ph, co2_in, co2star, dco2star, dco2surf, dpco2
      real pco2surf, scco2, piston_vel, avgflxc, fa, calday, f, sco2
      real o2sat, o2sato, o2surf, piston_o2, cfc11ccn, cfc12ccn, wt
      real sccfc, piston_cfc, sol_cfc, cfcsat, ao, tarea, tdc14ccn
      real h_r, d, f1, f2, f3, f4, f5
# if !defined uvic_embm_annual
      real cosz(is:ie,js:je)
# endif
      real dmsk(is:ie,js:je)

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

!     xconv is constant to convert piston_vel from cm/hr -> cm/s
!     here it is 100.*a*xconv (100 => m to cm, a=0.337, xconv=1/3.6e+05)
      xconv = 33.7/3.6e+05
# if defined uvic_carbon
      phlo = 6.
      phhi = 10.
      sit_in = 7.6875e-03 !mol/m^3
      pt_in = 0.5125e-3   !mol/m^3
      atmpres = 1.0       !atm
      co2flx = 0.
!     fa is used in converting g carbon cm-2 => ppmv CO2
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      fa = 1./(4.138e-7*rhoatm*shc)
# endif
# if defined uvic_mom
#  if defined uvic_newflxref

!-----------------------------------------------------------------------
!     calculate new global average sea surface flux references
!-----------------------------------------------------------------------
      dmsk(:,:) = 1.
      where (kmt(:,:) .eq. 0) dmsk(:,:) = 0.
      if (isss .ne. 0) then
        call areaavg (sbc(1,1,isss), dmsk, socn)
        socn = socn + 0.035
      endif
      if (issdic .ne. 0) call areaavg (sbc(1,1,issdic), dmsk, dicocn)
      if (issalk .ne. 0) call areaavg (sbc(1,1,issalk), dmsk, alkocn)
      if (isso2 .ne. 0) call areaavg (sbc(1,1,isso2), dmsk, o2ocn)
      if (isspo4 .ne. 0) call areaavg (sbc(1,1,isspo4), dmsk, po4ocn)
      if (issphyt .ne. 0) call areaavg (sbc(1,1,issphyt), dmsk, phytocn)
      if (isszoop .ne. 0) call areaavg (sbc(1,1,isszoop), dmsk, zoopocn)
      if (issdetr .ne. 0) call areaavg (sbc(1,1,issdetr), dmsk, detrocn)
      if (issno3 .ne. 0) call areaavg (sbc(1,1,issno3), dmsk, no3ocn)
      if (issdiaz .ne. 0) call areaavg (sbc(1,1,issdiaz), dmsk, diazocn)
      if (issc14 .ne. 0) call areaavg (sbc(1,1,issc14), dmsk, c14ocn)
#  endif
#  if defined uvic_cfc || defined uvic_cfc11 || defined uvic_cfc12
!     always update cfc flux references since initial is usually zero
      if (isscfc11 .ne. 0) call areaavg (sbc(1,1,isscfc11), dmsk
     &, cfc11ocn)
      if (isscfc12 .ne. 0) call areaavg (sbc(1,1,isscfc12), dmsk
     &, cfc12ocn)
#  endif
# endif

!-----------------------------------------------------------------------
!     zero totals for new accumulation
!-----------------------------------------------------------------------
      atatm = 0.
      flux(:,:,:) = 0.
# if defined llnl_plume
      subflux(:,:,:) = 0.
# endif
# if defined uvic_convect_brine
      cbf(:,:,:) = 0.
      cba(:,:,:) = 0.
# endif
      sbc(:,:,ihflx) = 0.
      sbc(:,:,isflx) = 0.
      sbc(:,:,iro) = 0.
# if defined uvic_mtlm
      sbc(:,:,iat) = 0.
      sbc(:,:,irh) = 0.
      sbc(:,:,ipr) = 0.
      sbc(:,:,ips) = 0.
      sbc(:,:,iaws) = 0.
      sbc(:,:,iswr) = 0.
# endif
# if defined uvic_carbon
      sbc(:,:,idicflx) = 0.
#  if defined uvic_carbon_14
      sbc(:,:,ic14flx) = 0.
#  endif
# endif
# if defined uvic_alk
      sbc(:,:,ialkflx) = 0.
# endif
# if defined uvic_o2
      sbc(:,:,io2flx) = 0.
# endif
# if defined uvic_npzd
      sbc(:,:,ipo4flx) = 0.
#  if defined uvic_npzd_vflux
      sbc(:,:,iphytflx) = 0.
      sbc(:,:,izoopflx) = 0.
      sbc(:,:,idetrflx) = 0.
#  endif
#  if defined uvic_nitrogen
      sbc(:,:,ino3flx) = 0.
#   if defined uvic_npzd_vflux
      sbc(:,:,idiazflx) = 0.
#   endif
#  endif
# endif
# if defined uvic_cfc11
      sbc(:,:,icfc11flx) = 0.
# endif
# if defined uvic_cfc12
      sbc(:,:,icfc12flx) = 0.
# endif
# if defined uvic_embm_solardata || defined uvic_embm_solardata_transient

!-----------------------------------------------------------------------
!     set solar constant
!-----------------------------------------------------------------------
      call solardata
# endif
# if defined uvic_embm_co2data || defined uvic_embm_co2data_transient

!-----------------------------------------------------------------------
!     set co2 concentration
!-----------------------------------------------------------------------
      call co2data
# endif
# if defined uvic_carbon_14
#  if defined uvic_embm_c14data || defined uvic_embm_c14data_transient

!-----------------------------------------------------------------------
!     set c14 concentration
!-----------------------------------------------------------------------
      call c14data
#   if defined uvic_embm_c14data
      tdc14ccn = 0.
      tarea = 0.
#   endif
#  endif
# endif
# if defined uvic_cfc11 || defined uvic_cfc12

!-----------------------------------------------------------------------
!     set CFC concentration
!-----------------------------------------------------------------------
      call cfcdata
# endif
# if !defined uvic_embm_annual

!-----------------------------------------------------------------------
!     update insolation for the current day
!-----------------------------------------------------------------------
!     subroutine decl is expecting a 365.25 day year
      calday = dayoyr*365.25/yrlen
      call decl (calday, eccen, obliq, mvelp, lambm0, sindec, eccf)
      i = (ie-is+1)*(je-js+1)
      call zenith (i, c0, daylen, daylen, tlat, tlon, sindec, cosz)
      solins(:,:) = solarconst*eccf*cosz(:,:)
# endif

!-----------------------------------------------------------------------
!     update any atmospheric data
!-----------------------------------------------------------------------
      call atmos
# if defined uvic_embm_awind

!-----------------------------------------------------------------------
!     calculate winds with new feedback
!-----------------------------------------------------------------------
      call add_awind (is, ie, js, je)
# endif

!-----------------------------------------------------------------------
!     calculate freezing point of sea water using UNESCO (1983)
!-----------------------------------------------------------------------

      do j=jsp1,jem1
        do i=isp1,iem1
          if (tmsk(i,j) .ge. 0.5) then
            sss = 1000.0*sbc(i,j,isss) + 35.0
            frzpt(i,j) = -.0575*sss + 1.71e-3*sss**1.5 - 2.155e-4*sss**2
# if defined uvic_mom
#  if defined uvic_carbon || defined uvic_o2 || defined uvic_cfc11 || defined uvic_cfc12
            sst = sbc(i,j,isst)
#   if defined uvic_ice
            ao = 1. - aice(i,j,1)
#   else
            ao = 1.
#   endif
#  endif
#  if defined uvic_carbon

!-----------------------------------------------------------------------
!           calculate ocean carbon fluxes
!-----------------------------------------------------------------------
#   if defined uvic_alk
            ta_in = sbc(i,j,issalk)
#   else
            ta_in = 2.36775*sss/(socn*1000.)
#   endif
            dic_in = sbc(i,j,issdic)
#   if defined uvic_carbon_co2_2d
            co2_in = at(i,j,2,ico2)
#   else
            co2_in = co2ccn
#   endif
            call co2calc_SWS (sst, sss, dic_in, ta_in, pt_in, sit_in
     &,                       phlo, phhi, ph, co2_in, atmpres, co2star
     &,                       dco2star, pCO2surf, dpco2)
!           Schmidt number for CO2
            scco2 = 2073.1 - 125.62*sst + 3.6276*sst**2
     &            - 0.043219*sst**3
            piston_vel = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *((scco2/660.)**(-0.5))
!           dic in umol cm-3 or (mol m-3) => flux in umol cm-2 s-1
            sbc(i,j,idicflx) = piston_vel*dco2star
#   if defined uvic_carbon_co2_2d
!           convert from umol cm-2 s-1 => ppmv s-1
            flux(i,j,ico2) = sbc(i,j,idicflx)*12.e-6*fa + flux(i,j,ico2)
#   endif
#   if defined uvic_carbon_14

!-----------------------------------------------------------------------
!           calculate ocean c14 fluxes
!-----------------------------------------------------------------------
#    if defined uvic_embm_c14data
            if (tlat(i,j) .gt. 20.) then
              dc14ccn = dc14ccnn
            elseif (tlat(i,j) .lt. -20.) then
              dc14ccn = dc14ccns
            else
              dc14ccn = dc14ccne
            endif
            tarea = tarea + dxt(i)*dyt(j)*cst(j)
            tdc14ccn = tdc14ccn + dc14ccn*dxt(i)*dyt(j)*cst(j)
#    endif

            sbc(i,j,ic14flx) = piston_vel*((dco2star + co2star)
     &                         *(1 + dc14ccn*0.001)*rstd
     &                       - co2star*sbc(i,j,issc14)/sbc(i,j,issdic))
#   endif
#  endif
#  if defined uvic_o2

!-----------------------------------------------------------------------
!           calculate ocean oxygen fluxes
!-----------------------------------------------------------------------
!           Schmidt number for O2
            sco2 = 1638.0 - 81.83*sst + 1.483*sst**2 - 0.008004*sst**3
!           piston velocity for O2
            piston_o2 = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *(sco2/660.0)**(-0.5)
!           oxygen saturation concentration [mol/m^3]
            f1 = alog((298.15 - sst)/(273.15 + sst))
            f2 = f1*f1
            f3 = f2*f1
            f4 = f3*f1
            f5 = f4*f1
            o2sat = exp (2.00907 + 3.22014*f1 + 4.05010*f2
     &             + 4.94457*f3 - 2.56847E-1*f4 + 3.88767*f5
     &             + sss*(-6.24523e-3 - 7.37614e-3*f1 - 1.03410e-2*f2
     &             - 8.17083E-3*f3) - 4.88682E-7*sss*sss)
!           Convert from ml/l to mol/m^3
            o2sat = o2sat/22391.6*1000.0
            sbc(i,j,io2flx) = piston_o2*(o2sat - sbc(i,j,isso2))
#  endif
#  if defined uvic_cfc11

!-----------------------------------------------------------------------
!           calculate ocean CFC11 fluxes
!-----------------------------------------------------------------------
            if (tlat(i,j) .gt. 10.) then
              cfc11ccn = cfc11ccnn
            elseif (tlat(i,j) .lt. -10.) then
              cfc11ccn = cfc11ccns
            else
              wt = (tlat(i,j) + 10.)/20.
              cfc11ccn = cfc11ccnn*wt + cfc11ccns*(1. - wt)
            endif
!           Schmidt number for CFC11
            sccfc = 3501.8 -210.31*sst + 6.1851*sst**2 -0.07513*sst**3
!           piston velocity for CFC
            piston_cfc = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                 *(sccfc/660.0)**(-0.5)
!           cfc saturation concentration [mol/m^3]
            f1 = (sst + 273.16)*0.01
            d = (0.091459 - 0.0157274*f1)*f1 - 0.142382
            sol_cfc = exp(-229.9261 + 319.6552/f1 + 119.4471*alog(f1)
     $              - 1.39165*f1*f1 + sss*d )
!           conversion from mol/(l * atm) to mol/(m3 * pptv)

            cfcsat =  1.0e-12 *1000.*sol_cfc*cfc11ccn
            sbc(i,j,icfc11flx) = piston_cfc*(cfcsat - sbc(i,j,isscfc11))
#  endif
#  if defined uvic_cfc12

!-----------------------------------------------------------------------
!           calculate ocean CFC12 fluxes
!-----------------------------------------------------------------------
            if (tlat(i,j) .gt. 10.) then
              cfc12ccn = cfc12ccnn
            elseif (tlat(i,j) .lt. -10.) then
              cfc12ccn = cfc12ccns
            else
              wt = (tlat(i,j) + 10.)/20.
              cfc12ccn = cfc12ccnn*wt + cfc12ccns*(1. - wt)
            endif
!           Schmidt number for CFC12
            sccfc = 3845.4 -228.95*sst + 6.1908*sst**2 -0.067430*sst**3
 !          piston velocity for CFC12
            piston_cfc = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                 *(sccfc/660.0)**(-0.5)
!           cfc saturation concentration [mol/m^3]
            f1 = (sst + 273.16)*0.01
            d = (0.091015 - 0.0153924*f1)*f1 - 0.143566
            sol_cfc = exp(-218.0971 + 298.9702/f1 + 113.8049*alog(f1)
     $              - 1.39165*f1*f1 + sss*d )
!           conversion from mol/(l * atm) to mol/(m3 * pptv)
            cfcsat =  1.0e-12 *1000.*sol_cfc*cfc12ccn
            sbc(i,j,icfc12flx) = piston_cfc*(cfcsat - sbc(i,j,isscfc12))
#  endif
# endif
# if defined uvic_carbon && defined uvic_mtlm

          else
!-----------------------------------------------------------------------
!           calculate land carbon fluxes
!-----------------------------------------------------------------------
#  if defined uvic_carbon_co2_2d
!           convert from kg m-2 s-1 => ppmv s-1
            flux(i,j,ico2) = (sbc(i,j,inpp) - sbc(i,j,isr))*0.1*fa
     &                     + flux(i,j,ico2)
#  else
!           convert from kg m-2 s-1 => umol cm-2 s-1
            sbc(i,j,idicflx) = (sbc(i,j,inpp) - sbc(i,j,isr))*0.1/12.e-6
#  endif
#  if defined uvic_carbon_14
!           use the carbon flux scaled by rstd for c14
            sbc(i,j,ic14flx) = (sbc(i,j,inpp) - sbc(i,j,isr))
     &                         *rstd*0.1/12.e-6
#  endif
# endif
          endif
        enddo
      enddo
# if defined uvic_carbon

!-----------------------------------------------------------------------
!     set boundary conditions for carbon
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,idicflx), imt, jmt)
#  if defined uvic_carbon_co2_2d
      call setbcx (flux(1,1,ico2), imt, jmt)
#  elif defined uvic_carbon_coupled
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,idicflx), dmsk, avgflxc)
      co2ccn = co2ccn - (co2flx + avgflxc)*segtim*12.e-6*daylen*fa
#  endif
#  if defined uvic_carbon_14
#   if defined uvic_carbon_14_coupled
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,ic14flx), dmsk, avgflxc)
      c14ccn = c14ccn + (c14prod - avgflxc)*segtim*12.e-6*daylen*fa
!     calculate dc14ccn from c14ccn and co2ccn
      dc14ccn = 1000.*(c14ccn/co2ccn/rstd - 1.)
#   endif

!-----------------------------------------------------------------------
!     set boundary conditions for c14
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,ic14flx), imt, jmt)
#   if defined uvic_embm_c14data
      if (tarea .gt. 0) dc14ccn = tdc14ccn/tarea
#   endif
#  endif
# endif
# if defined uvic_o2

!-----------------------------------------------------------------------
!     set boundary conditions for oxygen
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,io2flx), imt, jmt)
# endif
# if defined uvic_cfc11

!-----------------------------------------------------------------------
!     set boundary conditions for CFC11
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,icfc11flx), imt, jmt)
# endif
# if defined uvic_cfc12

!-----------------------------------------------------------------------
!     set boundary conditions for CFC12
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,icfc12flx), imt, jmt)
# endif

!-----------------------------------------------------------------------
!     calculate CO2 forcing
!-----------------------------------------------------------------------
      call co2forc
# if defined uvic_embm_adv_q || defined uvic_embm_adv_t

!-----------------------------------------------------------------------
!     set flags to calculate new coefficients
!-----------------------------------------------------------------------
      newcoef(1,1) = .true.
      newcoef(2,1) = .true.
      newcoef(1,2) = .true.
      newcoef(2,2) = .true.
# endif
# if defined uvic_embm_cropdata_transient

!-----------------------------------------------------------------------
!     update boundary conditions over vegetation
!-----------------------------------------------------------------------
      call gvsbc
# endif
# if defined time_averages

!-----------------------------------------------------------------------
!     zero time averages if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. timavgperts) call ta_embm_snap (is, ie, js, je, 0)
# endif
# if defined time_step_monitor

!-----------------------------------------------------------------------
!     zero time step integrals if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. tsiperts) call ta_embm_tsi (0)
# endif
#endif

      return
      end
