       subroutine embmout (is, ie, js, je)
!=======================================================================
!     output routine for energy-moisture balance model

!     input:
!       is, ie, js, je = starting and ending indicies for i and j

!     based on code by: A. Fanning and M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "param.h"
# include "calendar.h"
# include "csbc.h"
# include "atm.h"
# include "solve.h"
# include "coord.h"
# include "grdvar.h"
# if defined uvic_ice
#  include "ice.h"
# endif
# if defined uvic_ice_evp
#  include "evp.h"
# endif
# if defined uvic_mtlm
#  include "mtlm.h"
# endif
# include "cembm.h"
# include "iounit.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
# include "riv.h"
# if defined uvic_mom
#  include "cregin.h"
# else
      integer mskhr(imt,jmt)
# endif

# if defined uvic_embm_running_average
      real tmpij(2:imt-1,2:jmt-1), tmpi(2:imt-1), tmpj(2:jmt-1)
# endif
      real tmp_at(imt,jmt,nat)
# if defined uvic_ice_cpts
      real tmp_Ts(imt,jmt,ncat), tmp_heff(imt,jmt,ncat)
      real tmp_A(imt,jmt,ncat), tmp_hseff(imt,jmt,ncat)
# endif
# if defined uvic_mtlm
      real sm(imt,jmt), st(imt,jmt), hs(imt,jmt), ro(imt,jmt)
      real rntatsl, rntatil
# endif

      character(120) :: fname

      integer i, ie, is, id_xt, id_yt, iou, j, je, js, n
      integer ndx, ntrec, it(10), ib(10), ic(10)

      real avgper, ca, time, wt, wtp1
      real p_alb(is:ie,js:je), a_alb(is:ie,js:je), s_alb(is:ie,js:je)
      real sat(is:ie,js:je)
      real c100, c400, C2K, tmp

      c100 = 100.
      c400 = 400.
      C2K = 273.15

# if !defined uvic_mom
      mskhr(:,:) = 0.0
# endif

# if defined uvic_embm_running_average

      if (snapts .or. timavgts .or. eorun) then

!-----------------------------------------------------------------------
!       write running average air temperature and humidity
!-----------------------------------------------------------------------

        call embmbc (rtbar)
        fname = new_file_name ("rtbar.nc")
        call opennew (fname, iou)
        call redef (iou)
        call defdim ('xt', iou, imt-2, id_xt)
        call defdim ('yt', iou, jmt-2, id_yt)
        it(:) = id_xt
        call defvar ('xt', iou, 1, it, c0, c0, 'X', 'F'
     &,   'longitude of the t grid', 'grid_longitude', 'degrees_east')
        it(:) = id_yt
        call defvar ('yt', iou, 1, it, c0, c0, 'Y', 'F'
     &,   'latitude of the t grid', 'grid_latitude', 'degrees_north')
        it(1) = id_xt
        call defvar ('slat', iou, 2, it, c100, c400
     &,   '', 'F', 'running annual average sea level air temperature'
     &,   '', 'K')
        call enddef (iou)
        ib(:) = 1
        ic(:) = imt-2
        tmpi(2:imt-1) = xt(2:imt-1)
        call putvara ('xt', iou, imt-2, ib, ic, tmpi, c1, c0)
        ic(:) = jmt-2
        tmpj(2:jmt-1) = yt(2:jmt-1)
        call putvara ('yt', iou, jmt-2, ib, ic, tmpj, c1, c0)
        ib(:) = 1
        ic(:) = 1
        ic(1) = imt-2
        ic(2) = jmt-2
        tmpij(2:imt-1,2:jmt-2) = rtbar(2:imt-1,2:jmt-2)
        call putvara ('slat', iou, imt*jmt, ib, ic, tmpij, c1, -C2K)
        call closefile (iou)

      endif
# endif

!-----------------------------------------------------------------------
!     write atmospheric diagnostics
!-----------------------------------------------------------------------

# if defined time_step_monitor
      if (tsits .and. ntatia .ne. 0) then

        call ta_embm_tsi (2)

        if (iotsi .eq. stdout .or. iotsi .lt. 0) then
          write (*,'(1x, a3, i7, 1x, a32, 3(a,1pe13.6))')
     &      'ts=',itt, stamp, ' iterations =', tai_maxit
     &,     ' TAbar=', tai_sat, ' QAbar=', tai_shum
        endif

        if (iotsi .ne. stdout) then
          avgper = dtatm*ntatia/daylen
          time = relyr + year0
#  if defined uvic_mtlm
          rntatil = 0.
          if (ntatil .gt. 0) rntatil = 1./float(ntatil)
          tai_hsno = tai_hsno + tai_LYING_SNOW*1000.*rntatil/rhosno
#  endif
          call def_tsi
          call def_tsi_embm (fname)
          call embm_tsi_out (fname, avgper, time, stamp, tai_sat
     &,                      tai_shum, tai_precip, tai_evap, tai_ohice
     &,                      tai_oaice, tai_hsno, tai_lhice, tai_laice
     &,                      tai_co2ccn, tai_dc14ccn, tai_cfc11ccn
     &,                      tai_cfc12ccn, tai_maxit, tai_sst, tai_sss
     &,                      tai_nsat, tai_ssat, tai_nshum, tai_sshum
     &,                      tai_nprecip, tai_sprecip, tai_nevap
     &,                      tai_sevap, tai_nohice, tai_sohice
     &,                      tai_noaice, tai_soaice, tai_nhsno
     &,                      tai_shsno, tai_nlhice, tai_slhice
     &,                      tai_nlaice, tai_slaice, ntrec)
        endif

        call ta_embm_tsi (0)

      endif

# endif
# if defined time_averages
      if (timavgts .and. ntatsa .ne. 0) then

!-----------------------------------------------------------------------
!       write atmospheric time averaged data
!-----------------------------------------------------------------------

!       calculate average values

        call ta_embm_snap (is, ie, js, je, 2)

!       write time averaged data

#  if defined uvic_mtlm
        rntatsl = 0.
        if (ntatsl .gt. 0) rntatsl = 1./float(ntatsl)
        call unloadland (POINTS, TA_M, imt, jmt, land_map, sm)
        call unloadland (POINTS, TA_TSTAR_GB, imt, jmt, land_map, st)
        call unloadland (POINTS, TA_LYING_SNOW, imt, jmt, land_map, hs)
        call unloadland (POINTS, TA_SURF_ROFF, imt, jmt, land_map, ro)
#  endif
        do j=js,je
          do i=is,ie
            sat(i,j) =  ta_at(i,j,isat) - elev(i,j)*rlapse
#  if defined uvic_embm_icedata
     &                - hicel(i,j,2)*rlapse
#  endif
            if (ta_solins(i,j) .gt. 0) then
!             calculate incoming swr reaching the surface
              tmp = ta_solins(i,j) - ta_arswr(i,j) - ta_absin(i,j)
              s_alb(i,j) = 1. - ta_dnswr(i,j)/tmp
              a_alb(i,j) = ta_arswr(i,j)/ta_solins(i,j)
!             calculate outgoing swr leaving the atmosphere
              tmp = tmp - ta_dnswr(i,j) - ta_absout(i,j)
              p_alb(i,j) = (ta_arswr(i,j) + tmp)/ta_solins(i,j)
            else
              s_alb(i,j) = -1.
              a_alb(i,j) = -1.
              p_alb(i,j) = -1.
            endif
#  if defined uvic_mtlm
            if (land_map(i,j) .eq. 0) then
              sm(i,j) = ta_soilm(i,j)
              st(i,j) = ta_surf(i,j)
            else
!             convert from kg m-2 to cm (converted back later)
              sm(i,j) = sm(i,j)*rntatsl*0.1
!             convert from K to C (converted back later)
              st(i,j) = st(i,j)*rntatsl - 273.15
!             convert from kg/m2/s to g/cm2/s (converted back later)
              ta_runoff(i,j) = ro(i,j)*rntatsl*0.1
#   if defined uvic_ice
      !       convert from kg/m2 to cm (converted back later)
              ta_hsno(i,j,1) = hs(i,j)*rntatsl*0.1/rhosno
#   endif
            endif
#  endif
          enddo
        enddo

        avgper = dtatm*ntatsa/daylen
        time = relyr + year0
        call def_tavg
        call def_tavg_embm (fname)
        call embm_snap_out (fname, is, ie, js, je, imt, jmt, nat, ncat
     &,   xt, yt, xu, yu, dxt, dyt, dxu, dyu, avgper, time, stamp, mapat
     &,   ta_at, sat, ta_precip, ta_evap, ta_disch, ta_outlwr, ta_uplwr
     &,   ta_upsens, ta_dnswr, ta_upltnt, ta_solins, p_alb, a_alb, s_alb
     &,   elev, ta_psno, ta_ws, ta_runoff
#  if defined uvic_embm_adv_q
     &,   ta_wx_q, ta_wy_q
#  endif
#  if defined uvic_embm_adv_t
     &,   ta_wx_t, ta_wy_t
#  endif
#  if defined uvic_embm_awind
     &,   ta_awx, ta_awy, ta_rtbar, ta_apress
#  endif
#  if defined uvic_mtlm
     &,   sm, st
#  else
     &,   ta_soilm, ta_surf
#  endif
#  if defined uvic_ice
     &,   ta_tice, ta_hice, ta_aice, ta_hsno
#  endif
#  if defined uvic_ice_evp
     &,   ta_uice, ta_vice, ta_xint, ta_yint
#  endif
     &,   tlat, tlon, ulat, ulon
#  if defined uvic_save_flxadj
     &,   ta_flxadj(is,js,1), ta_flxadj(is,js,2)
#  endif
#  if defined uvic_embm_save_diff
     &,   ta_dn(is,js,1), ta_de(is,js,1)
#  endif
#  if defined uvic_embm_icedata
     &,   ta_aicel(is,js), ta_hicel(is,js)
#  endif
#  if defined uvic_carbon_co2_2d
     &,   ta_flxco2(is,js)
#  endif
     &,   tmsk, mskhr, nriv, ntrec)

        write (*,'(a,i5,a,a,a,a)') '=> Atm time means #'
     &,   ntrec, ' written to ',trim(fname),' on ', stamp

!       zero time average accumulators

        call ta_embm_snap (is, ie, js, je, 0)

      endif

# endif
# if defined snapshots
      if (snapts) then

!-----------------------------------------------------------------------
!       write atmospheric snapshot
!-----------------------------------------------------------------------

#  if defined uvic_ice_cpts
        ndx = 1
        if (idx .eq. 1) ndx = 2
#  endif

!       write snapshot

        tmp_at(is:ie,js:je,1:nat) = at(is:ie,js:je,2,1:nat)
#  if defined uvic_ice_cpts
        tmp_Ts(is:ie,js:je,1:ncat) = Ts(is:ie,js:je,ndx,1:ncat)
        tmp_heff(is:ie,js:je,1:ncat) = heff(is:ie,js:je,ndx,1:ncat)
        tmp_A(is:ie,js:je,1:ncat) = A(is:ie,js:je,ndx,1:ncat)
        tmp_hseff(is:ie,js:je,1:ncat) = hseff(is:ie,js:je,ndx,1:ncat)
#  endif

#  if defined uvic_mtlm
        call unloadland (POINTS, M, imt, jmt, land_map, sm)
        call unloadland (POINTS, TSTAR_GB, imt, jmt, land_map, st)
        call unloadland (POINTS, LYING_SNOW, imt, jmt, land_map, hs)
        call unloadland (POINTS, SURF_ROFF, imt, jmt, land_map, ro)
#  endif
        do j=js,je
          do i=is,ie
            sat(i,j) =  at(i,j,2,isat) - elev(i,j)*rlapse
#  if defined uvic_embm_icedata
     &                - hicel(i,j,2)*rlapse
#  endif
            if (solins(i,j) .gt. 0) then
!             calculate incoming swr reaching the surface
              tmp = solins(i,j)*sbc(i,j,iaca)*pass
              s_alb(i,j) = 1. - dnswr(i,j)/tmp
              a_alb(i,j) = 1. - sbc(i,j,iaca)
!             calculate outgoing swr leaving the atmosphere
              tmp = tmp - dnswr(i,j)
              p_alb(i,j) = a_alb(i,j) + tmp/solins(i,j)
            else
              s_alb(i,j) = -1.
              a_alb(i,j) = -1.
              p_alb(i,j) = -1.
            endif
#  if defined uvic_ice_cpts
            if (tmsk(i,j) .ge. 0.5) then
              Ts(i,j,ndx,1) = tice(i,j)
              heff(i,j,ndx,1) = hice(i,j,2)
              A(i,j,ndx,1) = aice(i,j,2)
              hseff(i,j,ndx,1) = hsno(i,j,2)
            endif
#  endif
#  if defined uvic_mtlm
            if (land_map(i,j) .eq. 0) then
              sm(i,j) = soilm(i,j,2)
              st(i,j) = surf(i,j)
#   if defined uvic_ice && !defined uvic_ice_cpts
              hs(i,j) = hsno(i,j,2)
#   endif
            else
!             convert from kg m-2 to cm (converted back later)
              sm(i,j) = sm(i,j)*0.1
!             convert from K to C (converted back later)
              st(i,j) = st(i,j) - 273.15
!             convert from kg/m2 to cm (converted back later)
#   if defined uvic_ice_cpts
              tmp_hseff(i,j,1) = hs(i,j)*0.1/rhosno
#   else
              hs(i,j) = hs(i,j)*0.1/rhosno
#   endif
            endif
#  endif
          enddo
        enddo

        avgper = 0.
        time = relyr + year0
        call def_snap
        call def_snap_embm (fname)
        call embm_snap_out (fname, is, ie, js, je, imt, jmt, nat, ncat
     &,   xt, yt, xu, yu, dxt, dyt, dxu, dyu, avgper, time, stamp, mapat
     &,   tmp_at, sat, precip, evap, disch, outlwr, uplwr, upsens, dnswr
     &,   upltnt, solins, p_alb, a_alb, s_alb, elev, psno
     &,   sbc(is,js,iws), sbc(is,js,iro)
#  if defined uvic_embm_adv_q
     &,   sbc(is,js,iwxq), sbc(is,js,iwyq)
#  endif
#  if defined uvic_embm_adv_t
     &,   sbc(is,js,iwxt), sbc(is,js,iwyt)
#  endif
#  if defined uvic_embm_awind
     &,   awx, awy, rtbar, apress
#  endif
#  if defined uvic_mtlm
     &,   sm, st
#  else
     &,   soilm(is,js,2), surf
#  endif
#  if defined uvic_ice_cpts
     &,   tmp_Ts, tmp_heff, tmp_A, tmp_hseff
#  elif defined uvic_ice
#   if defined uvic_mtlm
     &,   tice(is,js), hice(is,js,2), aice(is,js,2), hs(is,js)
#   else
     &,   tice(is,js), hice(is,js,2), aice(is,js,2), hsno(is,js,2)
#   endif
#  endif
#  if defined uvic_ice_evp
     &,   uice, vice, xint, yint
#  endif
     &,   tlat, tlon, ulat, ulon
#  if defined uvic_save_flxadj
     &,   flxadj(is,js,1), flxadj(is,js,2)
#  endif
#  if defined uvic_embm_save_diff
     &,   dn(is,js,1), de(is,js,1)
#  endif
#  if defined uvic_embm_icedata
     &,   aicel(is,js,2), hicel(is,js,2)
#  endif
#  if defined uvic_carbon_co2_2d
     &,   flux(is,js,ico2)
#  endif
     &,   tmsk, mskhr, nriv, ntrec)

        write (*,'(a,i5,a,a,a,a)') '=> Atm snapshot #'
     &,   ntrec, ' written to ',trim(fname),' on ', stamp

#  if defined uvic_ice_cpts

        do j=js,je
          do i=is,ie
            if (tmsk(i,j) .ge. 0.5) then
              Ts(i,j,ndx,1) = 0.
              heff(i,j,ndx,1) = 0.
              A(i,j,ndx,1) = 0.
              hseff(i,j,ndx,1) = 0.
            endif
          enddo
        enddo

#  endif
      endif

# endif
!-----------------------------------------------------------------------
!       write atmospheric restart
!-----------------------------------------------------------------------

      if (eoyear) call eoyr (1, imt, 1, jmt)

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_embm (0, fname)
          call embm_rest_out (fname, is, ie, js, je)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_embm (1, fname)
          call embm_rest_out (fname, is, ie, js, je)
        endif
      endif
#endif

      return
      end

#if defined uvic_embm && defined time_averages
      subroutine ta_embm_snap (is, ie, js, je, m)

!=======================================================================
!     atmospheric data time averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       m = flag (0 = zero, 1 = accumulate, 2 = write)

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "param.h"
# include "cembm.h"
# include "atm.h"
# include "riv.h"
# if defined uvic_ice
#  include "ice.h"
# endif
# include "csbc.h"

      integer i, ie, is, j, je, js, k, m, n, ndx

      real rntatsa

!-----------------------------------------------------------------------
!     time averaged data
!-----------------------------------------------------------------------

      if (m .eq. 0.) then
!       zero
        ntatsa = 0
        ta_at(:,:,:) = 0.
        ta_solins(:,:) = 0.
        ta_arswr(:,:) = 0.
        ta_dnswr(:,:) = 0.
        ta_absin(:,:) = 0.
        ta_absout(:,:) = 0.
        ta_uplwr(:,:) = 0.
        ta_upsens(:,:) = 0.
        ta_upltnt(:,:) = 0.
        ta_outlwr(:,:) = 0.
        ta_precip(:,:) = 0.
        ta_psno(:,:) = 0.
        ta_evap(:,:) = 0.
        ta_disch(:,:) = 0.
        ta_ws(:,:) = 0.
        ta_runoff(:,:) = 0.
# if defined uvic_embm_adv_q
        ta_wx_q(:,:) = 0.
        ta_wy_q(:,:) = 0.
# endif
# if defined uvic_embm_adv_t
        ta_wx_t(:,:) = 0.
        ta_wy_t(:,:) = 0.
# endif
# if defined uvic_embm_awind
        ta_awx(:,:) = 0.
        ta_awy(:,:) = 0.
        ta_rtbar(:,:) = 0.
        ta_apress(:,:) = 0.
# endif
        ta_soilm(:,:) = 0.
        ta_surf(:,:) = 0.
# if defined uvic_ice
        ta_tice(:,:,:) = 0.
        ta_hice(:,:,:) = 0.
        ta_aice(:,:,:) = 0.
        ta_hsno(:,:,:) = 0.
#  if defined uvic_ice_evp
        ta_uice(:,:) = 0.
        ta_vice(:,:) = 0.
        ta_pice(:,:) = 0.
        ta_xint(:,:) = 0.
        ta_yint(:,:) = 0.
#  endif
# endif
# if defined uvic_save_flxadj
        ta_flxadj(:,:,:) = 0.
# endif
# if defined uvic_embm_save_diff
        ta_dn(:,:,:) = 0.
        ta_de(:,:,:) = 0.
# endif
# if defined uvic_embm_icedata
        ta_aicel(:,:) = 0.
        ta_hicel(:,:) = 0.
# endif
# if defined uvic_carbon && defined uvic_carbon_co2_2d
        ta_flxco2(:,:) = 0.
# endif
# if !defined uvic_mom
        ta_sst(:,:) = 0.
        ta_sss(:,:) = 0.
        ta_hflux(:,:) = 0.
        ta_sflux(:,:) = 0.
        ta_taux(:,:) = 0.
        ta_tauy(:,:) = 0.
# endif
        ta_psum(:) = 0.
      elseif (m .eq. 1) then
!       accumulate
        ntatsa = ntatsa + 1
# if defined uvic_ice_cpts
        ndx = 1
        if (idx .eq. 1) ndx = 2
# endif
        do n=1,nat
          ta_at(:,:,n) = ta_at(:,:,n) + at(:,:,2,n)
        enddo
        ta_solins(:,:) = ta_solins(:,:) + solins(:,:)
        ta_arswr(:,:) = ta_arswr(:,:) + solins(:,:)
     &                  *(1. - sbc(:,:,iaca))
        ta_dnswr(:,:) = ta_dnswr(:,:) + dnswr(:,:)
        ta_absin(:,:) = ta_absin(:,:) + solins(:,:)*sbc(:,:,iaca)
     &                  *scatter
        ta_absout(:,:) = ta_absout(:,:) + (solins(:,:)*sbc(:,:,iaca)
     &                   *pass - dnswr(:,:))*scatter
        ta_uplwr(:,:) = ta_uplwr(:,:) + uplwr(:,:)
        ta_upsens(:,:) = ta_upsens(:,:) + upsens(:,:)
        ta_upltnt(:,:) = ta_upltnt(:,:) + upltnt(:,:)
        ta_outlwr(:,:) = ta_outlwr(:,:) + outlwr(:,:)
        ta_precip(:,:) = ta_precip(:,:) + precip(:,:)
        ta_psno(:,:) = ta_psno(:,:) + psno(:,:)
        ta_evap(:,:) = ta_evap(:,:) + evap(:,:)
        ta_disch(:,:) = ta_disch(:,:) + disch(:,:)
        ta_ws(:,:) = ta_ws(:,:) + sbc(:,:,iws)
        ta_runoff(:,:) = ta_runoff(:,:) + runoff(:,:)
# if defined uvic_embm_adv_q
        ta_wx_q(:,:) = ta_wx_q(:,:) + sbc(:,:,iwxq)
        ta_wy_q(:,:) = ta_wy_q(:,:) + sbc(:,:,iwyq)
# endif
# if defined uvic_embm_adv_t
        ta_wx_t(:,:) = ta_wx_t(:,:) + sbc(:,:,iwxt)
        ta_wy_t(:,:) = ta_wy_t(:,:) + sbc(:,:,iwyt)
# endif
# if defined uvic_embm_awind
        ta_awx(:,:) = ta_awx(:,:) + awx(:,:)
        ta_awy(:,:) = ta_awy(:,:) + awy(:,:)
        ta_rtbar(:,:) = ta_rtbar(:,:) + rtbar(:,:)
        ta_apress(:,:) = ta_apress(:,:) + apress(:,:)
# endif
        ta_soilm(:,:) = ta_soilm(:,:) + soilm(:,:,2)
        ta_surf(:,:) = ta_surf(:,:) + surf(:,:)
# if defined uvic_ice
#  if defined uvic_ice_cpts
        do k=1,ncat
          ta_hsno(:,:,k) = ta_hsno(:,:,k) + hseff(:,:,ndx,k)
          ta_tice(:,:,k) = ta_tice(:,:,k) + Ts(:,:,ndx,k)
          ta_hice(:,:,k) = ta_hice(:,:,k) + heff(:,:,ndx,k)
          ta_aice(:,:,k) = ta_aice(:,:,k) + A(:,:,ndx,k)
        enddo
#  else
        ta_tice(:,:,1) = ta_tice(:,:,1) + tice(:,:)
        ta_hice(:,:,1) = ta_hice(:,:,1) + hice(:,:,2)
        ta_aice(:,:,1) = ta_aice(:,:,1) + aice(:,:,2)
        ta_hsno(:,:,1) = ta_hsno(:,:,1) + hsno(:,:,2)
#  endif
#  if defined uvic_ice_evp
        ta_uice(:,:) = ta_uice(:,:) + uice(:,:)
        ta_vice(:,:) = ta_vice(:,:) + vice(:,:)
        ta_pice(:,:) = ta_pice(:,:) + pice(:,:)
        ta_xint(:,:) = ta_xint(:,:) + xint(:,:)
        ta_yint(:,:) = ta_yint(:,:) + yint(:,:)
#  endif
# endif
# if defined uvic_save_flxadj
        ta_flxadj(:,:,:) = ta_flxadj(:,:,:) + flxadj(:,:,:)
# endif
# if defined uvic_embm_save_diff
        ta_dn(:,:,:) = ta_dn(:,:,:) + dn(:,:,:)
        ta_de(:,:,:) = ta_de(:,:,:) + de(:,:,:)
# endif
# if defined uvic_embm_icedata
        ta_aicel(:,:) = ta_aicel(:,:) + aicel(:,:,2)
        ta_hicel(:,:) = ta_hicel(:,:) + hicel(:,:,2)
# endif
# if defined uvic_carbon && defined uvic_carbon_co2_2d
        ta_flxco2(:,:) = ta_flxco2(:,:) + flux(:,:,ico2)
# endif
# if !defined uvic_mom
        ta_sst(:,:) = ta_sst(:,:) + sbc(:,:,isst)
        ta_sss(:,:) = ta_sss(:,:) + sbc(:,:,isss)
        ta_hflux(:,:) = ta_hflux(:,:) + sbc(:,:,ihflx)
        ta_sflux(:,:) = ta_sflux(:,:) + sbc(:,:,isflx)
        ta_taux(:,:) = ta_taux(:,:) + sbc(:,:,itaux)
        ta_tauy(:,:) = ta_tauy(:,:) + sbc(:,:,itauy)
# endif
        ta_psum(:) = ta_psum(:) + psum(:)
      elseif (m .eq. 2 .and. ntatsa .ne. 0) then
!       average
        rntatsa = 1./float(ntatsa)
        ta_at(:,:,:) = ta_at(:,:,:)*rntatsa
        ta_solins(:,:) = ta_solins(:,:)*rntatsa
        ta_arswr(:,:) = ta_arswr(:,:)*rntatsa
        ta_dnswr(:,:) = ta_dnswr(:,:)*rntatsa
        ta_absin(:,:) = ta_absin(:,:)*rntatsa
        ta_absout(:,:) = ta_absout(:,:)*rntatsa
        ta_uplwr(:,:) = ta_uplwr(:,:)*rntatsa
        ta_upsens(:,:) = ta_upsens(:,:)*rntatsa
        ta_upltnt(:,:) = ta_upltnt(:,:)*rntatsa
        ta_outlwr(:,:) = ta_outlwr(:,:)*rntatsa
        ta_precip(:,:) = ta_precip(:,:)*rntatsa
        ta_psno(:,:) = ta_psno(:,:)*rntatsa
        ta_evap(:,:) = ta_evap(:,:)*rntatsa
        ta_disch(:,:) = ta_disch(:,:)*rntatsa
        ta_ws(:,:) = ta_ws(:,:)*rntatsa
        ta_runoff(:,:) = ta_runoff(:,:)*rntatsa
# if defined uvic_embm_adv_q
        ta_wx_q(:,:) = ta_wx_q(:,:)*rntatsa
        ta_wy_q(:,:) = ta_wy_q(:,:)*rntatsa
# endif
# if defined uvic_embm_adv_t
        ta_wx_t(:,:) = ta_wx_t(:,:)*rntatsa
        ta_wy_t(:,:) = ta_wy_t(:,:)*rntatsa
# endif
# if defined uvic_embm_awind
        ta_awx(:,:) = ta_awx(:,:)*rntatsa
        ta_awy(:,:) = ta_awy(:,:)*rntatsa
        ta_rtbar(:,:) = ta_rtbar(:,:)*rntatsa
        ta_apress(:,:) = ta_apress(:,:)*rntatsa
# endif
        ta_soilm(:,:) = ta_soilm(:,:)*rntatsa
        ta_surf(:,:) = ta_surf(:,:)*rntatsa
# if defined uvic_ice
        ta_hsno(:,:,:) = ta_hsno(:,:,:)*rntatsa
        ta_tice(:,:,:) = ta_tice(:,:,:)*rntatsa
        ta_hice(:,:,:) = ta_hice(:,:,:)*rntatsa
        ta_aice(:,:,:) = ta_aice(:,:,:)*rntatsa
#  if defined uvic_ice_evp
        ta_uice(:,:) = ta_uice(:,:)*rntatsa
        ta_vice(:,:) = ta_vice(:,:)*rntatsa
        ta_pice(:,:) = ta_pice(:,:)*rntatsa
        ta_xint(:,:) = ta_xint(:,:)*rntatsa
        ta_yint(:,:) = ta_yint(:,:)*rntatsa
#  endif
# endif
# if defined uvic_save_flxadj
        ta_flxadj(:,:,1) = ta_flxadj(:,:,1)*rntatsa
        ta_flxadj(:,:,2) = ta_flxadj(:,:,2)*rntatsa
# endif
# if defined uvic_embm_save_diff
        ta_dn(:,:,:) = ta_dn(:,:,:)*rntatsa
        ta_de(:,:,:) = ta_de(:,:,:)*rntatsa
# endif
# if defined uvic_embm_icedata
        ta_aicel(:,:) = ta_aicel(:,:)*rntatsa
        ta_hicel(:,:) = ta_hicel(:,:)*rntatsa
# endif
# if defined uvic_carbon && defined uvic_carbon_co2_2d
        ta_flxco2(:,:) = ta_flxco2(:,:)*rntatsa
# endif
# if !defined uvic_mom
        ta_sst(:,:) = ta_sst(:,:)*rntatsa
        ta_sss(:,:) = ta_sss(:,:)*rntatsa
        ta_hflux(:,:) = ta_hflux(:,:)*rntatsa
        ta_sflux(:,:) = ta_sflux(:,:)*rntatsa
        ta_taux(:,:) = ta_taux(:,:)*rntatsa
        ta_tauy(:,:) = ta_tauy(:,:)*rntatsa
# endif
        ta_psum(:) = ta_psum(:)*rntatsa
      endif

      return
      end
#endif

#if defined uvic_embm && defined time_step_monitor
      subroutine ta_embm_tsi (m)

!=======================================================================
!     atmospheric data time integral averaging

!     input:
!       m = flag (0 = zero, 1 = accumulate, 2 = write)

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "param.h"
# include "csbc.h"
# include "grdvar.h"
# include "cembm.h"
# include "atm.h"
# if defined uvic_ice
# include "ice.h"
# endif
# include "solve.h"

      integer maxit, m, n
      real rntatia, sum, tmp, sat(imt,jmt), dmsk(imt,jmt)

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

       if (m .eq. 0.) then
!       zero
        ntatia = 0
        tai_maxit = 0.
        tai_sat  = 0.
        tai_shum = 0.
        tai_co2ccn = 0.
        tai_dc14ccn = 0.
        tai_cfc11ccn = 0.
        tai_cfc12ccn = 0.
        tai_precip = 0.
        tai_evap = 0.
        tai_ohice = 0.
        tai_oaice = 0.
        tai_hsno = 0.
        tai_lhice = 0.
        tai_laice = 0.
        tai_sst = 0.
        tai_sss = 0.
        tai_nsat = 0.
        tai_ssat = 0.
        tai_nshum = 0.
        tai_sshum = 0.
        tai_nprecip = 0.
        tai_sprecip = 0.
        tai_nevap = 0.
        tai_sevap = 0.
        tai_nohice = 0.
        tai_sohice = 0.
        tai_noaice = 0.
        tai_soaice = 0.
        tai_nhsno = 0.
        tai_shsno = 0.
        tai_nlhice = 0.
        tai_slhice = 0.
        tai_nlaice = 0.
        tai_slaice = 0.
      elseif (m .eq. 1) then
!       accumulate
        ntatia = ntatia + 1
        sat(:,:) = at(:,:,2,isat) - elev(:,:)*rlapse
# if defined uvic_embm_icedata
     &            - hicel(:,:,2)*rlapse
# endif
        maxit = 0
!       set data mask
        dmsk(:,:) = 1.
# if defined uvic_embm_explicit
        maxit = ns
# else
        do n=1,nat
          if (itout(n) .gt. maxit) maxit = itout(n)
        enddo
# endif
        tai_maxit = tai_maxit + maxit
        call areaavg (sat, dmsk, tmp)
        tai_sat  = tai_sat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_shum = tai_shum + tmp
# if defined uvic_carbon_co2_2d
        call areaavg (at(1,1,2,ico2), dmsk, tmp)
        tai_co2ccn = tai_co2ccn + tmp
# else
        tai_co2ccn = tai_co2ccn + co2ccn
# endif
# if defined uvic_carbon && defined uvic_carbon_14
        tai_dc14ccn = tai_dc14ccn + dc14ccn
# endif
# if defined uvic_cfc11
        tai_cfc11ccn = tai_cfc11ccn + 0.5*(cfc11ccnn + cfc11ccns)
# endif
# if defined uvic_cfc11
        tai_cfc12ccn = tai_cfc12ccn + 0.5*(cfc12ccnn + cfc12ccns)
# endif
        call areaavg (precip, dmsk, tmp)
        tai_precip = tai_precip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_evap = tai_evap + tmp
# if defined uvic_ice
#  if defined uvic_ice_cpts
        sum = 0.
        do n=1,ncat
           call areatot (hseff(:,:,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_hsno = tai_hsno + sum
!       mask for ocean
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        sum = 0.
        do n=1,ncat
           call areatot (heff(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_ohice = tai_ohice + sum
        sum = 0.
        do n=1,ncat
           call areatot (A(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_oaice = tai_oaice + sum
#  else
        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_hsno = tai_hsno + tmp
!       mask for ocean
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_ohice = tai_ohice + tmp
        call areatot (aice(1,1,2), dmsk, tmp)
        tai_oaice = tai_oaice + tmp
#  endif
#  if defined uvic_embm_icedata
!       mask for land
        dmsk(:,:) = 1.
        where (tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.
        call areatot (hicel(1,1,2), dmsk, tmp)
        tai_lhice = tai_lhice + tmp
        call areatot (aicel(1,1,2), dmsk, tmp)
        tai_laice = tai_laice + tmp
#  endif
# endif
!       mask for ocean
        dmsk(:,:) = 1.
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areaavg (sbc(1,1,isst), dmsk, tmp)
        tai_sst = tai_sst + tmp
        call areaavg (sbc(1,1,isss), dmsk, tmp)
        tai_sss = tai_sss + tmp
# if defined uvic_tai_ns
!       northern hemisphere mask
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0.) dmsk(:,:) = 0.
        call areaavg (sat, dmsk, tmp)
        tai_nsat  = tai_nsat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_nshum  = tai_nshum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_nprecip  = tai_nprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_nevap  = tai_nevap + tmp
#  if defined uvic_ice
#   if defined uvic_ice_cpts
        sum = 0.
        do n=1,ncat
           call areatot (hseff(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_nhsno = tai_nhsno + sum
!       mask for northern hemisphere ocean
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        sum = 0.
        do n=1,ncat
           call areatot (heff(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_nohice = tai_nohice + sum
        sum = 0.
        do n=1,ncat
           call areatot (A(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_noaice = tai_noaice + sum
#   else
        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_nhsno = tai_nhsno + tmp
!       mask for northern hemisphere ocean
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_nohice = tai_nohice + tmp
        call areatot (aice(1,1,2), dmsk, tmp)
        tai_noaice = tai_noaice + tmp

#   endif
#   if defined uvic_embm_icedata
!       mask for northern hemisphere land
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0. .or. tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.
        call areatot (hicel(1,1,2), dmsk, tmp)
        tai_nlhice = tai_nlhice + tmp
        call areatot (aicel(1,1,2), dmsk, tmp)
        tai_nlaice = tai_nlaice + tmp
#   endif
#  endif

!       southern hemisphere mask
        dmsk(:,:) = 1.
        where (tlat(:,:) .gt. 0.) dmsk(:,:) = 0.
        call areaavg (sat, dmsk, tmp)
        tai_ssat  = tai_ssat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_sshum  = tai_sshum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_sprecip  = tai_sprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_sevap  = tai_sevap + tmp
#  if defined uvic_ice
#   if defined uvic_ice_cpts
        sum = 0.
        do n=1,ncat
           call areatot (hseff(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_shsno = tai_shsno + sum
!       mask for southern hemisphere ocean
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        sum = 0.
        do n=1,ncat
           call areatot (heff(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_sohice = tai_sohice + sum
        sum = 0.
        do n=1,ncat
           call areatot (A(1,1,2,n), dmsk, tmp)
           sum = sum + tmp
        enddo
        tai_soaice = tai_soaice + sum
#   else
        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_shsno = tai_shsno + tmp
!       mask for southern hemisphere ocean
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_sohice = tai_sohice + tmp
        call areatot (aice(1,1,2), dmsk, tmp)
        tai_soaice = tai_soaice + tmp
#   endif
#   if defined uvic_embm_icedata
!       mask for southern hemisphere land
        dmsk(:,:) = 1.
        where (tlat(:,:) .gt. 0. .or. tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.
        call areatot (hicel(1,1,2), dmsk, tmp)
        tai_slhice = tai_slhice + tmp
        call areatot (aicel(1,1,2), dmsk, tmp)
        tai_slaice = tai_slaice + tmp
#   endif
#  endif
# endif
      elseif (m .eq. 2 .and. ntatia .ne. 0) then
!       average
        rntatia = 0.
        if (ntatia .gt. 0.) rntatia = 1./float(ntatia)
        tai_maxit = tai_maxit*rntatia
        tai_sat  = tai_sat*rntatia
        tai_shum = tai_shum*rntatia
        tai_co2ccn = tai_co2ccn*rntatia
        tai_dc14ccn = tai_dc14ccn*rntatia
        tai_cfc11ccn = tai_cfc11ccn*rntatia
        tai_cfc12ccn = tai_cfc12ccn*rntatia
        tai_precip = tai_precip*rntatia
        tai_evap = tai_evap*rntatia
        tai_ohice = tai_ohice*rntatia
        tai_oaice = tai_oaice*rntatia
        tai_hsno = tai_hsno*rntatia
        tai_lhice = tai_lhice*rntatia
        tai_laice = tai_laice*rntatia
        tai_sst = tai_sst*rntatia
        tai_sss = tai_sss*rntatia
        tai_nsat = tai_nsat*rntatia
        tai_ssat = tai_ssat*rntatia
        tai_nshum = tai_nshum*rntatia
        tai_sshum = tai_sshum*rntatia
        tai_nprecip = tai_nprecip*rntatia
        tai_sprecip = tai_sprecip*rntatia
        tai_nevap = tai_nevap*rntatia
        tai_sevap = tai_sevap*rntatia
        tai_nohice = tai_nohice*rntatia
        tai_sohice = tai_sohice*rntatia
        tai_noaice = tai_noaice*rntatia
        tai_soaice = tai_soaice*rntatia
        tai_nhsno = tai_nhsno*rntatia
        tai_shsno = tai_shsno*rntatia
        tai_nlhice = tai_nlhice*rntatia
        tai_slhice = tai_slhice*rntatia
        tai_nlaice = tai_nlaice*rntatia
        tai_slaice = tai_slaice*rntatia
      endif

      return
      end
#endif
