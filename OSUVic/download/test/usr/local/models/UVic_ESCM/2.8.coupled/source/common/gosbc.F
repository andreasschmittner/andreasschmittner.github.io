      subroutine gosbc

!=======================================================================
!     calculate the average fluxes for next ocean time step

!     based on code by: A. Fanning and M. Eby
!=======================================================================

      implicit none

#if defined uvic_mom || defined uvic_embm
# include "param.h"
# include "csbc.h"
# include "grdvar.h"
# if defined uvic_embm
#  include "cembm.h"
#  include "atm.h"
# endif
# if defined uvic_ice
#  include "ice.h"
# endif
# if defined uvic_mtlm
#  include "mtlm.h"
# endif
# if defined uvic_fwa
#  include "fwa.h"
#  include "levind.h"
#  include "cregin.h"
#  include "tmngr.h"
# endif

      integer i, j, nc

      real f1, f1a, f1l, fh, fs, fwcflx, fwaflx, fzoop, time
      real area, tarea, tsflx, rsocn, tmp

# if defined uvic_mom && defined uvic_embm
      f1 = 1./atatm
      fh = 2.389e-8/atatm
      fs = -socn/atatm

!-----------------------------------------------------------------------
!     calculate average net fluxes. convert heat flux to cal/cm**2/s
!     from kW/m**2 and fresh water flux (cm/s) to an apparent salt
!     flux (g/cm**2/s) using global ocean average salinity, socn
!-----------------------------------------------------------------------
      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .ge. 0.5) then
#  if defined llnl_plume
#   if defined llnl_plume_heat
            subflux(i,j,1) = flux(i,j,isat)
#   endif
#   if defined llnl_plume_salt
            subflux(i,j,2) = flux(i,j,ishum)
#   endif
            if (subflux(i,j,1) .gt. 0.) subflux(i,j,1) = 0.
            if (subflux(i,j,2) .gt. 0.) subflux(i,j,2) = 0.
#   if defined llnl_plume_all_heat
            subflux(i,j,1) = flux(i,j,isat)
#   endif
#   if defined llnl_plume_all_salt
            subflux(i,j,2) = flux(i,j,ishum)
#   endif
            flux(i,j,isat) = flux(i,j,isat) - subflux(i,j,1)
            subflux(i,j,1) = fh*subflux(i,j,1)
            flux(i,j,ishum) = flux(i,j,ishum) - subflux(i,j,2)
            subflux(i,j,2) = fs*subflux(i,j,2)
#  endif
#  if defined uvic_convect_brine
            cba0(i,j) = 0.
            do nc=0,ncat
              cba(i,j,nc) = f1*cba(i,j,nc)
              if (cba(i,j,nc) .gt. 0.) then
                cbf(i,j,nc) = fs*cbf(i,j,nc)/cba(i,j,nc)
                cba0(i,j) = cba0(i,j) + cba(i,j,nc)
              else
                cbf(i,j,nc) = 0.
                cba(i,j,nc) = 0.
              endif
            enddo
            if (cba0(i,j) .gt. 1.) then
              if (cba0(i,j) .gt. 1.000001) then
                print*, "==> Warning: ice area > 1: ", cba0(i,j)
                endif
              cba0(i,j) = 1.
            endif
            cba0(i,j) = 1. - cba0(i,j)
#  endif
            sbc(i,j,ihflx) = sbc(i,j,ihflx) + fh*flux(i,j,isat)
!           add virtual fluxes of salinity
            sbc(i,j,isflx) = sbc(i,j,isflx) + fs*flux(i,j,ishum)

          else
#  if defined llnl_plume
            subflux(i,j,1) = 0.
            subflux(i,j,2) = 0.
#  endif
#  if defined uvic_convect_brine
            cbf(i,j,:) = 0.
            cba(i,j,:) = 0.
            cba0(i,j) = 1.
#  endif
            sbc(i,j,ihflx) = 0.
            sbc(i,j,isflx) = 0.
          endif
          if (umsk(i,j) .ge. 0.5) then
#  if defined uvic_ice_evp || defined uvic_embm_awind
            sbc(i,j,itaux) = f1*flux(i,j,nat+1)
            sbc(i,j,itauy) = f1*flux(i,j,nat+2)
#  endif
          else
            sbc(i,j,itaux) = 0.
            sbc(i,j,itauy) = 0.
          endif
        enddo
      enddo

      call setbcx (sbc(1,1,ihflx), imt, jmt)
      call setbcx (sbc(1,1,isflx), imt, jmt)
      call setbcx (sbc(1,1,itaux), imt, jmt)
      call setbcx (sbc(1,1,itauy), imt, jmt)
#  if defined uvic_convect_brine
      do nc=0,ncat
        call setbcx (cbf(1,1,nc), imt, jmt)
        call setbcx (cba(1,1,nc), imt, jmt)
      enddo
      call setbcx (cba0, imt, jmt)
#  endif
# endif

!-----------------------------------------------------------------------
!     update boundary conditions from the land model
!     do this now instead of in gasbc so fields can be written out
!-----------------------------------------------------------------------

      f1l = 0.
      f1a = 0.
# if defined uvic_embm
      if (atatm .ne. 0.) f1a = 1.0/atatm
# endif
# if defined uvic_mtlm
      if (atlnd .ne. 0.) f1l = 1.0/atlnd
      do j=2,jmtm1
        do i=2,imtm1
          if (land_map(i,j) .ne. 0) then
            sbc(i,j,iro) = sbc(i,j,iro)*f1l
            sbc(i,j,isca) = sbc(i,j,isca)*f1l
            sbc(i,j,ievap) = sbc(i,j,ievap)*f1l
            sbc(i,j,ilwr) = sbc(i,j,ilwr)*f1l
            sbc(i,j,isens) = sbc(i,j,isens)*f1l
#  if defined uvic_carbon
            sbc(i,j,inpp) =  sbc(i,j,inpp)*f1l
            sbc(i,j,isr) =  sbc(i,j,isr)*f1l
#  endif
          else
            sbc(i,j,iro) = sbc(i,j,iro)*f1a
            sbc(i,j,ievap) = 0.
            sbc(i,j,ilwr) = 0.
            sbc(i,j,isens) = 0.
#  if defined uvic_mtlm && defined uvic_carbon
            sbc(i,j,inpp) = 0.
            sbc(i,j,isr) = 0.
#  endif
          endif
        enddo
      enddo
      call setbcx (sbc(1,1,isca), imt, jmt)
      call setbcx (sbc(1,1,ievap), imt, jmt)
      call setbcx (sbc(1,1,ilwr), imt, jmt)
      call setbcx (sbc(1,1,isens), imt, jmt)
#  if defined uvic_carbon
      call setbcx (sbc(1,1,inpp), imt, jmt)
      call setbcx (sbc(1,1,isr), imt, jmt)
#  endif
# else
      sbc(:,:,iro) = sbc(:,:,iro)*f1a
# endif
      call setbcx (sbc(1,1,iro), imt, jmt)

# if defined uvic_embm
!-----------------------------------------------------------------------
!     zero diagnostic for river discharge and call river model
!-----------------------------------------------------------------------
      disch(:,:) = 0.
      call rivmodel

# endif
# if defined uvic_mom && defined uvic_embm
#  if defined uvic_fwa
!-----------------------------------------------------------------------
!     add additional fresh water flux anomaly
!-----------------------------------------------------------------------

      time = relyr + year0
      if (time .ge. fwayri .and. time .le. fwayrf) then

        fwcflx = 0.
        if (areafwa .gt. 0) then
!         flux is in Sv (1e6 m3 s-1)
          fwaflx = fwaflxi + (time - fwayri)*fwarate
!         convert to flux in g salt cm-2 s-1
          fwaflx = -socn*fwaflx*1.e12/areafwa
        endif

        if (compensate) then

          fwcflx = 0.
          if (areafwc .gt. 0) fwcflx = -fwaflx*areafwa/areafwc
          if (mrfwa .gt. 0 .and. mrfwa .le. nhreg) then
            do j=2,jmtm1
              if (j .ge. jsfwa .and. j .le. jefwa) then
                do i=2,imtm1
                  if (kmt(i,j) .gt. 0) then
                    if (mskhr(i,j) .eq. mrfwa) then
                      sbc(i,j,isflx) = sbc(i,j,isflx) + fwaflx
                    else
                      sbc(i,j,isflx) = sbc(i,j,isflx) + fwcflx
                    endif
                  endif
                enddo
              else
                do i=2,imtm1
                  if (kmt(i,j) .gt. 0) then
                    sbc(i,j,isflx) = sbc(i,j,isflx) + fwcflx
                  endif
                enddo
              endif
            enddo
          else
            do j=2,jmtm1
              if (j .ge. jsfwa .and. j .le. jefwa) then
                do i=2,imtm1
                  if (kmt(i,j) .gt. 0) then
                    if ((i .ge. isfwa1 .and. i .le. iefwa1) .or.
     &                (i .ge. isfwa2 .and. i .le. iefwa2)) then
                      sbc(i,j,isflx) = sbc(i,j,isflx) + fwaflx
                    else
                      sbc(i,j,isflx) = sbc(i,j,isflx) + fwcflx
                    endif
                  endif
                enddo
              else
                do i=2,imtm1
                  if (kmt(i,j) .gt. 0) then
                    sbc(i,j,isflx) = sbc(i,j,isflx) + fwcflx
                  endif
                enddo
              endif
            enddo
          endif

        else

          if (mrfwa .gt. 0 .and. mrfwa .le. nhreg) then
            do j=jsfwa,jefwa
              do i=isfwa1,iefwa1
                if (mskhr(i,j) .eq. mrfwa) then
                  sbc(i,j,isflx) = sbc(i,j,isflx) + fwaflx
                endif
              enddo
            enddo
          else
            do j=jsfwa,jefwa
              do i=isfwa1,iefwa1
                if (kmt(i,j) .gt. 0) then
                  sbc(i,j,isflx) = sbc(i,j,isflx) + fwaflx
                endif
              enddo
              do i=isfwa2,iefwa2
                if (kmt(i,j) .gt. 0) then
                  sbc(i,j,isflx) = sbc(i,j,isflx) + fwaflx
                endif
              enddo
            enddo
          endif

        endif

        call setbcx (sbc(1,1,isflx), imt, jmt)

      endif

#  endif
#  if defined uvic_carbon || defined uvic_alk || defined uvic_o2 || defined uvic_npzd || defined uvic_cfc11 || defined uvic_cfc12
!-----------------------------------------------------------------------
!     add normalized virtual fluxes to other tracers
!-----------------------------------------------------------------------
      tarea = 0.
      tsflx = 0.
      rsocn = 1./socn
      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .ge. 0.5) then
            area = dxt(i)*dyt(j)*cst(j)
            tarea = tarea + area
            tsflx = tsflx + sbc(i,j,isflx)*area
          endif
        enddo
      enddo
      tsflx = tsflx/tarea
      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .ge. 0.5) then
            tmp = (sbc(i,j,isflx) - tsflx)*rsocn
#   if defined uvic_carbon
            sbc(i,j,idicflx) = sbc(i,j,idicflx) + dicocn*tmp
#    if defined uvic_carbon_14
            sbc(i,j,ic14flx) = sbc(i,j,ic14flx) + c14ocn*tmp
#    endif
#   endif
#   if defined uvic_alk
            sbc(i,j,ialkflx) = sbc(i,j,ialkflx) + alkocn*tmp
#   endif
#   if defined uvic_o2
            sbc(i,j,io2flx) = sbc(i,j,io2flx) + o2ocn*tmp
#   endif
#   if defined uvic_npzd
            sbc(i,j,ipo4flx) = sbc(i,j,ipo4flx) + po4ocn*tmp
#    if defined uvic_npzd_vflux
            sbc(i,j,iphytflx) = sbc(i,j,iphytflx) + phytocn*tmp
            sbc(i,j,izoopflx) = sbc(i,j,izoopflx) + zoopocn*tmp
            sbc(i,j,idetrflx) = sbc(i,j,idetrflx) + detrocn*tmp
#    endif
#    if defined uvic_nitrogen
            sbc(i,j,ino3flx) = sbc(i,j,ino3flx) + no3ocn*tmp
#     if defined uvic_npzd_vflux
            sbc(i,j,idiazflx) = sbc(i,j,idiazflx) + diazocn*tmp
#     endif
#    endif
#   endif
#   if defined uvic_cfc11
            sbc(i,j,icfc11flx) = sbc(i,j,icfc11flx) + cfc11ocn*tmp
#   endif
#   if defined uvic_cfc12
            sbc(i,j,icfc12flx) = sbc(i,j,icfc12flx) + cfc12ocn*tmp
#   endif
          endif
        enddo
      enddo
#   if defined uvic_carbon
      call setbcx (sbc(1,1,idicflx), imt, jmt)
#    if defined uvic_carbon_14
      call setbcx (sbc(1,1,ic14flx), imt, jmt)
#    endif
#   endif
#   if defined uvic_alk
      call setbcx (sbc(1,1,ialkflx), imt, jmt)
#   endif
#   if defined uvic_o2
      call setbcx (sbc(1,1,io2flx), imt, jmt)
#   endif
#   if defined uvic_npzd
      call setbcx (sbc(1,1,ipo4flx), imt, jmt)
#    if defined uvic_npzd_vflux
      call setbcx (sbc(1,1,iphytflx), imt, jmt)
      call setbcx (sbc(1,1,izoopflx), imt, jmt)
      call setbcx (sbc(1,1,idetrflx), imt, jmt)
#    endif
#    if defined uvic_nitrogen
      call setbcx (sbc(1,1,ino3flx), imt, jmt)
#     if defined uvic_npzd_vflux
      call setbcx (sbc(1,1,idiazflx), imt, jmt)
#     endif
#    endif
#   endif
#   if defined uvic_cfc11
      call setbcx (sbc(1,1,icfc11flx), imt, jmt)
#   endif
#   if defined uvic_cfc12
      call setbcx (sbc(1,1,icfc12flx), imt, jmt)
#   endif
#  endif
# endif
#endif

      return
      end
