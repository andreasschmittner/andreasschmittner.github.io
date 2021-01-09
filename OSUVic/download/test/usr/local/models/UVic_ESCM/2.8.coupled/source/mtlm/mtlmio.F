      subroutine mtlmout (is, ie, js, je)
!-----------------------------------------------------------------------
!     Output routine for the mtlm

!     based on code by: K. Meissner and M. Eby
!-----------------------------------------------------------------------
#if defined uvic_mtlm

      implicit none

# include "param.h"
# include "calendar.h"
# include "coord.h"
# include "grdvar.h"
# include "mtlm.h"
# include "csbc.h"
# include "iounit.h"
# include "switch.h"
# include "tmngr.h"

      character(120) :: fname

      integer is, ie, js, je, ntrec

      real avgper, timeout

# if defined time_step_monitor
      if (tsits .and. ntatil .ne. 0) then

        call ta_mtlm_tsi (2)

        if (iotsi .ne. stdout) then
#  if defined uvic_mtlm_segday
          avgper = segtim*dtlnd*ntatil/daylen
#  else
          avgper = dtlnd*ntatil/daylen
#  endif
          timeout = relyr + year0
          call def_tsi
          call def_tsi_mtlm (fname)
          call mtlm_tsi_out (fname, avgper, timeout, stamp, tai_CS
     &,     tai_RESP_S, tai_CV, tai_NPP, tai_GPP, tai_HT, tai_LAI
     &,     tai_LYING_SNOW, tai_TSOIL, tai_TSTAR, tai_M, tai_ET, ntrec)
        endif

        call ta_mtlm_tsi (0)

      endif
# endif
# if defined time_averages
      if (timavgts .and. ntatsl .ne. 0) then

!-----------------------------------------------------------------------
!       write atmospheric time averaged data
!-----------------------------------------------------------------------

!       calculate average values

        call ta_mtlm_snap (2)

!       write time averaged data

#  if defined uvic_mtlm_segday
        avgper = segtim*dtlnd*ntatsl/daylen
#  else
        avgper = dtlnd*ntatsl/daylen
#  endif
        timeout = relyr + year0
        call def_tavg
        call def_tavg_mtlm (fname)
        call mtlm_snap_out (fname, is, ie, js, je, imt, jmt
     &,   POINTS, NPFT, NTYPE, xt, yt, xu, yu, dxt, dyt, dxu, dyu
     &,   avgper, timeout, stamp, land_map, ta_TS1, ta_CS, ta_RESP_S
     &,   ta_FRAC, ta_GPP, ta_NPP, ta_HT, ta_LAI, ta_C_VEG
# if !defined uvic_embm
     &,   ta_TSTAR_GB, ta_ALBLAND, ta_ET, ta_M
# endif
     &,   tlat, tlon, ulat, ulon, ntrec)

        write (*,'(a,i5,a,a,a,a)') '=> Lnd time means #'
     &,   ntrec, ' written to ',trim(fname),' on ', stamp

!       zero time average accumulators

        call ta_mtlm_snap (0)

      endif

# endif

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_mtlm (0, fname)
          call mtlm_rest_out (fname, is, ie, js, je)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_mtlm (1, fname)
          call mtlm_rest_out (fname, is, ie, js, je)
        endif
      endif

      return
      end

      subroutine ta_mtlm_snap (i)

!=======================================================================
!     land data time averaging

!     input:
!       i = flag (0 = zero, 1 = accumulate, 2 = write)

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "size.h"
# include "mtlm.h"
# include "csbc.h"

      integer i, n

      real rntatsl

!-----------------------------------------------------------------------
!     time averaged data
!-----------------------------------------------------------------------

      if (i .eq. 0.) then
!       zero
        ntatsl = 0
        ta_TS1(:) = 0.
        ta_TSTAR_GB(:) = 0.
# if !defined uvic_embm
        ta_ALBLAND(:) = 0.
        ta_ET(:) = 0.
# endif
        ta_M(:) = 0.
        ta_CS(:) = 0.
        ta_RESP_S(:) = 0.
        ta_GPP(:,:) = 0.
        ta_NPP(:,:) = 0.
        ta_HT(:,:) = 0.
        ta_LAI(:,:) = 0.
        ta_C_VEG(:,:) = 0.
        ta_LYING_SNOW(:) = 0.
        ta_SURF_ROFF(:) = 0.
        ta_FRAC(:,:) = 0.
      elseif (i .eq. 1) then
!       accumulate
        ntatsl = ntatsl + 1
        ta_TS1(:) = ta_TS1(:) + TS1(:)
        ta_TSTAR_GB(:) = ta_TSTAR_GB(:) + TSTAR_GB(:)
# if !defined uvic_embm
        ta_ALBLAND(:) = ta_ALBLAND(:) + ALBLAND(:)
        ta_ET(:) = ta_ET(:) + ET(:)
# endif
        ta_M(:) = ta_M(:) + M(:)
        ta_CS(:) = ta_CS(:) + CS(:)
        ta_RESP_S(:) = ta_RESP_S(:) + RESP_S(:)
        do n=1,NPFT
          ta_GPP(:,n) = ta_GPP(:,n) + GPP(:,n)*FRAC(:,n)
          ta_NPP(:,n) = ta_NPP(:,n) + NPP(:,n)*FRAC(:,n)
          ta_HT(:,n) = ta_HT(:,n) + HT(:,n)*FRAC(:,n)
          ta_LAI(:,n) = ta_LAI(:,n) + LAI(:,n)*FRAC(:,n)
          ta_C_VEG(:,n) = ta_C_VEG(:,n) + C_VEG(:,n)*FRAC(:,n)
        enddo
        ta_LYING_SNOW(:) = ta_LYING_SNOW(:) + LYING_SNOW(:)
        ta_SURF_ROFF(:) = ta_SURF_ROFF(:) + SURF_ROFF(:)
        ta_FRAC(:,:) = ta_FRAC(:,:) + FRAC(:,:)

      elseif (i .eq. 2 .and. ntatsl .ne. 0) then
!       average
        rntatsl = 1./float(ntatsl)
        ta_TS1(:) = ta_TS1(:)*rntatsl
        ta_TSTAR_GB(:) = ta_TSTAR_GB(:) *rntatsl
# if !defined uvic_embm
        ta_ALBLAND(:) = ta_ALBLAND(:)*rntatsl
        ta_ET(:) = ta_ET(:)*rntatsl
# endif
        ta_M(:) = ta_M(:)*rntatsl
        ta_CS(:) = ta_CS(:)*rntatsl
        ta_RESP_S(:) = ta_RESP_S(:)*rntatsl
        ta_GPP(:,:) = ta_GPP(:,:)*rntatsl
        ta_NPP(:,:) = ta_NPP(:,:)*rntatsl
        ta_HT(:,:) = ta_HT(:,:)*rntatsl
        ta_LAI(:,:) = ta_LAI(:,:)*rntatsl
        ta_C_VEG(:,:) = ta_C_VEG(:,:)*rntatsl
        ta_LYING_SNOW(:) = ta_LYING_SNOW(:)*rntatsl
        ta_SURF_ROFF(:) = ta_SURF_ROFF(:)*rntatsl
        ta_FRAC(:,:) = ta_FRAC(:,:)*rntatsl
      endif

      return
      end

      subroutine ta_mtlm_tsi (i)

!=======================================================================
!     atmospheric data time integral averaging

!     input:
!       i = flag (0 = zero, 1 = accumulate, 2 = write)

!     based on code by: M. Eby
!=======================================================================

      implicit none

# include "size.h"
# include "csbc.h"
# include "mtlm.h"

      integer i, n

      real rntatil, data(imt,jmt), dmsk(imt,jmt), wt(imt, jmt), tmp

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

       if (i .eq. 0.) then
!       zero
        ntatil = 0
        tai_CS = 0
        tai_RESP_S = 0
        tai_CV = 0
        tai_NPP = 0
        tai_GPP = 0
        tai_HT = 0
        tai_LAI = 0
        tai_LYING_SNOW = 0
        tai_TSOIL = 0
        tai_TSTAR = 0
        tai_M = 0
        tai_ET = 0
      elseif (i .eq. 1) then
!       set data mask
        dmsk(:,:) = 1.
        where (land_map(:,:) .eq. 0) dmsk(:,:) = 0.
!       accumulate
        ntatil = ntatil + 1
        call unloadland (POINTS, CS, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_CS = tai_CS + tmp*1.e-4
        call unloadland (POINTS, RESP_S, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_RESP_S = tai_RESP_S + tmp*1.e-4
        call unloadland (POINTS, CV, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_CV = tai_CV + tmp*1.e-4
        do n=1,NPFT
          call unloadland (POINTS, FRAC(1,n), imt, jmt, land_map, wt)
          call unloadland (POINTS, NPP(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areatot (data, dmsk, tmp)
!         convert area to m2
          tai_NPP = tai_NPP + tmp*1.e-4
          call unloadland (POINTS, GPP(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areatot (data, dmsk, tmp)
!         convert area to m2
          tai_GPP = tai_GPP + tmp*1.e-4
          call unloadland (POINTS, HT(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areaavg (data, dmsk, tmp)
          tai_HT = tai_HT + tmp
          call unloadland (POINTS, LAI(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areaavg (data, dmsk, tmp)
          tai_LAI = tai_LAI + tmp
        enddo
        call unloadland (POINTS, LYING_SNOW, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_LYING_SNOW = tai_LYING_SNOW + tmp*1.e-4
        call unloadland (POINTS, TSOIL, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_TSOIL = tai_TSOIL + tmp
        call unloadland (POINTS, TSTAR_GB, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_TSTAR = tai_TSTAR + tmp
        call unloadland (POINTS, M, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_M = tai_M + tmp
        call unloadland (POINTS, ET, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_ET = tai_ET + tmp
      elseif (i .eq. 2 .and. ntatil .ne. 0) then
!       average
        rntatil = 1./float(ntatil)
        tai_CS = tai_CS*rntatil
        tai_RESP_S = tai_RESP_S*rntatil
        tai_CV = tai_CV*rntatil
        tai_NPP = tai_NPP*rntatil
        tai_GPP = tai_GPP*rntatil
        tai_HT = tai_HT*rntatil
        tai_LAI = tai_LAI*rntatil
        tai_LYING_SNOW = tai_LYING_SNOW*rntatil
        tai_TSOIL = tai_TSOIL*rntatil
        tai_TSTAR = tai_TSTAR*rntatil
        tai_M = tai_M*rntatil
        tai_ET = tai_ET*rntatil
      endif
#endif

      return
      end

      subroutine unloadland (ld, dl, id, jd, map, dij)

      integer i, id, j, jd, l, ld, map(id,jd)
      real dl(ld), dij(id,jd)

      dij(:,:) = 0.
      do j=1,jd
        do i=1,id
          l = map(i,j)
          if (l .ge. 1 .and. l .le. ld) dij(i,j) = dl(l)
        enddo
      enddo

      return
      end

      subroutine loadland (ld, dl, id, jd, map, dij)

      integer i, id, j, jd, l, ld, map(id,jd)
      real dl(ld), dij(id,jd)

      dl(:) = 0.
      do j=1,jd
        do i=1,id
          l = map(i,j)
          if (l .ge. 1 .and. l .le. ld) dl(l) = dij(i,j)
        enddo
      enddo

      return
      end
