      subroutine globalsum (ind)

!=======================================================================
!     calculate the global sums of heat and fresh water

!     input:
!       ind = index for calculation (1=initial, 2=current, 3=final)

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_global_sums
# if !defined uvic_mom

      implicit none
# endif

# include "param.h"
# if defined uvic_embm
#  include "atm.h"
#  include "cembm.h"
# endif
# if defined uvic_ice && defined uvic_embm
#  include "ice.h"
# endif
# include "coord.h"
# include "grdvar.h"
# include "levind.h"
# include "csbc.h"
# include "iounit.h"
# if defined uvic_mom
#  include "mw.h"
# endif
# if defined uvic_mtlm
#  include "mtlm.h"
#  include "switch.h"
# endif
# if defined uvic_npzd
#  include "npzd.h"
# endif
# include "tmngr.h"

      character(120) :: name, title, new_file_name, file_stamp
      save name
      data name /' '/

      integer i, iou, id_time, ind, j, jrow, k, l, n, nc, ntrec
      integer nyear, nmonth, nday, nhour, nmin, nsec, it(10)

      logical defined

      real area, cstdyt, dtoih1, dtoih2
      real dtah1, dtah2, dtaf1, dtaf2, dtac1, dtac2
      real dtsh1, dtsh2, dtsf1, dtsf2, dtsc1, dtsc2
      real dtih1, dtih2, dtif1, dtif2, dtic1, dtic2
      real dtlh1, dtlh2, dtlf1, dtlf2, dtlc1, dtlc2
      real dtoh1, dtof1, dtoh2, dtof2, dtoc1, dtoc2
      real dth1,  dth2,  dtf1,  dtf2,  dtc1,  dtc2
      real tah, taf, tac, tsh, tsf, tsc, tih, tif, tic
      real tlh, tlf, tlc, toh, tof, toc, th, tf, tc
      real RFTIME, time, tmp, c1e20

      common /global_r/ tah(3), taf(3), tac(3), tsh(3), tsf(3), tsc(3)
      common /global_r/ tih(3), tif(3), tic(3), tlh(3), tlf(3), tlc(3)
      common /global_r/ toh(3), tof(3), toc(3), dtoih1

!-----------------------------------------------------------------------
!     sum heat, water, carbon in atmosphere, snow, ice, land, ocean
!-----------------------------------------------------------------------

      c1e20 = 1.e20
      if (ind .eq. 1) dtoih1 = 0.0

!     zero accumulators

      do n=min(3,ind*2-1),3
        tah(n) = 0.
        taf(n) = 0.
        tac(n) = 0.
        tsh(n) = 0.
        tsf(n) = 0.
        tsc(n) = 0.
        tih(n) = 0.
        tif(n) = 0.
        tic(n) = 0.
        tlh(n) = 0.
        tlf(n) = 0.
        tlc(n) = 0.
        toh(n) = 0.
        tof(n) = 0.
        toc(n) = 0.
      enddo

# if !defined uvic_embm
      dtoih = 0.0
# endif
# if defined uvic_mtlm
#  if defined uvic_mtlm_segday
      RFTIME = DAY_TRIF*segtim/DAY_YEAR
      atlnd = segtim*SEC_DAY
#  else
      RFTIME = DAY_TRIF/DAY_YEAR
#  endif
# endif

!     calculate totals

      L = 0
      do jrow=2,jmtm1
# if defined uvic_mom
        if (wide_open_mw) then
          j = jrow
        else
          j = jmw
          call getrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
        endif
# endif
        cstdyt = cst(jrow)*dyt(jrow)
        do i=2,imtm1
          area = cstdyt*dxt(i)
          if (kmt(i,jrow) .eq. 0) L = L + 1

!         atmosphere
# if defined uvic_embm
          tah(3) = tah(3) + at(i,jrow,2,isat)*area
          taf(3) = taf(3) + at(i,jrow,2,ishum)*area
#  if defined uvic_carbon_co2_2d
          tac(3) = tac(3) + at(i,jrow,2,ico2)*area
#  else
          tac(3) = tac(3) + co2ccn*area
#  endif
# endif

!         snow and ice
# if defined uvic_ice && defined uvic_embm
#  if defined uvic_ice_cpts
          do nc=1,ncat
            tsf(3) = tsf(3) + hseff(i,jrow,2,nc)*area
            tif(3) = tif(3) + heff(i,jrow,2,nc)*area
          enddo
          do layer=1,ntilay
            tih(3) = tih(3) + E(i,jrow,2,layer)*area
          enddo
#  else
          tsf(3) = tsf(3) + hsno(i,jrow,2)*area
          tif(3) = tif(3) + hice(i,jrow,2)*area
#  endif
# endif
# if defined uvic_mtlm
          if (kmt(i,jrow) .eq. 0)
     &      tsf(3) = tsf(3) + LYING_SNOW(L)*area*0.1/rhosno
# endif

!         land
          tlf(3) = tlf(3) + soilm(i,jrow,2)*area*1.0e-3
# if defined uvic_mtlm
          if (kmt(i,jrow) .eq. 0) then
            tlh(3) = tlh(3) + (TS1(L)*HCAP_SOIL*ROOTDEP
!           add in most recent fluxes not felt by atmosphere yet
     &             + (sbc(i,jrow,ievap)*vlocn + sbc(i,jrow,isens)
     &             + sbc(i,jrow,ilwr))*atlnd*1.e-3)*area*1.e-4
            tlf(3) = tlf(3) + (M(L) + MNEG(L)
!           add back in most recent fluxes not felt by atmosphere yet
     &             + sbc(i,jrow,ievap)*atlnd*10.)*area*1.e-4
            tlc(3) = tlc(3) + (CV(L) + CS(L)
#  if defined uvic_carbon
!           add in most recent fluxes not felt by atmosphere yet
     &             + (sbc(i,jrow,isr) - sbc(i,jrow,inpp))*atlnd
!           add in driving fluxes not felt by triffid yet
     &             - RESP_S_DR(L)*RFTIME
#  endif
     &               )*area*1.e-4
            do N=1,NPFT
              tlc(3) = tlc(3) + NPP_DR(L,N)*FRAC(L,N)*RFTIME*area*1.0e-4
            enddo
          endif
# endif

!         ocean
# if defined uvic_mom
          do k=1,km
            vol = dzt(k)*area
            toh(3) = toh(3) + t(i,k,j,itemp,taup1)*vol
            tof(3) = tof(3) + t(i,k,j,isalt,taup1)*vol
#  if defined uvic_carbon
            toc(3) = toc(3) + t(i,k,j,idic,taup1)*vol
#   if defined uvic_npzd
            toc(3) = toc(3) + t(i,k,j,iphyt,taup1)*vol*redctn
            toc(3) = toc(3) + t(i,k,j,izoop,taup1)*vol*redctn
            toc(3) = toc(3) + t(i,k,j,idetr,taup1)*vol*redctn
#    if defined uvic_nitrogen
            toc(3) = toc(3) + t(i,k,j,idiaz,taup1)*vol*redctn
#    endif
#   endif
#  endif
          enddo
# endif
        enddo
      enddo

!     convert units to Joules and kilograms

!     atmosphere
# if defined uvic_embm
      tah(3) = (taf(3)*rhoatm*shq*vlocn+tah(3)*cpatm*rhoatm*sht)*1.0e-7
      taf(3) = taf(3)*rhoatm*shq*1.0e-3
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      tac(3) = tac(3)*4.138e-7*rhoatm*shc*1.0e-3
# endif

!     snow and ice
# if defined uvic_ice && defined uvic_embm
      tsh(3) = -tsf(3)*rhosno*flice*1.0e-7
      tsf(3) = tsf(3)*rhosno*1.0e-3
#  if defined uvic_ice_cpts
      tih(3) = tih(3)*1.0e-7   ! cpts ice energy of melting is negative
#  else
      tih(3) = -tif(3)*rhoice*flice*1.0e-7
#  endif
      tif(3) = tif(3)*rhoice*1.0e-3
# endif

!     land
      tlh(3) = tlh(3)
      tlf(3) = tlf(3)
      tlc(3) = tlc(3)

!     ocean
      toh(3) = toh(3)/0.2389
# if defined uvic_embm
      tof(3) = -tof(3)/(socn*1.e3)
# else
      tof(3) = -tof(3)/34.9
# endif
!     12e-9 kg/umol
      toc(3) = toc(3)*12.e-9

!     calculate differences from start of run

      dtah1 = tah(3) - tah(1)
      dtaf1 = taf(3) - taf(1)
      dtac1 = tac(3) - tac(1)
      dtsh1 = tsh(3) - tsh(1)
      dtsf1 = tsf(3) - tsf(1)
      dtsc1 = tsc(3) - tsc(1)
      dtih1 = tih(3) - tih(1)
      dtif1 = tif(3) - tif(1)
      dtic1 = tic(3) - tic(1)
      dtlh1 = tlh(3) - tlh(1)
      dtlf1 = tlf(3) - tlf(1)
      dtlc1 = tlc(3) - tlc(1)
      dtoh1 = toh(3) - toh(1)
      dtof1 = tof(3) - tof(1)
      dtoc1 = toc(3) - toc(1)
      dtoih1 = dtoih1 + dtoih

!     calculate differences from last time step

      dtah2 = tah(3) - tah(2)
      dtaf2 = taf(3) - taf(2)
      dtac2 = tac(3) - tac(2)
      dtsh2 = tsh(3) - tsh(2)
      dtsf2 = tsf(3) - tsf(2)
      dtsc2 = tsc(3) - tsc(2)
      dtih2 = tih(3) - tih(2)
      dtif2 = tif(3) - tif(2)
      dtic2 = tic(3) - tic(2)
      dtlh2 = tlh(3) - tlh(2)
      dtlf2 = tlf(3) - tlf(2)
      dtlc2 = tlc(3) - tlc(2)
      dtoh2 = toh(3) - toh(2)
      dtof2 = tof(3) - tof(2)
      dtoc2 = toc(3) - toc(2)
      dtoih2 = dtoih

      dtoih = 0.0
      do n=ind,2
        tah(n) = tah(3)
        taf(n) = taf(3)
        tac(n) = tac(3)
        tsh(n) = tsh(3)
        tsf(n) = tsf(3)
        tsc(n) = tsc(3)
        tih(n) = tih(3)
        tif(n) = tif(3)
        tic(n) = tic(3)
        tlh(n) = tlh(3)
        tlf(n) = tlf(3)
        tlc(n) = tlc(3)
        toh(n) = toh(3)
        tof(n) = tof(3)
        toc(n) = toc(3)
      enddo

!     write differences or totals

      if (ind .eq. 2) then

        dth1 = dtah1 + dtsh1 + dtih1 + dtlh1 + dtoh1 + dtoih1
        dtf1 = dtaf1 + dtsf1 + dtif1 + dtlf1 + dtof1
        dtc1 = dtac1 + dtsc1 + dtic1 + dtlc1 + dtoc1
        dth2 = dtah2 + dtsh2 + dtih2 + dtlh2 + dtoh2 + dtoih2
        dtf2 = dtaf2 + dtsf2 + dtif2 + dtlf2 + dtof2
        dtc2 = dtac2 + dtsc2 + dtic2 + dtlc2 + dtoc2

        if (iotsi .ne. stdout) then
          time = relyr + year0
          defined = .true.
          ntrec = 1
          if (name .eq. ' ') then
            defined = .false.
            ntrec = 0
            name = file_stamp ('ts_intgrls',stamp,'.nc')
            name = new_file_name (name)
          endif
          call opennext (name, time, ntrec, iou)
          if (.not. defined .or. ntrec .eq. 0) then
            call redef (iou)
            call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
            call putatttext (iou, 'global', 'experiment_name', expnam)
            call putatttext (iou, 'global', 'run_stamp', runstamp)
            call putatttext (iou, 'global', 'time_unit', timunit)
            call defdim ('time', iou, 0, id_time)
            it(:) = id_time
            call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &,       'time since initial condition', 'time', trim(timunit))
            call defvar ('year', iou, 1, it, c0, c0, ' '
     &,       'F', 'year', ' ',' ')
            call defvar ('month', iou, 1, it, c0, c0, ' '
     &,       'F', 'month', ' ',' ')
            call defvar ('day', iou, 1, it, c0, c0, ' '
     &,       'F', 'day', ' ',' ')
            call defvar ('hour', iou, 1, it, c0, c0, ' '
     &,       'F', 'hour', ' ',' ')
            call defvar ('minute', iou, 1, it, c0, c0, ' '
     &,       'F', 'minute', ' ',' ')
            call defvar ('second', iou, 1, it, c0, c0, ' '
     &,       'F', 'second', ' ',' ')
            call defvar ('dtheat', iou, 1, it, -c1e20
     &,       c1e20, ' ', 'F', 'heat conservation error', ' ', 'J')
            call defvar ('dtwater', iou, 1, it, -c1e20
     &,       c1e20, ' ', 'F', 'water conservation error', ' ', 'kg')
            call defvar ('dtcarbon', iou, 1, it, -c1e20
     &,       c1e20, ' ', 'F', 'carbon conservation error', ' ', 'kg')
            call enddef (iou)
            if (ntrec .eq. 0 ) ntrec = 1
          endif

          call putvars ('time', iou, ntrec, time, c1, c0)
          call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
          tmp = nyear
          call putvars ('year', iou, ntrec, tmp, c1, c0)
          tmp = nmonth
          call putvars ('month', iou, ntrec, tmp, c1, c0)
          tmp = nday
          call putvars ('day', iou, ntrec, tmp, c1, c0)
          tmp = nhour
          call putvars ('hour', iou, ntrec, tmp, c1, c0)
          tmp = nmin
          call putvars ('minute', iou, ntrec, tmp, c1, c0)
          tmp = nsec
          call putvars ('second', iou, ntrec, tmp, c1, c0)
          call putvars ('dtheat', iou, ntrec, dth1, c1, c0)
          call putvars ('dtwater', iou, ntrec, dtf1, c1, c0)
          call putvars ('dtcarbon', iou, ntrec, dtc1, c1, c0)
        endif
        call closefile (iou)

        if (iotsi .eq. stdout .or. iotsi .lt. 0) then
          write (*,'(a,a)') 'Changes in heat, fresh water and carbon '
     &,     'from start of run and last calculation'
          write (*,'(a,e22.14,a,e22.14,a,e22.14,a
     &,                e22.14,a,e22.14,a,e22.14,a)')
     &      ' d atm    ',dtah1, ' J ',dtaf1,' kg ',dtac1,' kg'
     &,                  dtah2, ' J ',dtaf2,' kg ',dtac2,' kg'
     &,     ' d snow   ',dtsh1, ' J ',dtsf1,' kg ',dtsc1,' kg'
     &,                  dtsh2, ' J ',dtsf2,' kg ',dtsc2,' kg'
     &,     ' d ice    ',dtih1, ' J ',dtif1,' kg ',dtic1,' kg'
     &,                  dtih2, ' J ',dtif2,' kg ',dtic2,' kg'
     &,     ' d lnd    ',dtlh1, ' J ',dtlf1,' kg ',dtlc1,' kg'
     &,                  dtlh2, ' J ',dtlf2,' kg ',dtlc2,' kg'
     &,     ' d ocn    ',dtoh1, ' J ',dtof1,' kg ',dtoc1,' kg'
     &,                  dtoh2, ' J ',dtof2,' kg ',dtoc2,' kg'
     &,     ' d out-in ',dtoih1,' J ',0.0,  ' kg ',0.0,  ' kg'
     &,                  dtoih2,' J ',0.0,  ' kg ',0.0,  ' kg'
     &,     ' d total  ',dth1,  ' J ',dtf1, ' kg ',dtc1, ' kg'
     &,                  dth2,  ' J ',dtf2, ' kg ',dtc2, ' kg'
        endif

      else

        th = tah(ind) + tsh(ind) + tih(ind) + tlh(ind) + toh(ind)
        tf = taf(ind) + tsf(ind) + tif(ind) + tlf(ind) + tof(ind)
        tc = tac(ind) + tsc(ind) + tic(ind) + tlc(ind) + toc(ind)
        write (*,'(/,a,a)') 'Total heat (in Joules referenced to 0 C'
     &,   ' and no ice or snow) and fresh water (in kg)'
# if defined uvic_embm
        write (*,'(a,a)') 'Total ocean fresh water is the equivalent'
     &,   ' difference from the ocean volume referenced to socn'
# else
        write (*,'(a,a)') 'Total ocean fresh water is the equivalent'
     &,   ' difference from the ocean volume referenced to 34.9 ppt'
# endif
        write (*,'(a,e22.14,a,e22.14,a,e22.14,a)')
     &    ' t atm    ', tah(ind),' J ',taf(ind),' kg ',tac(ind),' kg'
     &,   ' t snow   ', tsh(ind),' J ',tsf(ind),' kg ',tsc(ind),' kg'
     &,   ' t ice    ', tih(ind),' J ',tif(ind),' kg ',tic(ind),' kg'
     &,   ' t lnd    ', tlh(ind),' J ',tlf(ind),' kg ',tlc(ind),' kg'
     &,   ' t ocn    ', toh(ind),' J ',tof(ind),' kg ',toc(ind),' kg'
     &,   ' t total  ', th,      ' J ',tf,      ' kg ',tc,      ' kg'

      endif
#endif

      return
      end
