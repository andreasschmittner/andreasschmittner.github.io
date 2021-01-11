      program UVic_ESCM

!=======================================================================

!     UNIVERSITY OF VICTORIA EARTH SYSTEM CLIMATE MODEL

!     A climate model developed by researchers in the Climate Research
!     Group, in the School of Earth and Ocean Sciences, located at the
!     University of Victoria, Victoria, B.C., Canada.

!     Many people have contributed to the development of this code
!     and attempts are made to indicate originators of code where
!     possible or appropriate. Please direct problems or questions
!     to the code contact person at:
!     http://climate.uvic.ca/climate-lab/model

!     Requirements:

!     Standard fortran 90 is used

!     Disclaimer:

!     The UVic Earth System Climate Model (UVic_ESCM) is a climate
!     modeling research tool developed at the University of Victoria.
!     Others may use it freely but we assume no responsibility for
!     problems or incorrect use. It is left to the user to ensure that
!     a particular configuration is working correctly.
!=======================================================================

!     This is the main driver. Integration time is divided into a number
!     of equal time segments and SBC are held fixed for each time
!     segment. When coupling, SBC are supplied each time segment (the
!     coupling period) and held fixed for that period.

!     based on code by: R. C. Pacanowski, A. Rosati and M. Eby
!=======================================================================

#include "param.h"
#if defined uvic_embm
# include "atm.h"
# include "cembm.h"
#endif
#include "coord.h"
#include "csbc.h"
#include "iounit.h"
#include "levind.h"
#if defined uvic_mom
# include "mw.h"
# include "task_on.h"
#endif
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"

      logical   mk_out
      namelist /uvic/ mk_out

      print*, '== UNIVERSITY OF VICTORIA EARTH SYSTEM CLIMATE MODEL =='
      print*, '                                            '

#if defined uvic_mk_out
!-----------------------------------------------------------------------
!     mk options
!-----------------------------------------------------------------------

      mk_out = .false.
      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, uvic, end=100)
100   continue
      call relunit (ioun)
      if (mk_out) then
!       print mk options
        print*, ' '
        print*, 'Mk Options are:'
        include "uvic_mk.h"
        print*, ' '
      endif

#endif
!-----------------------------------------------------------------------
!     initialize i/o units
!-----------------------------------------------------------------------

      call ioinit

!-----------------------------------------------------------------------
!     setup file renaming
!-----------------------------------------------------------------------

      call file_names

!-----------------------------------------------------------------------
!     Initialize S.B.C. indices
!-----------------------------------------------------------------------

      sbc(:,:,:) = 0.0
      mapsbc(:) = " "
      itaux = 0
      itauy = 0
      iws = 0
      iaca = 0
      isca = 0
      ihflx = 0
      isflx = 0
      isst = 0
      isss = 0
      iro = 0
      iwa = 0
      iwxq = 0
      iwyq = 0
      iwxt = 0
      iwyt = 0
      ipsw = 0
      isu = 0
      isv = 0
      igu = 0
      igv = 0
      issdic = 0
      idicflx = 0
      issalk = 0
      ialkflx = 0
      isso2 = 0
      io2flx = 0
      isspo4 = 0
      ipo4flx = 0
      issphyt = 0
      iphytflx = 0
      isszoop = 0
      izoopflx = 0
      issdetr = 0
      idetrflx = 0
      issno3 = 0
      ino3flx = 0
      issdiaz = 0
      idiazflx = 0
      issc14 = 0
      ic14flx = 0
      isscfc11 = 0
      icfc11flx = 0
      isscfc12 = 0
      icfc12flx = 0
      iat = 0
      irh = 0
      ipr = 0
      ips = 0
      iaws = 0
      iswr = 0
      ilwr = 0
      isens = 0
      ievap = 0
      idtr = 0
      isr = 0
      inpp = 0

      m = 1
      call set (itaux, m, mapsbc(m), 'taux', m)
      call set (itauy, m, mapsbc(m), 'tauy', m)
      call set (iws, m, mapsbc(m), 'ws', m)
      call set (iaca, m, mapsbc(m), 'a_calb', m)
      call set (isca, m, mapsbc(m), 's_calb', m)
      call set (ihflx, m, mapsbc(m), 'hflx', m)
      call set (isflx, m, mapsbc(m), 'sflx', m)
      call set (isst, m, mapsbc(m), 'sst', m)
      call set (isss, m, mapsbc(m), 'sss', m)
      call set (iro, m, mapsbc(m), 'ro', m)
#if defined uvic_embm_awind && defined uvic_embm
      call set (iwa, m, mapsbc(m), 'wa', m)
#endif
#if defined uvic_embm_adv_q && defined uvic_embm
      call set (iwxq, m, mapsbc(m), 'wx_q', m)
      call set (iwyq, m, mapsbc(m), 'wy_q', m)
#endif
#if defined uvic_embm_adv_t && defined uvic_embm
      call set (iwxt, m, mapsbc(m), 'wx_t', m)
      call set (iwyt, m, mapsbc(m), 'wy_t', m)
#endif
#if defined shortwave
      call set (ipsw, m, mapsbc(m), 'psw', m)
#endif
#if defined uvic_ice_evp
      call set (isu, m, mapsbc(m), 'su', m)
      call set (isv, m, mapsbc(m), 'sv', m)
      call set (igu, m, mapsbc(m), 'gu', m)
      call set (igv, m, mapsbc(m), 'gv', m)
#endif
#if defined uvic_carbon
      call set (issdic, m, mapsbc(m), 'ssdic', m)
      call set (idicflx, m, mapsbc(m), 'dicflx', m)
# if defined uvic_carbon_14
      call set (issc14, m, mapsbc(m), 'ssc14', m)
      call set (ic14flx, m, mapsbc(m), 'c14flx', m)
# endif
#endif
#if defined uvic_alk
      call set (issalk, m, mapsbc(m), 'ssalk', m)
      call set (ialkflx, m, mapsbc(m), 'alkflx', m)
#endif
#if defined uvic_o2
      call set (isso2, m, mapsbc(m), 'sso2', m)
      call set (io2flx, m, mapsbc(m), 'o2flx', m)
#endif
#if defined uvic_npzd
      call set (isspo4, m, mapsbc(m), 'sspo4', m)
      call set (ipo4flx, m, mapsbc(m), 'po4flx', m)
# if defined uvic_npzd_vflux
      call set (issphyt, m, mapsbc(m), 'ssphyt', m)
      call set (iphytflx, m, mapsbc(m), 'phytflx', m)
      call set (isszoop, m, mapsbc(m), 'sszoop', m)
      call set (izoopflx, m, mapsbc(m), 'zoopflx', m)
      call set (issdetr, m, mapsbc(m), 'ssdetr', m)
      call set (idetrflx, m, mapsbc(m), 'detrflx', m)
# endif
# if defined uvic_nitrogen
      call set (issno3, m, mapsbc(m), 'ssno3', m)
      call set (ino3flx, m, mapsbc(m), 'no3flx', m)
#  if defined uvic_npzd_vflux
      call set (issdiaz, m, mapsbc(m), 'ssdiaz', m)
      call set (idiazflx, m, mapsbc(m), 'diazflx', m)
#  endif
# endif
#endif
#if defined uvic_cfc11
      call set (isscfc11, m, mapsbc(m), 'sscfc11', m)
      call set (icfc11flx, m, mapsbc(m), 'cfc11flx', m)
#endif
#if defined uvic_cfc12
      call set (isscfc12, m, mapsbc(m), 'sscfc12', m)
      call set (icfc12flx, m, mapsbc(m), 'cfc12flx', m)
#endif
#if defined uvic_mtlm
      call set (iat, m, mapsbc(m), 'at', m)
      call set (irh, m, mapsbc(m), 'rh', m)
      call set (ipr, m, mapsbc(m), 'pr', m)
      call set (ips, m, mapsbc(m), 'ps', m)
      call set (iaws, m, mapsbc(m), 'aws', m)
      call set (iswr, m, mapsbc(m), 'swr', m)
      call set (ilwr, m, mapsbc(m), 'lwr', m)
      call set (isens, m, mapsbc(m), 'sens', m)
      call set (ievap, m, mapsbc(m), 'evap', m)
      call set (idtr, m, mapsbc(m), 'dtr', m)
#endif
#if defined uvic_mtlm && defined uvic_carbon
      call set (isr, m, mapsbc(m), 'sr', m)
      call set (inpp, m, mapsbc(m), 'npp', m)
#endif

      if ( m-1 .gt. numsbc) then
        print*, '==> Error: increase numsbc in csbc.h to ', m-1
        stop '=>UVic_ESCM'
      endif

#if defined uvic_mom
!-----------------------------------------------------------------------
!     Initialize ocean tracer names
!-----------------------------------------------------------------------

      do n=1,nt
        mapt(n) = " "
      enddo
      itemp = 0
      isalt = 0
      idic = 0
      ic14 = 0
      icfc11 = 0
      icfc12 = 0
      io2 = 0
      ialk = 0
      ipo4 = 0
      ino3 = 0
      iphyt = 0
      izoop = 0
      idetr = 0
      idiaz = 0

      m = 1
      call set (itemp, m, mapt(m), 'temp', m)
      call set (isalt, m, mapt(m), 'salt', m)
# if defined uvic_carbon
      call set (idic, m, mapt(m), 'dic', m)
#  if defined uvic_carbon_14
      call set (ic14, m, mapt(m), 'c14', m)
#  endif
# endif
# if defined uvic_cfc11
      call set (icfc11, m, mapt(m), 'cfc11', m)
# endif
# if defined uvic_cfc12
      call set (icfc12, m, mapt(m), 'cfc12', m)
# endif
# if defined uvic_o2
      call set (io2, m, mapt(m), 'o2', m)
# endif
# if defined uvic_alk
      call set (ialk, m, mapt(m), 'alk', m)
# endif
# if defined uvic_npzd
      call set (ipo4, m, mapt(m), 'po4', m)
      call set (iphyt, m, mapt(m), 'phyt', m)
      call set (izoop, m, mapt(m), 'zoop', m)
      call set (idetr, m, mapt(m), 'detr', m)
#  if defined uvic_nitrogen
      call set (ino3, m, mapt(m), 'no3', m)
      call set (idiaz, m, mapt(m), 'diaz', m)
#  endif
# endif
      if ( m-1 .gt. nt) then
        print*, '==> Error: increase nt for tracers in size.h'
        stop '=>UVic_ESCM'
      endif

!-----------------------------------------------------------------------
!     Initialize ocean source tracer names, must have equivalent tracer
!-----------------------------------------------------------------------

      mapst(:) = " "
      itrc(:) = 0

      m = 1
# if defined uvic_carbon && defined uvic_carbon_14
      call set (isc14, m, mapst(m), 'sc14', m)
      itrc(ic14) = m-1
# endif
# if defined uvic_npzd
#  if defined uvic_carbon
      call set (isdic, m, mapst(m), 'sdic', m)
      itrc(idic) = m-1
#  endif
#  if defined uvic_o2
      call set (iso2, m, mapst(m), 'so2', m)
      itrc(io2) = m-1
#  endif
#  if defined uvic_alk
      call set (isalk, m, mapst(m), 'salk', m)
      itrc(ialk) = m-1
#  endif
      call set (ispo4, m, mapst(m), 'spo4', m)
      itrc(ipo4) = m-1
      call set (isphyt, m, mapst(m), 'sphyt', m)
      itrc(iphyt) = m-1
      call set (iszoop, m, mapst(m), 'szoop', m)
      itrc(izoop) = m-1
      call set (isdetr, m, mapst(m), 'sdetr', m)
      itrc(idetr) = m-1
#  if defined uvic_nitrogen
      call set (isno3, m, mapst(m), 'sno3', m)
      itrc(ino3) = m-1
      call set (isdiaz, m, mapst(m), 'sdiaz', m)
      itrc(idiaz) = m-1
#  endif
# endif
      if ( m-1 .gt. nt) then
        print*, '==> Error: increase nsrc for tracer sources in size.h'
        stop '=>UVic_ESCM'
      endif

#endif
#if defined uvic_embm
!-----------------------------------------------------------------------
!     Initialize atmosphere tracer names
!-----------------------------------------------------------------------

      do n=1,nat
        mapat(n) = " "
      enddo
      isat = 0
      ishum = 0
      ico2 = 0

      m = 1
      call set (isat, m, mapat(m), 'sat', m)
      call set (ishum, m, mapat(m), 'shum', m)
# if defined uvic_carbon && defined uvic_carbon_co2_2d
      call set (ico2, m, mapat(m), 'co2', m)
# endif

      if ( m-1 .gt. nat) then
        print*, '==> Error: increase nat in size.h'
        stop '=>UVic_ESCM'
      endif
#endif

!-----------------------------------------------------------------------
!     do the introductory ocean setup once per run
!-----------------------------------------------------------------------

      call setocn

!-----------------------------------------------------------------------
!     do the introductory atmosphere setup once per run
!-----------------------------------------------------------------------

      write (stdout,'(/a36/)') ' ==> Note: the atmos setup follows:'

!     "setatm" must do the following:
!       1) set up the atmospheric S.B.C. grid definition
!       2) define the atmosphere land/sea mask
!       3) set the atmosphere time step "dtatm" {seconds}

#if defined uvic_embm
      call setembm
#else
      call setatm

!     when the MOM S.B.C. come from a dataset, force the segment time
!     and atmospheric time step to one MOM time step. This will force
!     the number of segments to one and the number of time steps per
!     segment to represent the length of the run in days.

      dtatm = dtocn
      segtim = dtocn*secday
#endif
#if defined uvic_mtlm

!-----------------------------------------------------------------------
!     do the introductory land setup once per run
!-----------------------------------------------------------------------
      call setmtlm
#endif

!-----------------------------------------------------------------------
!     compute the number of ocean time steps "numots" for this run and
!     the number of ocean time steps per ocean segment "ntspos".
!     compute the number of atmos time steps "numats" for this run and
!     the number of atmos time steps per atmos segment "ntspas".
!     divide the integration time "days" into "numseg" segments.
!     each will be length "segtim" days. Surface boundary conditions
!     are supplied every "segtim" days.
!-----------------------------------------------------------------------

      numots = nint(rundays/(dtocn*secday))
      ntspos = nint(segtim/(dtocn*secday))
      numats = nint(rundays/(dtatm*secday))
      ntspas = nint(segtim/(dtatm*secday))
      numseg = numots/ntspos
#if defined uvic_mtlm
# if defined uvic_mtlm_segday
      if (segtim .gt. 1.) then
        ntspls = nint(c1/(dtlnd*secday))
      else
        ntspls = nint(segtim/(dtlnd*secday))
      endif
# else
      ntspls = nint(segtim/(dtlnd*secday))
# endif
#endif

#if !defined uvic_embm && defined uvic_mom
!     load the tracers (SST & SSS) for each row "j". (zero on land)
!     load from the MW if fully opened otherwise load from disk

      if (wide_open_mw) then
        do j=1,jmt
          do i=1,imt
# if !defined uvic_replacst
            if (isst .ne. 0) sbc(i,j,isst) = t(i,1,j,itemp,taup1)
            if (isss .ne. 0) sbc(i,j,isss) = t(i,1,j,isalt,taup1)
# endif
# if defined uvic_carbon
            if (issdic .ne. 0) sbc(i,j,issdic) = t(i,1,j,idic,taup1)
#  if defined uvic_carbon_14
            if (issc14 .ne. 0) sbc(i,j,issc14) = t(i,1,j,ic14,taup1)
#  endif
# endif
# if defined uvic_alk
            if (issalk .ne. 0) sbc(i,j,issalk) = t(i,1,j,ialk,taup1)
# endif
# if defined uvic_o2
            if (isso2 .ne. 0) sbc(i,j,isso2) = t(i,1,j,io2,taup1)
# endif
# if defined uvic_npzd
            if (isspo4 .ne. 0) sbc(i,j,isspo4) = t(i,1,j,ipo4,taup1)
            if (issphyt .ne. 0) sbc(i,j,issphyt) = t(i,1,j,iphyt,taup1)
            if (isszoop .ne. 0) sbc(i,j,isszoop) = t(i,1,j,izoop,taup1)
            if (issdetr .ne. 0) sbc(i,j,issdetr) = t(i,1,j,idetr,taup1)
#  if defined uvic_nitrogen
            if (issno3 .ne. 0) sbc(i,j,issno3) = t(i,1,j,ino3,taup1)
            if (issdiaz .ne. 0) sbc(i,j,issdiaz) = t(i,1,j,idiaz,taup1)
#  endif
# endif
# if defined uvic_cfc11
            if (isscfc11 .ne. 0) sbc(i,j,isscfc11) =
     &        t(i,1,j,icfc11,taup1)
# endif
# if defined uvic_cfc12
            if (isscfc12 .ne. 0) sbc(i,j,isscfc12) =
     &        t(i,1,j,icfc12,taup1)
# endif
          enddo
        enddo
      else
        do j=1,jmt
# if !defined uvic_replacst
          if (isst .ne. 0) call getst (j, sbc(1,1,isst), itemp)
          if (isss .ne. 0) call getst (j, sbc(1,1,isss), isalt)
# endif
# if defined uvic_carbon
          if (issdic .ne. 0) call getst (j, sbc(1,1,issdic), idic)
#  if defined uvic_carbon_14
          if (issc14 .ne. 0) call getst (j, sbc(1,1,issc14), ic14)
#  endif
# endif
# if defined uvic_alk
          if (issalk .ne. 0) call getst (j, sbc(1,1,issalk), ialk)
# endif
# if defined uvic_o2
          if (isso2 .ne. 0) call getst (j, sbc(1,1,isso2), io2)
# endif
# if defined uvic_npzd
            if (isspo4 .ne. 0) call getst (j, sbc(1,1,isspo4), ipo4)
            if (issphyt .ne. 0) call getst (j, sbc(1,1,issphyt), iphyt)
            if (isszoop .ne. 0) call getst (j, sbc(1,1,isszoop), izoop)
            if (issdetr .ne. 0) call getst (j, sbc(1,1,issdetr), idetr)
#  if defined uvic_nitrogen
            if (issno3 .ne. 0) call getst (j, sbc(1,1,issno3), ino3)
            if (issdiaz .ne. 0) call getst (j, sbc(1,1,issdiaz), idiaz)
#  endif
# endif
# if defined uvic_cfc11
            if (isscfc11 .ne. 0) call getst (j, sbc(1,1,isscfc11)
     &,       icfc11)
# endif
# if defined uvic_cfc12
            if (isscfc12 .ne. 0) call getst (j, sbc(1,1,isscfc12)
     &,       icfc12)
# endif
        enddo
      endif

#endif
!-----------------------------------------------------------------------
!     check for consistency in the S.B.C. setup
!-----------------------------------------------------------------------

      call chkcpl

#if defined uvic_global_sums
!-----------------------------------------------------------------------
!     get global sums at the start of the run
!-----------------------------------------------------------------------

      dtoih = 0.
      call globalsum (1)

#endif
!-----------------------------------------------------------------------
!     S T A R T    S E G M E N T    L O O P
!-----------------------------------------------------------------------

      do n=1,numseg

!-----------------------------------------------------------------------
!       get the atmospheric S.B.C.
!-----------------------------------------------------------------------

        call gasbc (1, imt, 1, jmt)

!-----------------------------------------------------------------------
!       call the atmospheric model once for each time step until one
!       segment of "segtim" days is complete. hold atmos S.B.C. fixed
!       during each segment and predict average S.B.C. for ocean
!-----------------------------------------------------------------------

        do loop=1,ntspas
#if defined uvic_embm
          call embm (1, imt, 1, jmt)

# if !defined uvic_mom
#  if defined uvic_mtlm
          call mtlmout (1, imt, 1, jmt)
#  endif
          if (tsits .and. iotsi .ne. stdout .and. iotsi .gt. 0) then
            write (*,'(1x, a3, i7, 1x, a32)') 'ts=',itt, stamp
          endif
# endif
#else
          call atmos
#endif
        enddo

#if defined uvic_mtlm
!-----------------------------------------------------------------------
!       get land S.B.C.s
!-----------------------------------------------------------------------

        call glsbc

!----------------------------------------------------------------------
!       call the land-surface and vegetation  model once for each time
!       step until one segment of "segtim" days is complete.
!       hold land S.B.C. fixed
!       during each segment and predict average S.B.C. for the EMBM
!-----------------------------------------------------------------------

        do loop=1,ntspls
           call mtlm
        enddo

#endif
!-----------------------------------------------------------------------
!       get ocean S.B.C.s
!-----------------------------------------------------------------------

        call gosbc

#if defined uvic_mom
!-----------------------------------------------------------------------
!       call the ocean model once for each time step until one
!       segment of "segtim" days is complete. hold ocean S.B.C. fixed
!       during each segment and predict average S.B.C. for atmos
!-----------------------------------------------------------------------

        do loop=1,ntspos
          call mom
# if defined uvic_embm
          call embmout (1, imt, 1, jmt)
# endif
# if defined uvic_mtlm
          call mtlmout (1, imt, 1, jmt)
# endif
          if (tsits .and. iotsi .ne. stdout .and. iotsi .gt. 0) then
            write (*,'(1x, a3, i7, 1x, a32)') 'ts=',itt, stamp
          endif
        enddo

#endif
#if defined uvic_global_sums
!-----------------------------------------------------------------------
!       write change in global sums for heat and fresh water
!-----------------------------------------------------------------------

        if (tsits) call globalsum (2)

#endif
      enddo

#if defined uvic_global_sums
!-----------------------------------------------------------------------
!     get global sums at the end of the run
!-----------------------------------------------------------------------

      call globalsum (3)

#endif
!-----------------------------------------------------------------------
!     E N D    S E G M E N T    L O O P
!-----------------------------------------------------------------------

      print*, ' ==>  UVIC_ESCM integration is complete.'

      call release_all

      stop
      end

      subroutine chkcpl
      logical errorc
#include "param.h"
#include "csbc.h"
#include "switch.h"

!-----------------------------------------------------------------------
!     do consistency checks before allowing model to continue
!-----------------------------------------------------------------------

      errorc = .false.
      write (stdout,*) ' '
      write (stdout,*) '    (checking S.B.C. setup)'

      if (dtatm .eq. c0) then
          write (stdout,9000)
     & '==> Error: the atmospheric time step must be set in "setatm"  '
          errorc = .true.
          dtatm = 1.e-6
      endif
!      critv = 1.e-6
      critv = 1.e-4
      if (segtim .ne. c0) then
        r1 = rundays/segtim
      else
        r1 = 0.5
      endif
      r2 = segtim/(dtocn*secday)
      r3 = segtim/(dtatm*secday)
      if (segtim .eq. c0) then
          write (stdout,9000)
     & '==> Error: coupling period "segtim" must be specified         '
          errorc = .true.
      elseif (abs(r1-nint(r1)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of segments       '
     &,'    "segtim"  within "rundays" (the length of the run)        '
          errorc = .true.
      elseif (abs(r2-nint(r2)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of density time   '
     &,'    steps "dtocn"  within "segtim" (the segment time)         '
          errorc = .true.
      elseif (abs(r3-nint(r3)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of atmos time     '
     &,'    steps "dtatm"  within "segtim" (the segment time)         '
          errorc = .true.
      endif
#if defined uvic_mom && defined uvic_embm && defined restorst
        write (stdout,9000)
     & '==> Warning: restoring to surface tracers (restorst)          '
     &,'    with coupled model                                        '
#endif
#if defined uvic_mom && defined uvic_embm && defined uvic_replacst
        write (stdout,9000)
     & '==> Warning: replacing surface tracers (uvic_replacst)        '
     &,'    with coupled model                                        '
#endif
#if defined restorst && defined uvic_replacst
        write (stdout,9000)
     & '==> Warning: both restoring and replacing surface tracers     '
     &,'    (restorst, uvic_replacst ). Replacing takes precidence    '
#endif
      write (stdout,*) '    (End of S.B.C. checks) '
      write (stdout,*) ' '
      if (errorc) stop '=>chkcpl'

9000  format (/,(1x,a80))
      return
      end

      subroutine set (index, num, name, text, inc)

!-----------------------------------------------------------------------
!     increment counter, set index and text
!-----------------------------------------------------------------------

      character(*) :: name, text

      name = text
      index = num
      inc = index + 1

      return
      end

      subroutine getst (jrow, ocnout, ntabc)

#if defined uvic_mom
!-----------------------------------------------------------------------
!     read surface tracers from disk row "jrow"
!-----------------------------------------------------------------------

# include "param.h"
# include "iounit.h"
# include "mw.h"
# include "tmngr.h"

      dimension ocnout(imt,jmt)

      call getrow (latdisk(taup1disk), nslab, jrow
     &,          u(1,1,jmw,1,taup1), t(1,1,jmw,1,taup1))
      do i=1,imt
        ocnout(i,jrow) = t(i,1,jmw,ntabc,taup1)
      enddo
#endif

      return
      end
