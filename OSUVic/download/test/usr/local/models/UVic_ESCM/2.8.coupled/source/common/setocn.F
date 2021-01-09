      subroutine setocn

!-----------------------------------------------------------------------
!     set up everything which must be done only once per run
!-----------------------------------------------------------------------
#if defined uvic_embm && !defined uvic_mom

      implicit none

# include "param.h"
# include "accel.h"
# include "calendar.h"
# include "coord.h"
# include "csbc.h"
# include "cembm.h"
# include "cprnts.h"
# include "emode.h"
# include "grdvar.h"
# include "index.h"
# include "iounit.h"
# include "levind.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"

      integer i, iotraj, ioun, j, n, num_processors
      integer isot1, ieot1, isot2, ieot2, jsot, jeot, ksot, keot, mrot

      logical annlev, annlevobc, initpt

      real ah, ahbkg, ahbi, aidif, am, ambi, cflate, cflats, cfldpe
      real cfldps, cflone, cflons, dtatms, kappa_h, kappa_m, maxcfl
      real runstep, senep, snapde, snaple, snapls, spnep
#else
# include "param.h"
# include "accel.h"
# include "calendar.h"
# if defined sponges
#  include "sponge.h"
# else
      logical annlev
# endif
# if defined neptune
#  include "cnep.h"
# endif
# if defined obc
#  include "cobc.h"
#  if defined orlanski
      logical annlevobc
#  else
#   include "obc_data.h"
#  endif
# else
      logical annlevobc
# endif
# include "coord.h"
# if defined fourfil || defined firfil
#  include "cpolar.h"
# endif
# include "cprnts.h"
# include "cregin.h"
# include "csbc.h"
# if defined shortwave
#  include "cshort.h"
# endif
# include "csnap.h"
# include "ctmb.h"
# include "diag.h"
# include "docnam.h"
# if defined uvic_embm
#  include "cembm.h"
# endif
# include "emode.h"
# include "grdvar.h"
# include "hmixc.h"
# include "index.h"
# include "iounit.h"
# include "isleperim.h"
# if defined isopycmix || defined isoneutralmix
#  include "isopyc.h"
# endif
# include "levind.h"
# include "mw.h"
# if defined trajectories
#  include "ptraj.h"
# else
      logical initpt
# endif
# include "scalar.h"
# include "stab.h"
# include "state.h"
# include "switch.h"
# include "task_on.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vmixc.h"
# if defined uvic_tbt
#  include "tbt.h"
# endif
# if defined fourfil || defined firfil
      dimension kmz(imt,jmt)
# endif
# if defined uvic_tidal_kv
#  include "tidal_kv.h"
# endif
#endif
#if defined uvic_fwa
# include "fwa.h"
#else
      integer isfwa1, iefwa1, isfwa2, iefwa2, jsfwa, jefwa, mrfwa
      logical compensate
      real fwaflxi, fwayri, fwayrf, fwarate
#endif
#if defined uvic_npzd
# include "npzd.h"
#else
      real alpha, kw, kc, abio, bbio, cbio, k1, nup, gamma1, gbio
      real epsbio, nuz, gamma2, nud0, LFe, wd0, par, dtnpzd
      real redctn, ki, redptn, capr, dcaco3, nupt0, redotn
      real jdiar
#endif

      character (120) :: fname, new_file_name, logfile

      integer ib(10), ic(10)

      logical error, vmixset, hmixset, uvic_mk, exists

      namelist /contrl/ init, runlen, rununits, restrt, initpt
     &,                 num_processors, runstep
      namelist /tsteps/ dtts, dtuv, dtsf, dtatm, dtatms, namix, segtim
     &,                 daylen
      namelist /riglid/ mxscan, sor, tolrsf, tolrsp, tolrfs
      namelist /mixing/ am, ah, ahbkg, ambi, ahbi, kappa_m, kappa_h
     &,                 cdbot, spnep, senep, aidif, ncon, nmix, eb
     &,                 acor, dampts, dampdz, annlev, annlevobc
      namelist /diagn/  tsiint, tsiper, tavgint, itavg, tmbint, tmbper
     &,                 itmb, stabint, zmbcint, glenint, trmbint, itrmb
     &,                 vmsfint, gyreint,igyre, extint, prxzint, trajint
     &,                 exconvint, dspint, dspper, snapint, snapls
     &,                 snaple, snapde, timavgint, timavgper, cmixint
     &,                 prlat, prslon, prelon, prsdpt, predpt
     &,                 slatxy, elatxy, slonxy, elonxy
     &,                 cflons, cflone, cflats, cflate, cfldps, cfldpe
     &,                 maxcfl, xbtint, xbtper, crossint, densityint
     &,                 fctint, tyzint, restint, tbtint, tbtper
      namelist /io/     expnam, logfile, runstamp, iotavg, iotmb, iotrmb
     &,                 ioglen, iovmsf, iogyre, ioprxz, ioext, iodsp
     &,                 iotsi, iozmbc, iotraj, ioxbt, uvic_mk, isot1
     &,                 ieot1, isot2, ieot2, jsot, jeot, ksot, keot
     &,                 mrot
      namelist /ictime/ year0, month0, day0, hour0, min0, sec0, ryear
     &,                 rmonth, rday, rhour, rmin, rsec, refrun, refinit
     &,                 refuser, eqyear, eqmon, monlen, init_time
     &,                 init_time_in, init_time_out, omega
      namelist /npzd/   alpha, kw, kc, abio, bbio, cbio, k1, nup, gamma1
     &,                 gbio, epsbio, nuz, gamma2, nud0, LFe, wd0, par
     &,                 dtnpzd, redctn, ki, redptn, capr, dcaco3
     &,                 nupt0, redotn, jdiar
      namelist /fwa/    isfwa1, iefwa1, isfwa2, iefwa2, jsfwa, jefwa
     &,                 mrfwa, fwaflxi, fwayri, fwayrf, fwarate
     &,                 compensate

!-----------------------------------------------------------------------
!     physical constants
!-----------------------------------------------------------------------

      rho0 = 1.035
      rho0r = c1/rho0
      grav = 980.6
      radius = 6370.0e5
      pi = atan(1.0) * 4.0

!-----------------------------------------------------------------------
!     initialize variables
!-----------------------------------------------------------------------

!     set defaults for namelist contrl
      init = .true.
      runlen = 365.0
      rununits = 'days'
      restrt = .true.
      initpt =.true.
      num_processors = 1
      runstep = -1.0

!     set defaults for namelist tsteps
      dtts = 43200.0
      dtuv = 600.0
      dtsf = 600.0
      dtatm = 43200.0
      dtatms = 43200.0
      namix = 16
      segtim = 1.0
      daylen = 86400.0

!     set defaults for namelist riglid
      mxscan = 300
      sor = 1.60
      tolrsf = 5.0e8
      tolrsp = 1.0e4
      tolrfs = 1.0e4

      am = 2.0e9
      ah = 2.0e7
      ahbkg = 0.
      ambi = 1.0e23
      ahbi = 5.0e22
      kappa_m = 10.0
      kappa_h = 1.0
      cdbot = 1.3e-3
      spnep = 3.0e5
      senep = 12.0e5
#if defined implicitvmix || defined isopycmix || defined isoneutralmix || defined redi_diffusion
      aidif = 1.0
#else
      aidif = 0.0
#endif
      ncon = 1
      nmix = 16
      eb = .false.
      acor = 0.0
      do n=1,nt
        dampts(n) = 50.0
        dampdz(n) = 26.575e2
      enddo
      annlev = .false.
      annlevobc = .false.

!     set defaults for namelist mixing
      tsiint = 1.0
      tsiper = 1.0
      tavgint = -36500.0
      itavg = .true.
      tmbint = -36500.0
      tmbper = 365.0
      itmb = .true.
      stabint = -36500.0
      zmbcint = -36500.0
      glenint = -36500.0
      trmbint = -36500.0
      itrmb = .true.
      vmsfint = -36500.0
      gyreint = -36500.0
      igyre = .true.
      extint = -36500.0
      prxzint = -36500.0
      trajint = -36500.0
      exconvint = -36500.0
      dspint = -36500.0
      dspper = -365.0
      snapint = 36500.0
      snapls = -90.0
      snaple = 90.0
      snapde = 5000.0e2
      timavgint = 36500.0
      timavgper = 365.0
      cmixint = -36500.0
      do n=1, nlatpr
        prlat(n)  = 100.0
        prslon(n) = 0.0
        prelon(n) = 0.0
        prsdpt(n) = 0.0
        predpt(n) = 6000.0e2
        if (n. le. 4) then
          prslon(n) = 180.0
          prelon(n) = 250.0
        endif
      enddo
      prlat(1) = -60.0
      prlat(2) =   0.0
      prlat(3) =  27.0
      prlat(4) =  55.0
      slatxy   =  -90.0
      elatxy   =   90.0
      slonxy   =   3.0
      elonxy   = 357.0
      cflons = 0.0
      cflone = 360.0
      cflats = -90.0
      cflate = 90.0
      cfldps = 0.0
      cfldpe = 6000.0e2
      maxcfl = 3
      xbtint = -36500.0
      xbtper = -365.0
      crossint = 365000.0
      densityint = -36500.0
      fctint = -36500.0
      tyzint = -36500.0
      restint = 36500.0
      tbtint = -36500.0
      tbtper = -365.0

!     set defaults for namelist io
      expnam = ' '
      logfile = 'machine.log'
      runstamp = ' '
      restrt = .false.
      iotavg = -1
      iotmb = -1
      iotrmb = -1
      ioglen = -1
      iovmsf = -1
      iogyre = -1
      ioprxz = -1
      ioext = -1
      iodsp = -1
      iotsi = -1
      iozmbc = -1
      ioxbt = -1
      uvic_mk = .false.
      isot1 = 2
      ieot1 = imtm1
      isot2 = 2
      ieot2 = 1
      jsot = 2
      jeot = jmtm1
      ksot = 1
      keot  = km
      mrot = 0

!     set defaults for namelist ictime
      year0 = 1
      month0 = 1
      day0 = 1
      hour0 = 0
      min0 = 0
      sec0 = 0
      ryear = 1
      rmonth = 1
      rday = 1
      rhour = 0
      rmin = 0
      rsec = 0
      refrun = .false.
      refinit = .true.
      refuser = .false.
      eqyear = .true.
      eqmon = .false.
      monlen = 30
      init_time = .false.
      init_time_in = .false.
      init_time_out = .false.
      omega = pi/43082.0

!     set defaults for namelist npzd
      alpha = 0.1  ! Initial slope P-I curve [(W/m^2)^(-1)/day]
      kw = 0.04  ! Light attenuation due to water [1/m]
      kc = 0.047  ! Light atten. by phytoplankton [1/(m*mmol/m^3)]
      ki = 5.0  ! Light attenuation through sea ice & snow
      abio = 0.11  ! a; Maximum growth rate parameter [1/day]
      bbio = 1.066  ! b
      cbio = 1.0  ! c [1/deg_C]
      k1 = 0.7  ! Half saturation constant for N uptake [mmol/m^3]
      nup = 0.025  ! Specific mortality rate (Phytoplankton) [1/day]
      nupt0 = 0.02  ! Specific mortality rate (Phytoplankton) [1/day]
      gamma1 = 0.925  ! gama1; Assimilation efficiency (zpk)
      gbio = 1.575  ! Maximum grazing rate [1/day]
      epsbio = 1.6  ! Prey capture rate [(mmol/m^3)^(-2)day^(-1)]
      nuz = 0.34  ! Quadratic mortality (zpk) [(mmol/m^3)^(-2)day^(-1)]
      gamma2 = 0.01  ! Excretion (zpk)
      nud0 = 0.048  ! remineralization rate [1/day]
      LFe = 1.0  ! Iron limitation
      wd0 = 6.0  ! Sinking speed of detritus at surface [m/day]
      par = 0.43  ! fraction of photosythetically active radiation
      dtnpzd = dtts  ! time step of biology [s]
      redctn = 7.  ! C/N Redfield ratio
      redptn = 1./16.  ! P/N Redfield ratio
      capr = 0.018  ! carbonate to carbon production ratio
      dcaco3 = 450000.0  ! remineralisation depth of calcite [cm]
      redotn = 10.6  ! O2/N Redfield ratio
      jdiar = 0.5  ! factor reducing the growth rate of diazotrophs

!     set defaults for namelist fwa
      isfwa1 = 1
      iefwa1 = imt
      isfwa2 = 0
      iefwa2 = 0
      jsfwa = 1
      jefwa = jmt
      mrfwa = 1
      fwaflxi = 0.
      fwayri = -1.e20
      fwayrf =  1.e20
      fwarate = 0.
      compensate = .false.
#if defined uvic_mom

!     set other stuff

      do k=1,km
        dtxcel(k) = 1.0
      enddo

# if defined uvic_tidal_kv
      gravrho0r = grav*rho0r
      zetar = 1./500.e2        ! inverse vertical decay scale, zeta = 500m
      ogamma = 0.2*rho0r*zetar ! Osborn, 1980's mixing efficiency, gamma
                               ! scaled by 1./(zeta*rho0)
# endif
# if defined coarse_grained_parallelism
      cmixts    = .false.
      crossts   = .false.
      densityts = .false.
      fctts     = .false.
      exconvts  = .false.
      snapts    = .false.
      timavgts  = .false.
      prxzts    = .false.

# endif
# if defined firfil || defined fourfil
!     set latitudes used in filtering of tracer and velocity fields

      rjfrst = -87.3
      rjft0  = -67.5
      rjft1  = -69.3
      rjft2  =  69.3
      rjfu0  = -68.4
      rjfu1  = -70.2
      rjfu2  =  70.2

# endif
# if defined firfil
      do j=1,jmtfil
        numflt(j) = 1
        numflu(j) = 1
      enddo

# endif
!      error in tracer conservation generated by wide_open_mw = .true.
!      if (jmw .eq. jmt) then
!        wide_open_mw = .true.
!      else
        wide_open_mw = .false.
!      endif

!     stability diagnostic

      call stabi

      error    = .false.
      vmixset  = .false.
      hmixset  = .false.

      visc_cbu_limit = 1.0e6
      diff_cbt_limit = 1.0e6
# if defined rigid_lid_surface_pressure

      alph     = 1.0
      gam      = 0.0
      theta    = 1.0
# endif
# if defined implicit_free_surface

      alph     = 0.3333333
      gam      = 0.3333333
      theta    = 0.5
# endif
#endif

!-----------------------------------------------------------------------
!     provide for change in above presets using "namelist"
!-----------------------------------------------------------------------

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, contrl, end=101)
101   continue
      runstep = float(int(runstep))
      write (stdout,contrl)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, tsteps, end=102)
102   continue
      write (stdout,tsteps)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, riglid, end=103)
103   continue
      write (stdout,riglid)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, mixing, end=104)
104   continue
      write (stdout,mixing)
      call relunit (ioun)
# if defined isopycmix || defined isoneutralmix
      ah = ahbkg
# endif
      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, diagn, end=105)
105   continue
      write (stdout,diagn)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, io, end=106)
106   continue
      write (stdout,io)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, ictime, end=107)
107   continue
      write (stdout,ictime)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, npzd, end=108)
108   continue
      write (stdout,npzd)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, fwa, end=109)
109   continue
      write (stdout,fwa)
      call relunit (ioun)

      dtocn = dtts

!     limit some variables
      isot1 = max(isot1, 1)
      ieot1 = min(ieot1, imt)
      isot2 = max(isot1, 1)
      ieot2 = min(ieot1, imt)
      jsot = max(jsot, 1)
      jeot = min(jeot, jmt)
      ksot = max(ksot, 1)
      keot = min(keot, km)
      tsiper = min(tsiper, tsiint)
      tbtper = min(tbtper, tbtint)
      if (runstep .gt. 0.0) runlen = min(runstep, runlen)
      if (init_time) then
        init_time_in = .true.
        init_time_out = .true.
      endif

      if (eqyear) then
        if (eqmon) then
          write (timunit,'(a,i10,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')
     &      'equal_month_year since', year0, '-', month0, '-', day0, ' '
     &,      hour0, ':', min0, ':', sec0, '.0'
         else
          write (timunit,'(a,i10,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')
     &      'common_year since' ,year0, '-', month0, '-', day0, ' '
     &,     hour0, ':', min0, ':', sec0, '.0'
        endif
      else
        write (timunit,'(a,i10,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')
     &    'Gregorian_year since', year0, '-', month0, '-', day0, ' '
     &,   hour0, ':', min0, ':', sec0, '.0'
      endif

!-----------------------------------------------------------------------
!     set runstamp
!-----------------------------------------------------------------------
      if (runstamp .eq. ' ') then
!       read runstamp from last line of logfile if it exists
        fname = new_file_name (logfile)
        inquire (file=trim(fname), exist=exists)
        if (exists) then
          call getunit (ioun, fname, 'f s r')
          do while (exists)
            read  (ioun, '(a120)', end=120) runstamp
          enddo
120       continue
          call relunit (ioun)
          print*, " "
          print*, "runstamp: ", trim(runstamp)
        endif
      endif

#if defined uvic_embm && !defined uvic_mom
!-----------------------------------------------------------------------
!     finish set up for a run without mom
!-----------------------------------------------------------------------

      call grids

!-----------------------------------------------------------------------
!     read land mask kmt for run without mom
!-----------------------------------------------------------------------

      call topog (kmt, map, xt, yt, zt, xu, yu, zw, imt, jmt, km)

!     construct depth arrays associated with "u" cells

      do j=1,jmt
        kmu(imt,j) = 0
      enddo
      do i=1,imt
        kmu(i,jmt) = 0
      enddo
      do j=1,jmt-1
        do i=1,imt-1
           kmu(i,j) = min (kmt(i,j), kmt(i+1,j), kmt(i,j+1)
     &,                    kmt(i+1,j+1))
        enddo
      enddo
# if defined cyclic
      do j=1,jmt
        kmu(imt,j) = kmu(2,j)
      enddo
# endif
# if defined symmetry
      do i=1,imt
        kmu(i,jmt) = kmu(i,jmt-2)
      enddo
# endif

      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,               imt, jmt, km, mnisle, maxipp
     &,               xu, yu, zw)

      return
      end
#else
!     user specified tracer names are placed into "trname" here.

      do m=1,nt
        trname(m) = '**unknown***'
      enddo
      trname(1) = 'potentl_temp'
      trname(2) = 'salinity    '

!-----------------------------------------------------------------------
!     get i/o units needed for prognostic variables
!-----------------------------------------------------------------------

      call getunitnumber (kflds)
      call getunitnumber (latdisk(1))
      call getunitnumber (latdisk(2))
# if defined coarse_grained_parallelism
      call getunitnumber (latdisk(3))

!     specify starting and ending latitude rows for each task

      trows = float(jmt-2)/num_processors
      write (stdout,'(/,10x,a)')
     & ' Assigning tasks as follows for coarse_grained_parallelism'
      do n = 1,num_processors
        jstask(n) = max(int((n-1)*trows - jextra + 1.0001),  1 - jextra)
        jetask(n) = min(int(n*trows + 1 + jextra + 1.0001),jmt + jextra)
        num_rcpt(n) = jetask(n) - jstask(n) + 1 - (jmw - ncrows)
        write (stdout,'(a,i4,a,i4,a,i4,a,i4)')
     & 'Processor # ',n, ': starting row =', jstask(n)
     &,': ending row =', jetask(n)
     &,': calculated rows per task =', num_rcpt(n)
      enddo
      write (stdout,*) ' '
# endif
# if defined uvic_tidal_kv

!-----------------------------------------------------------------------
!     read in 2d array of tidal energy dissipation rate
!     similar to that of Jayne and St. Laurent, GRL 2001
!     except here we use E0 = W/(A*Nlev)

!     no need to convert MOM native cgs
!     units lesson: W = kg m^2 /s^3, therefore, W/m^2 = kg/s^3
!     therefore multiply edr by N0*1e3 to get g/s^3, where N0 is some
!     reference value of N.
!     in tidal.nc the energy flux already includes the Nb, so it's only
!     multiplied by 1e3 to convert W/m^2 in g/s^3
!-----------------------------------------------------------------------

      edr(:,:)=0.
      ib(:) = 1
      ic(:) = imt
      ic(2) = jmt
      fname = new_file_name ("tidal.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        c1e3 = 1.e3
        call openfile (trim(fname), iou)
        call getvara ('dissipation', iou, imt*jmt, ib, ic, edr
     &,   c1e3, c0)
        call closefile (iou)
      else
        print*, "Warning => Can not find ",trim(fname)
      endif
# endif

!-----------------------------------------------------------------------
!     set up the grids in x (longitude), y (latitude), and z (depth)
!     corresponding to Arakawa "b" gird system
!-----------------------------------------------------------------------

      call grids

!-----------------------------------------------------------------------
!     compute density coefficients based on depth of grid points
!-----------------------------------------------------------------------

      call eqstate (zt, km, ro0, to, so, c, tmink, tmaxk, smink, smaxk)

      cksum = checksum(c, km, 9)
      write (stdout
     &,'(6x,"Checksum for density coefficients =",e14.7)') cksum

# if defined linearized_density
!     eliminate the pressure effect when using linearized density

      do k=1,km
        ro0(k) = c0
        to(k)  = c0
        so(k)  = c0
      enddo

# endif
# if defined linearized_advection || defined equatorial_thermocline
!     set the initial idealized stratification as a function of
!     temperature only. Salinity is fixed at 35 ppt (tbar(k,2)=0) but
!     passive tracers (if nt>2) are 1.0 at k=1 and 0.0 for k>1

      tzero = 7.5
      tone  = 10.0
      z0    = 30.0e2
      bigl  = 80.0e2
      write (stdout,'(//a/)')
     & 'Initial Tracer profile as a function of depth'
      do k=1,km
        tbarz(k,1) = tzero*(1.0-tanh((zt(k)-bigl)/z0)) +
     &                tone*(1.0-zt(k)/zt(km))
        tbarz(k,2) = c0
        if (nt .gt. 2) then
          do n=3,nt
            if (k .eq. 1) then
              tbarz(k,n) = c1
            else
              tbarz(k,n) = c0
            endif
          enddo
        endif
        write (stdout,'(1x,"k=",i3,", zt(k)=",f8.1
     &,        "cm,  tbarz(k,n),n=1,nt =",10e14.7)')
     &      k, zt(k),(tbarz(k,n),n=1,nt)
      enddo

# endif
!-----------------------------------------------------------------------
!     if the MW is not fully opened, then set time level indicators in
!     the MW ("tau-1" "tau" "tau+1") to constant values.
!-----------------------------------------------------------------------

      if (.not. wide_open_mw) then
        taum1 = -1
        tau   =  0
        taup1 = +1
      endif

!-----------------------------------------------------------------------
!     set prognostic quantities to either initial conditions or restart
!-----------------------------------------------------------------------

!     generate topographic levels "kmt" on "t" cells.

      call topog (kmt, map, xt, yt, zt, xu, yu, zw, imt, jmt, km)

!     construct depth arrays associated with "u" cells

      call depth_u (kmt, imt, jmt, zw, km, kmu, h, hr)
# if defined neptune

!     calculate neptune velocities

      call neptune
# endif

!     initialize two dimensional fields on disk

# if defined rigid_lid_surface_pressure || defined implicit_free_surface

!     initialize surface pressure fields in memory

      do jrow=1,jmt
        do i=1,imt
          ps(i,jrow,1)     = c0
          ps(i,jrow,2)     = c0
          pguess(i,jrow)   = c0
          ubar(i,jrow,1)   = c0
          ubar(i,jrow,2)   = c0
          ubarm1(i,jrow,1) = c0
          ubarm1(i,jrow,2) = c0
        enddo
      enddo
# endif
# if defined stream_function

!     initialize stream function fields in memory

      do n=1,2
        do jrow=1,jmt
          do i=1,imt
# if defined neptune

!           initialize to neptune streamfunction

            psi(i,jrow,n) = pnep(i,jrow)
            ptd(i,jrow)   = c0
# else
            psi(i,jrow,n) = c0
            ptd(i,jrow)   = c0
# endif
          enddo
        enddo
      enddo

!     initialize stream function guess fields on disk
!     block`s 1 & 2 are for the stream function guess field on disk

      do n=1,nkflds
        call oput (kflds, nwds, n, ptd(1,1))
      enddo
# endif

!     initialize all latitude rows

      call rowi

!     initialize time step counter = 0

      itt    = 0
      irstdy = 0
      msrsdy = 0

!     initialize all prognostic quantities from the restart

      if (.not. init) then
        fname = new_file_name ("restart_mom.nc")
        call mom_rest_in (fname, 1, imt, 1, jmt)
# if defined neptune

!       calculate neptune velocities

        call neptune
# endif

!       compute a topography checksum

        cksum = 0.0
        do jrow=1,jmt
          do i=1,imt
            cksum = cksum + i*jrow*kmt(i,jrow)
          enddo
        enddo
        write (stdout,*) ' "kmt" checksum = ', cksum
      endif

# if defined time_averages
!-----------------------------------------------------------------------
!     zero integrated time average accumulators
!-----------------------------------------------------------------------

      call ta_mom_tsi (0)

# endif

!-----------------------------------------------------------------------
!     initialize the time manager with specified initial conditions
!     time, user reference time, model time, and how long to integrate.
!-----------------------------------------------------------------------

      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,              ryear, rmonth, rday, rhour, rmin, rsec
     &,              irstdy, msrsdy
     &,              runlen, rununits, rundays, dtts)

# if defined stability_tests

!-----------------------------------------------------------------------
!     convert starting and ending longitudes for the stability tests
!     to nearest model grid points.
!-----------------------------------------------------------------------

      if (stabint .ge. c0) then
        iscfl  = max(indp (cflons, xt, imt), 2)
        cflons = xt(iscfl)
        iecfl  = min(indp (cflone, xt, imt), imt-1)
        cflone = xt(iecfl)
        jscfl  = max(indp (cflats, yt, jmt), 2)
        cflats = yt(jscfl)
        jecfl  = min(indp (cflate, yt, jmt), jmt-1)
        cflate = yt(jecfl)
        kscfl  = indp (cfldps, zt, km)
        cfldps = zt(kscfl)
        kecfl  = indp (cfldpe, zt, km)
        cfldpe = zt(kecfl)
      endif
# endif

# if defined restorst && !defined uvic_replacst

!-----------------------------------------------------------------------
!      damp surface tracers back to data. schematically, the restoring
!      term will be = (dampdz/(dampts*daylen)) * (data - t)
!      where dampdz is the thickness (cm) and dampts is the time
!      scale (days)
!-----------------------------------------------------------------------

      do n=1,nt
        write (stdout,'(/,1x,a,i2,a,a,/,1x,1pe14.7,a,1pe14.7,a,/)')
     &  ' Surface tracer #',n,' will be damped back to data using a'
     &, ' Newtonian restoring time scale of '
     &,  dampts(n),' days. and a level thickness =', dampdz(n),' cm.'
      enddo
# endif

# if defined fourfil || defined firfil

!     compute sin and cos values for vector correction before filtering

      fx =  1.0e-10
      fxa = dxt(1)/radius

      do i=2,imtm1
        fxb = fxa*float(i-2)
        spsin(i) = sin(fxb)
        spcos(i) = cos(fxb)
        if (abs(spsin(i)) .lt. fx) spsin(i) = c0
        if (abs(spcos(i)) .lt. fx) spcos(i) = c0
      enddo

      spsin(1) = c0
      spcos(1) = c0
      spsin(imt) = c0
      spcos(imt) = c0

!     set up model indices for filtering high latitudes

      jfrst  = indp (rjfrst, yt, jmt)
      jft0   = indp (rjft0, yt, jmt)
      jft1   = indp (rjft1, yt, jmt)
      jft2   = indp (rjft2, yt, jmt)
      jfu0   = indp (rjfu0, yu, jmt)
      jfu1   = indp (rjfu1, yu, jmt)
      jfu2   = indp (rjfu2, yu, jmt)
      jskpt  = jft2-jft1
      jskpu  = jfu2-jfu1
      njtbft = (jft1-jfrst+1)+(jmtm1-jft2+1)
      njtbfu = (jfu1-jfrst+1)+(jmtm1-jfu2+1)
      if (njtbft .gt. jmtfil .or. njtbfu .gt. jmtfil) then
        write (stdout,9599) njtbft, njtbfu
        write (stderr,9599) njtbft, njtbfu
        stop '>ocean 1'
      endif
9551  format (/' ==== start and end indices for Fourier filtering ====')
9552  format (' only 11 sets of indices fit across the page.',
     &       '  others will not be printed.'/)
9553  format (/,' == filtering indices for t,s ==')
9554  format (/,' == filtering indices for u,v ==')
9555  format (/,' == filtering indices for stream function ==')
9599  format (/,' error => jmtfil must be >= max(njtbft,njtbfu)',
     &        /,'          njtbft=',i8,' njtbfu=',i8)

#  if defined firfil

!     set "numflt" and "numflu" to filter more at higher latitudes
!     using "jft0" and "jfu0" as reference latitude rows.

      numfmx = imt * p25
      refcos = cst(jft0)
      write(stdout,9501) refcos, yt(jft0)
      do jrow=jfrst,jmt
        if ((jrow .le. jft1 .or. jrow .ge. jft2) .and. jrow .ge. jfrst)
     &    then
          jj = jrow - jfrst + 1
          if (jrow .ge. jft2) jj = jj - jskpt + 1
          numflt(jj) = max(1,int(refcos/cst(jrow)))
          if (numflt(jj) .gt. numfmx) numflt(jj) = numfmx
          write(stdout,9502) jrow, jj, numflt(jj), yt(jrow)
          if (jj .gt. jmtfil) then
            write (stdout,'(1x,a,i4,a,i4,a)')
     &      'Error: jj exceeds jmtfil. jj=',jj, ' jmtfil=',jmtfil
     &,     '.  increase jmtfil'
            stop "=>setocn"
          endif
        endif
      enddo
      refcos = csu(jfu0)
      write(stdout,9503) refcos, yu(jfu0)
      do jrow=jfrst,jmt
        if ((jrow .le. jfu1 .or. jrow .ge. jfu2) .and. jrow .ge. jfrst)
     &    then
          jj = jrow - jfrst + 1
          if (jrow .ge. jfu2) jj = jj - jskpu + 1
          numflu(jj) = max(1,int(refcos/csu(jrow)))
          if (numflu(jj) .gt. numfmx) numflu(jj) = numfmx
          write(stdout,9502) jrow, jj, numflu(jj), yu(jrow)
          if (jj .gt. jmtfil) then
            write (stdout,'(1x,a,i4,a,i4,a)')
     &      'Error: jj exceeds jmtfil. jj=',jj, ' jmtfil=',jmtfil
     &,     '.  increase jmtfil'
            stop "=>setocn"
          endif
        endif
      enddo

9501  format(/' firfil reference cosine for tracers =',e12.6,' (',
     &       f7.3,' deg)'/'     jrow      jj numflt(jj)  latitude')
9502  format(1x,3i8,6x,f7.3)
9503  format(/' firfil reference cosine for velocities =',e12.6,' (',
     &       f7.3,' deg)'/'     jrow      jj numflu(jj)  latitude')
#  endif
# endif

# if defined sponges && !defined equatorial_thermocline

!-----------------------------------------------------------------------
!     set latitude functions for Newtonian damping term in sponge layers
!     near artificial northern and southern boundaries in limited
!     domain models.
!     schematically: damping = - sponge(j) * (t(i,k,j) - levitus(i,k,j))
!     all related data is assumed to have been prepared using the
!     "sponge" routines in the MOM dataset.
!     disk resource factor of 13 is for 12 months + 1 annual mean
!-----------------------------------------------------------------------

      lrec = 4 + 4 + 2*jmt + 4*imt*km
      write (opt_sponge,'(a,1x,i8)') 'unformatted direct words =',lrec
      write (stdout,'(/a,1pg10.3,a)')
     & ' Sequential access disk resource for file "sponge.mom" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew1, 'sponge.mom'
     &,             'unformatted sequential rewind ieee')
      write (stdout,'(/a,1pg10.3,a)')
     & ' Direct access disk resource for file "sponges" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew2, 'sponges', opt_sponge)

      sum = 0.0
      cksum = 0.0
      do m=1,13
        read (ionew1)
        read (ionew1) sstamp, spdpm, im, kk, jm, j1, j2, mm
     &, spngs, spngn, spbuf1
        write (ionew2, rec=m) sstamp, spdpm, im, kk, jm, spngs, spngn
     &, spbuf1
        if (m .le. 12) then
          spgdpm(m) = spdpm
          sum = sum + spdpm
          tspng(m) = sum - 0.5*spdpm
        endif
        do n=1,4
          cksum = cksum + checksum (spbuf1(1,1,n), imt, km)
        enddo
      enddo
      print *,' checksum for sponge data = ',cksum
      if (annlev) then
        write (stdout,'(/a)')
     &  ' => Annual mean levitus data will be used for Newtonian sponge'
        do m=1,12
          write (ionew2, rec=m) sstamp, spdpm, im, kk, jm, spngs, spngn
     &,   spbuf1
        enddo
      else
        write (stdout,'(/a)')
     &  ' => Monthly levitus data will be used for Newtonian sponge'
      endif
      write (stdout,'(1a/)')  ' Newtonian damping zone setup:'
      do jrow=1,jmt
        if (spngs(jrow) .ne. c0) then
          write (stdout,'(a,i3,a,f7.2,a,e10.3,a,e10.3)') ' jrow='
     &,   jrow,', yt(jrow)=',yt(jrow), ', spngs(jrow) =',spngs(jrow)
     &,   ', Newtonian time scale (days)=',secday/spngs(jrow)
        endif
        if (spngn(jrow) .ne. c0) then
          write (stdout,'(a,i3,a,f7.2,a,e10.3,a,e10.3)') ' jrow='
     &,   jrow,', yt(jrow)=',yt(jrow), ', spngn(jrow) =',spngn(jrow)
     &,   ', Newtonian time scale (days)=',secday/spngn(jrow)
        endif
        if (spngn(jrow) .ne. c0 .and. spngs(jrow) .ne. c0) then
          write (stdout,'(/a/)')
     &    ' Error: Overlapping north and south sponges not allowed'
          stop '=>setocn'
        endif
      enddo
      write (stdout,*) ' '

!     set the current model time in days and initialize interpolation
!     information

      begtim = (realdays(initial) - 1.0) + realdays(imodeltime)

      if (.not.eqyear) then
        write (stdout,*) '=>Warning leap year being used with '
     &, ' climatological sponges?'
      endif
      iprev  = 1
      inext  = 2
      indxsp = 5
      method = 3
      call timeinterp(begtim, indxsp, tspng, spgdpm, 12, .true., method
     &,               inextd, iprevd, wprev, readsp, inext, iprev)

      read (ionew2, rec=iprevd) stprev, spdpmp, im, kk, jm, spngs
     &,    spngn, spbuf1
      do n=1,4
        do k=1,km
          do i=1,imt
            spbuf(i,k,n,iprev) = spbuf1(i,k,n)
          enddo
        enddo
      enddo
      read (ionew2, rec=inextd) stnext, spdpmn, im, kk, jm, spngs
     &,    spngn, spbuf1
      do n=1,4
        do k=1,km
          do i=1,imt
            spbuf(i,k,n,inext) = spbuf1(i,k,n)
          enddo
        enddo
      enddo

      write (stdout,'(2(/a,i3,1x,a,a,i2)/a,g14.7,1x,a,a,g14.7/)')
     & ' reading sponge record ', iprevd, stprev, ' into buffer ', iprev
     &,' reading sponge record ', inextd, stnext, ' into buffer ', inext
     &,' for day =', begtim, stamp, ' weight =',wprev
      write (stdout,'(/a,i2,a,i2/)') 'sponge is dataset index ',indxsp
     &,' for time interpolation using method #',method
      call relunit (ionew1)
      call relunit (ionew2)
# endif

# if defined obc
#  include "setocn_obc.inc"
# endif

# if defined shortwave

!-----------------------------------------------------------------------
!     Solar Shortwave energy penetrates below the ocean surface. Clear
!     water assumes energy partitions between two exponentials as
!     follows:

!     58% of the energy decays with a 35 cm e-folding scale
!     42% of the energy decays with a 23 m e-folding scale

!     if the thickness of the first ocean level "dzt(1)" is 50 meters,
!     then shortwave penetration wouldn't matter. however, for
!     dzt(1) = 10 meters, the effect can be significant. this can be
!     particularly noticeable in the summer hemisphere.

!     Paulson and Simpson (1977 Irradiance measurements in the upper
!                               ocean JPO 7, 952-956)
!     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
!                  A General Circulation Model for Upper Ocean
!                  Simulation (Rosati and Miyakoda JPO vol 18,Nov 1988)
!-----------------------------------------------------------------------

      write (stdout,'(/a/)')
     &' => Shortwave penetration is a double exponential as follows:'
      rpart  = 0.58
      efold1 = 35.0e0
      efold2 = 23.0e2
      rscl1  = 1.0/efold1
      rscl2  = 1.0/efold2

!     note that pen(0) is set 0.0 instead of 1.0 to compensate for the
!     shortwave part of the total surface flux in "stf(i,1)"

      pen(0) = c0

      do k=1,km
        swarg1 = -min(zw(k)*rscl1,70.0)
        swarg2 = -min(zw(k)*rscl2,70.0)
        pen(k) = rpart*exp(swarg1) + (c1-rpart)*exp(swarg2)
        divpen(k) = (pen(k-1) - pen(k))/dzt(k)
        write (stdout,9200) k, zw(k)*1.e-2, pen(k), zt(k)*1.e-2
     &,                     divpen(k)
      enddo
      write (stdout,*) ' '
9200  format (1x,'k=',i3,' zw(k)=',f8.2,'(m)  pen(k)=',e14.7
     &,       ' zt(k)=',f8.2,'(m)  divpen(k)=',e14.7)
# endif

!-----------------------------------------------------------------------
!     compute surface area and volume of ocean ("t" cells and "u" cells)
!     (note that areas are defined at each level)
!-----------------------------------------------------------------------

      do k=1,km
        tcella(k) = c0
        ucella(k) = c0
      enddo
      ocnp   = c0
      tcellv = c0
      ucellv = c0

!     this comment directive turns off autotasking on the YMP
!     for the following loop

!fpp$ noconcur l
      do jrow=2,jmtm1
        do i=2,imtm1
          tarea = cst(jrow)*dxt(i)*dyt(jrow)
          uarea = csu(jrow)*dxu(i)*dyu(jrow)
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) + tarea
            enddo
            tcellv = tcellv + tarea*zw(kmt(i,jrow))
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
          if (kmu(i,jrow) .gt. 0) then
            do k=1,kmu(i,jrow)
              ucella(k) = ucella(k) + uarea
            enddo
            ucellv = ucellv + uarea*zw(kmu(i,jrow))
          endif
        enddo
      enddo

      write (stdout,9341) tcella(1), tcellv, ucella(1), ucellv

!---------------------------------------------------------------------
!     set the horizontal regional masks and names to be used when
!     computing averages on the "t" grid in subroutine "region.F".
!     also set the vertical regional masks and names for use (along
!     with the horizontal ones) in term balance calculations for
!     tracers & momentum. For term balance calculations the number
!     of masks is the product of horizontal & vertical regions.
!---------------------------------------------------------------------

# if defined readrmsk

!    read in horizontal & vertical regional masks ("mskhr" & "mskvr")
!    and names ("hregnm" & "vregnm") on unit iormsk

      call getunit (iormsk, 'region_masks'
     &,             'formatted sequential rewind')
      call reg1st (iormsk, .true., .false., .false., .false., .true.)
      call relunit (iormsk)
# else

!     set up the horizontal regions:
!     specify "mskhr" and "hregnm" values rather than reading them in
!     from a file (arbitrarily defaulted below for 5 zonal bands).
!     note: there must be "nhreg" calls ... one for each horizontal
!     region. The form is:
!     call sethr (region number, starting lon, ending lon, starting lat,
!                 ending lat)
!     where the lon and lat limits are for the edges of the "t" cells
!     "sethr" will fit the region using the nearest model grid points

      do jrow=1,jmt
        do i=1,imt
          mskhr(i,jrow) = 0
        enddo
      enddo

      js = 1
      ej = 0.0
      do n=1, nhreg
        ej = ej + float(jmtm1)/float(nhreg)
        je = min(jmtm1, nint(ej))
        call sethr (n, xu(1), xu(imtm1), yu(js), yu(je))
        js = je
      enddo

      do jrow=1,jmt
        mskhr(1,jrow)   = 0
        mskhr(imt,jrow) = 0
      enddo

!     set up the vertical regions:
!     specify "mskvr" and "vregnm" values. also the starting & ending
!     levels for the vertical regions (arbitrarily defaulted for two
!     vertical regions). Note: there must be "nvreg" calls... one for
!     each vertical region.
!     The form is:
!     call setvr (region number, starting depth, ending depth)
!     where the depth limits are for the edges of the "t" cells
!     "setvr" will fit the region using the nearest model grid points

      ek = 0.0
      do n=1, nvreg
        ek = ek + float(km)/float(nvreg)
        ke = min(km, nint(ek))
        if (n .eq. 1) then
          call setvr (n, 0.0, zw(ke))
        else
          call setvr (n, zw(ks), zw(ke))
        endif
        ks = ke
      enddo

      do k=1,km
        mskvr(k) = 0
      enddo
      do n=1,nvreg
        ks = llvreg(n,1)
        ke = llvreg(n,2)
        do k=ks,ke
          mskvr(k) = n
        enddo
      enddo

!     print out the regional mask to a ascii file

      call getunit (iou, 'region_masks.out', 'f s r')
      call reg1st (iou, .true., .false., .false., .true., .false.)
      call relunit (iou)
# endif

!-----------------------------------------------------------------------
!     compute the surface area and volume of the ocean regions. index 0
!     refers to the sum of all regions.
!     (values used in subroutine region are done in terms of meters,
!     rather than centimetres)
!-----------------------------------------------------------------------

      areag = c0
      volgt  = c0

      do k=1,km
        volgk(k) = c0
      enddo

      do n=0,numreg
        areat(n) = c0
        areau(n) = c0
        volt(n)  = c0
        volu(n)  = c0
      enddo
      do mask=1,nhreg
        areab(mask) = c0
        volbt(mask) = c0
        do k=1,km
          volbk(mask,k) = c0
        enddo
      enddo

      do jrow=2,jmtm1
        do i=2,imtm1
          mask = mskhr(i,jrow)
          kz   = kmt(i,jrow)
          if (kz .gt. 0 .and. mask .gt. 0) then
            boxat = cst(jrow) * dxt(i) * dyt(jrow)
            if (kmu(i,jrow) .ne. 0) then
              boxau = csu(jrow) * dxu(i) * dyu(jrow)
            else
              boxau = c0
            endif
            areat(0) = areat(0) +  boxat
            areau(0) = areau(0) +  boxau
            areab(mask)  = areab(mask) + boxat * 1.e-4
            do k=1,kz
              volbk(mask,k) = volbk(mask,k) + boxat * dzt(k) * 1.e-6
              dvolt   = boxat * dzt(k)
              if (kmu(i,jrow) .ge. k) then
                dvolu   = boxau * dzt(k)
              else
                dvolu   = c0
              endif
              n = nhreg*(mskvr(k)-1) + mask
              if (n .gt. 0) then
                volt(0) = volt(0) +  dvolt
                volu(0) = volu(0) +  dvolu
                volt(n) = volt(n) +  dvolt
                volu(n) = volu(n) +  dvolu
                do nv=1,nvreg
                  ks = llvreg(nv,1)
                  if (k .eq. ks) then
                    areat(n) = areat(n) +  boxat
                    areau(n) = areau(n) +  boxau
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      enddo

      do mask=0,numreg
        if (areat(mask) .eq. c0) then
          rareat(mask) = c0
        else
          rareat(mask) = c1/areat(mask)
        endif

        if (areau(mask) .eq. c0) then
          rareau(mask) = c0
        else
          rareau(mask) = c1/areau(mask)
        endif

        if (volt(mask) .eq. c0) then
          rvolt(mask) = c0
        else
          rvolt(mask) = c1/volt(mask)
        endif

        if (volu(mask) .eq. c0) then
          rvolu(mask) = c0
        else
          rvolu(mask) = c1/volu(mask)
        endif
      enddo

      do mask=1,nhreg
        do k=1,km
          volbt(mask) = volbt(mask) + volbk(mask,k)
          volgk(k) = volgk(k) + volbk(mask,k)
        enddo
        areag = areag + areab(mask)
        volgt = volgt + volbt(mask)
      enddo

#  if defined tracer_averages
      if (iotavg .ne. stdout .or. iotavg .lt. 0) then
        call getunit (iou, 'tracer_avg.dta'
     &,               'unformatted sequential append ieee')
        call reg1st (iou, .false., .true., .true., .false., .false.)
        call relunit (iou)
      endif
#  endif
      if (iotavg .eq. stdout .or. iotavg .lt. 0) then
        call reg1st (stdout, .true., .true., .true., .false., .false.)
      endif

!     compute and print statistics for regions

      sum = c0
      do n=1,numreg
        sum = sum + volt(n)
      enddo
      sum    = 100.0*sum/tcellv
      pctocn = 100.0*ocnp/float((imt-2)*(jmt-2)*km)
      diffa  = 100.0 * (c1 - (tcella(1) - 10000.0*areag)/tcella(1))

      write (stdout,9342) diffa, numreg, sum, pctocn
9342  format ('  the horizontal regional masks cover',f8.3
     &, '% of the total ocean surface area.'/
     &, '  there are ', i6, ' regions over which tracer & '
     &, 'momentum balances will be computed',/,'  accounting for '
     &, f6.2,'% of the total ocean volume.'/
     &, 1x,f6.2,'% of the grid points lie within the ocean.'/)
#
!-----------------------------------------------------------------------
!     find all land mass perimeters for Poisson solvers
!-----------------------------------------------------------------------

      auto_kmt_changes = .false.
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)

!     set mask for land mass perimeters on which to perform calculations
!     imask(-n) = .false.  [no equations ever on dry land mass n]
!     imask(0)  = .true.   [equations at all mid ocean points]
!     imask(n)  = .true./.false [controls whether there will be
!                                equations on the ocean perimeter of
!                                land mass n]
!     note: land mass 1 is the northwest-most land mass
!           usually includes the "north pole", and at low resolutions,
!           the "main continent"
!     for the numbering of the other landmasses, see generated map(i,j)
!           by selecting option -Dshow_details

      do isle=-mnisle,mnisle
        if (isle .ge. 0 .and. isle .le. nisle) then
          imask(isle) = .true.
        else
          imask(isle) = .false.
        endif
      enddo
# if defined symmetry

!     do not perform island integrals for land masses whose perimeters
!     touch the equator in models symmetric about the equator

      do isle=1,nisle
        do n=1,nippts(isle)
          j = jperm(iofs(isle)+n)
          if (j .eq. jmt-1) then
            imask(isle) = .false.
          endif
        enddo
      enddo
# endif

!     user-specified changes to island mask
!       imask(1) = .true.
!       imask(2) = .true.

!     there are problems if imask is set .true. for a nonexistent
!     island.

!     print diagnostic information

      do isle=-mnisle,mnisle
        if (imask(isle)) then
          if (isle .eq. 0) then
            print '(a)','=> calculations enabled for mid ocean points'
          else
            print '(2a,i3)','=> calculations enabled for ocean ',
     &                      'perimeter of land mass',isle
          endif
        endif
      enddo
      do isle=0,nisle
        if (.not. imask(isle)) then
            print '(2a,i3)','=> calculations disabled for ocean ',
     &                      'perimeter of land mass',isle
        endif
      enddo

!     imain is the land mass on which dpsi is normalized to 0
!     if imain is 0, then dpsi is not normalized.
!     default value of imain is land mass with longest perimeter

      imain = min(2,nisle)
      do isle=1,nisle
        if (nippts(isle) .gt. nippts(imain)) then
          imain = isle
        endif
      enddo

!     if any island perimeter is not calculated, imain must be one such

      do isle=1,nisle
        if (.not.imask(isle)) then
          imain = isle
        endif
      enddo
# if defined symmetry

!     do not normalize dpsi when using symmetry about the equator

      imain = 0
# endif

      if (imain .gt. 0 .and. imain .le. nisle) then
        print '(a,i4)', 'dpsi normalized to zero on land mass',imain
      elseif (imain .eq. 0) then
        print *, 'no normalization on dpsi'
      else
        print *, 'ERROR: illegal value for choice of normalization ',
     &           'land mass, imain =', imain
      endif
      print *,' (user may set "imain" to any valid land mass number)'

!---------------------------------------------------------------------
!     compute checksum of density coefficients
!---------------------------------------------------------------------

      print *,' '
      call print_checksum (c(1,1), km, 9
     &,                   ' density coefficient checksum = ')

# if defined stream_function

!-----------------------------------------------------------------------
!     checksum the starting stream function.
!-----------------------------------------------------------------------

      call print_checksum (psi(1,1,1), imt, jmt
     &, ' checksum for psi(,,1) = ')
      call print_checksum (psi(1,1,2), imt, jmt
     &, ' checksum for psi(,,2) = ')
# endif

# if defined fourfil || defined firfil

!-----------------------------------------------------------------------
!     compute an array to indicate "interior" stream function grid cells
!-----------------------------------------------------------------------

      do jrow=1,jmt
        kmz(1,jrow)   = 0
        kmz(imt,jrow) = 0
      enddo

      do i=1,imt
        kmz(i,1)   = 0
        kmz(i,jmt) = 0
      enddo

      do jrow=2,jmtm1
        do i=2,imt
          kmz(i,jrow) = min(kmu(i-1,jrow-1), kmu(i,jrow-1)
     &,                     kmu(i-1,jrow), kmu(i,jrow))
        enddo
      enddo
#  if defined cyclic
      do jrow=1,jmt
        kmz(1,jrow) = kmz(imtm1,jrow)
      enddo
#  endif

!---------------------------------------------------------------------
!     find and print start & end indices for filtering
!---------------------------------------------------------------------

      write (stdout,9551)
      if (lsegf.gt.11) write (stdout,9552)
      write (stdout,9553)
      call findex (kmt, jmtfil, km, jft1, jft2, imt, istf, ietf)
      write (stdout,9554)
      call findex (kmu, jmtfil, km, jfu1, jfu2, imt, isuf, ieuf)
      write (stdout,9555)
#  if defined rigid_lid_surface_pressure || defined implicit_free_surface
      call findex (kmu, jmtfil, 1, jfu1, jfu2, imt, iszf, iezf)
#  endif
#  if defined stream_function
      call findex (kmz, jmtfil, 1, jft1, jft2, imt, iszf, iezf)
#  endif
# endif

!---------------------------------------------------------------------
!     print the timestep multipliers
!---------------------------------------------------------------------

      write (stdout,9601) (dtxcel(k),k=1,km)

!-----------------------------------------------------------------------
!     initialize various things
!-----------------------------------------------------------------------

      do jrow=1,jmt
        do i=1,imt
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
          uhat(i,jrow,1) = c0
          uhat(i,jrow,2) = c0
          divf(i,jrow) = c0
# endif
# if defined stream_function
          ztd(i,jrow) = c0
# endif
          zu(i,jrow,1)  = c0
          zu(i,jrow,2)  = c0
        enddo
      enddo

!     Coriolis factors

      do jrow=1,jmt
        do i=1,imt
          cori(i,jrow,1) = c2*omega*sin(ulat(i,jrow)/radian)
          cori(i,jrow,2) = -cori(i,jrow,1)
        enddo
      enddo

!     metric diffusion factors

# if defined consthmix
#  if defined biharmonic
      amix = sqrt(abs(ambi))
#  else
      amix = am
#  endif
      do jrow=1,jmt
        am3(jrow)   = amix*(c1-tng(jrow)*tng(jrow))/(radius**2)
        am4(jrow,1) = -amix*c2*sine(jrow)/(radius*csu(jrow)
     &                                     *csu(jrow))
        am4(jrow,2) = -am4(jrow,1)
      enddo
# endif

!     metric advection factors

      do jrow=1,jmt
        advmet(jrow,1) = tng(jrow)/radius
        advmet(jrow,2) = -advmet(jrow,1)
      enddo

!     diffusive flux through bottom of cells

      do j=jsmw,jemw
        do k=0,km
          do i=1,imt
            diff_fb(i,k,j) = c0
            adv_fb(i,k,j)  = c0
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     initialize diagnostics
!-----------------------------------------------------------------------

# if defined meridional_tracer_budget

!-----------------------------------------------------------------------
!     set basins and initialize arrays
!-----------------------------------------------------------------------

      if (tmbint .ge. c0) then

!       set all points to signify basin # 1

        do jrow=1,jmt
          do i=2,imtm1
            msktmb(i,jrow) = 1
          enddo
          msktmb(1,jrow)   = 0
          msktmb(imt,jrow) = 0
        enddo

!       divide "msktmb" into other basins (2 .. ntmbb) here if desired.

!       set all land points to basin # 0

        do jrow=1,jmt
          do i=1,imt
            if (kmt(i,jrow) .eq. 0) msktmb(i,jrow) = 0
          enddo
        enddo

        if (ntmbb .gt. 1) then
          write (stdout,*)
     &  ' Basin arrangement for Meridional Tracer Budget diagnostic'
          call iplot (msktmb, imt, imt, jmt)
        endif

!       write out the meridional tracer budget basin mask

        if ((iotmb .ne. stdout .or. iotmb .lt. 0) .and. itmb) then
          call getunit (iu, 'tracer_bud.dta'
     &,                'unformatted sequential append ieee')
          iotext =
     &  ' read (iotmb) imt, jmt, ((msktmb(i,j),i=1,imt),j=1,jmt)'
          write (iu) stamp, iotext, expnam
          write (iu) imt, jmt, ((msktmb(i,j),i=1,imt),j=1,jmt)
          call relunit (iu)
        endif
      endif

      numtmb = 0
      do mask=0,ntmbb
        do jrow=1,jmt
          smdvol(jrow,mask)  = c0
        enddo
      enddo
      do mask=0,ntmbb
        do n=1,nt
          do jrow=1,jmt
            tstor(jrow,n,mask) = c0
            tdiv(jrow,n,mask)  = c0
            tflux(jrow,n,mask) = c0
            tdif(jrow,n,mask)  = c0
            tsorc(jrow,n,mask) = c0
          enddo
        enddo
      enddo
# endif

# if defined bryan_lewis_vertical || defined bryan_lewis_horizontal

!-----------------------------------------------------------------------
!     initialize Bryan_Lewis tracer diffusion coefficients
!-----------------------------------------------------------------------

      call blmixi
# endif

# if defined trajectories

!-----------------------------------------------------------------------
!     initialize particle trajectories
!-----------------------------------------------------------------------

      if (initpt) call ptraji

# endif
# if defined time_averages

!-----------------------------------------------------------------------
!     initialize time mean "averaging" grid data
!-----------------------------------------------------------------------

      call avgi

# endif
# if defined xbts

!-----------------------------------------------------------------------
!     initialize XBT locations and averaging arrays
!-----------------------------------------------------------------------

      call xbti

# endif
# if defined uvic_tbt

!-----------------------------------------------------------------------
!     initialize tbt averaging arrays
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          do n=1,nt
            do k=1,kmt(i,j)
              do m=1,11
                tbt(i,j,k,n,m) = 0.
              enddo
            enddo
            tbtsf(i,j,n) = 0.
          enddo
        enddo
      enddo
      ntbtts = 0

# endif
# if defined ppvmix

!-----------------------------------------------------------------------
!     initialize pacanowski-philander vertical mixing scheme
!-----------------------------------------------------------------------

      call ppmixi (vmixset)
# endif
# if defined smagnlmix

!-----------------------------------------------------------------------
!     initialize Smagorinsky nonlinear horizontal mixing scheme
!-----------------------------------------------------------------------

      call smagnli (hmixset)
# endif
# if defined isopycmix || defined isoneutralmix

!-----------------------------------------------------------------------
!    initialize the isopycnal mixing
!-----------------------------------------------------------------------

      call isopi (error, am, ah)
#  if defined held_larichev

!-----------------------------------------------------------------------
!     initialize held_larichev horizontal mixing coefficients to
!     replace constant isopycnal mixing coefficient
!-----------------------------------------------------------------------

      call hlmixi
#  endif
# endif
# if defined uvic_npzd
#  if !defined uvic_embm
      print*, ' '
      print*, '==> Error: the npzd model needs short wave radiation'
      print*, '           from the embm. apply the uvic_embm option'
      print*, '           or modify the code and remove the stop.'
      print*, ' '
      stop '=>setocn'
#  endif
!     convert units of NPZD parameters to MOM units
      redctn = redctn*1.e-3
      redotn = redotn*1.e-3
      redotp = redotn/redptn
      redctp = redctn/redptn
      redntp = 1./redptn
      k1p   = k1*redptn
      kw = kw*1.e-2
      kc = kc*1.e-2
      ki = ki*1.e-2
      wd0 = wd0*1.e2
      alpha = alpha/daylen
      abio = abio/daylen
      nup = nup/daylen
      nupt0 = nupt0/daylen
      gbio = gbio/daylen
      epsbio = epsbio/daylen
      nuz = nuz/daylen
      gamma2 = gamma2/daylen
      nud0 = nud0/daylen
!     calculate sinking speed of detritus divided by grid width
      do k=1,km
!       linear increase wd0-200m with depth
        wd(k) = (wd0+4.0e-2*zt(k))/daylen/dzt(k)  ! [s-1]
        ztt(k) = -zt(k)+ dzt(k)/2.
        rkwz(k) = 1./(kw*dzt(k))
      enddo

#  if defined uvic_carbon || defined uvic_alk

!---------------------------------------------------------------------
!     calculate variables used in calcite remineralization
!---------------------------------------------------------------------

      rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
      rcab(1) = -1./dzt(1)
      do k=2,km
        rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
     &          + (exp(-zw(k-1)/dcaco3))/dzt(k)
        rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
      enddo
#  endif
# endif
# if defined uvic_fwa

!-----------------------------------------------------------------------
!     calculate areas for fresh water anomalies
!-----------------------------------------------------------------------

      areafwa = 0.
      areafwc = 0.
      if ( mrfwa .gt. 0 .and. mrfwa .le. nhreg) then
        do j=2,jmtm1
          if (j .ge. jsfwa .and. j .le. jefwa) then
            do i=2,imtm1
              if (kmt(i,j) .gt. 0) then
                if (mskhr(i,j) .eq. mrfwa) then
                  areafwa = areafwa + dxt(i)*dyt(j)*cst(j)
                else
                  areafwc = areafwc + dxt(i)*dyt(j)*cst(j)
                endif
              endif
            enddo
          else
            do i=2,imtm1
              if (kmt(i,j) .gt. 0) then
                areafwc = areafwc + dxt(i)*dyt(j)*cst(j)
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
     &            (i .ge. isfwa2 .and. i .le. iefwa2)) then
                  areafwa = areafwa + dxt(i)*dyt(j)*cst(j)
                else
                  areafwc = areafwc + dxt(i)*dyt(j)*cst(j)
                endif
              endif
            enddo
          else
            do i=2,imtm1
              if (kmt(i,j) .gt. 0) then
                areafwc = areafwc + dxt(i)*dyt(j)*cst(j)
              endif
            enddo
          endif
        enddo
      endif
# endif

!-----------------------------------------------------------------------
!     do all consistency checks last
!-----------------------------------------------------------------------

      call checks (error, vmixset, hmixset)

      return

9341  format (//,'  Global ocean statistics:'
     &,/,'  the total ocean surface area (t cells) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (t cells)       =',1pe15.8,'cm**3'
     &,/,'  the total ocean surface area (u cells) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (u cells)       =',1pe15.8,'cm**3')
9601  format(/,' "dtxcel(km)" tracer timestep multipliers:',/,10(f9.3))
      end

      subroutine depth_u (kmt, imt, jmt, zw, km, kmu, h, hr)

!=======================================================================
!     calculate depth arrays associated with "u" cells.

!     input:
!       kmt = number of oecan "t" cells from surface to bottom of ocean
!       imt = longitudinal dimension of arrays
!       jmt = latitudinal dimension of arrays
!       zw  = depth to bottom of "t" cells
!       km  = max number of depths

!     output:
!       kmu = number of ocean "u" cells from surface to bottom of ocean
!       h   = depth to ocean floor over "u" cells (cm)
!       hr  = reciprocal of "h"

!     based on code by: R. C. Pacanowski
!=======================================================================

      dimension kmt(imt,jmt), kmu(imt,jmt), h(imt,jmt), hr(imt,jmt)
      dimension zw(km)
# include "task_on.h"

!-----------------------------------------------------------------------
!     set some constants
!-----------------------------------------------------------------------

      c0 = 0.0
      c1 = 1.0

!-----------------------------------------------------------------------
!     compute number of vertical levels on the "u" grid
!-----------------------------------------------------------------------

      do jrow=1,jmt
        kmu(imt,jrow) = 0
      enddo

      do i=1,imt
        kmu(i,jmt) = 0
      enddo

      do jrow=1,jmt-1
        do i=1,imt-1
           kmu(i,jrow) = min (kmt(i,jrow), kmt(i+1,jrow), kmt(i,jrow+1)
     &,                       kmt(i+1,jrow+1))
        enddo
      enddo
# if defined obc_north
      do i=1,imt
        kmu(i,jmt) = kmu(i,jmt-1)
      enddo
# endif
# if defined obc_east
      do jrow=1,jmt
        kmu(imt,jrow) = kmu(imt-1,jrow)
      enddo
# endif
# if defined cyclic
      do jrow=1,jmt
        kmu(imt,jrow) = kmu(2,jrow)
      enddo
# endif
# if defined symmetry
      do i=1,imt
        kmu(i,jmt) = kmu(i,jmt-2)
      enddo
# endif

!---------------------------------------------------------------------
!     compute depths and reciprocal depths over "u" cells
!---------------------------------------------------------------------

      do jrow=1,jmt
        do i=1,imt
          hr(i,jrow) = c0
          h(i,jrow)  = c0
          if (kmu(i,jrow) .ne. 0) then
            hr(i,jrow) = c1/zw(kmu(i,jrow))
            h (i,jrow) = zw(kmu(i,jrow))
          endif
        enddo
      enddo

      return
      end

      subroutine rowi

!-----------------------------------------------------------------------
!     initialize prognostic quantities at "tau-1" and "tau"

!     inputs:

!     kmt  = number of vertical levels on "t" cells
!     yt   = latitudes of "t" points
!     zt   = depths of "t" points
!-----------------------------------------------------------------------

# include "param.h"
# include "coord.h"
# include "csbc.h"
# include "iounit.h"
# include "levind.h"
# include "mw.h"
# include "task_on.h"
# include "taskrows.h"
      character(120) :: fname, new_file_name
      integer iou
      integer ib(10), ic(10)
      logical exists
      real tmpik(imt,km), p1, c100, c1e3, dmsk(imt,jmt)
# if defined obctest
      real tempsouth(km),tempnorth(km),saltsouth(km),saltnorth(km)
      data tempsouth / 10.30,  7.32,  4.74,  2.50,  0.52, -0.65 /
      data tempnorth /  4.60, -0.54, -1.27, -1.10, -1.15, -1.20 /
      data saltsouth / 35.053, 35.113, 34.985, 34.914, 34.915, 34.912 /
      data saltnorth / 34.500, 34.782, 34.873, 34.894, 34.984, 34.889 /
# endif

!-----------------------------------------------------------------------
!     update pointers to tau-1, tau, & tau+1 data on disk.
!     for latitude rows they point to latdisk(1) or latdisk(2)
!     for 2D fields they point to records on kflds
!-----------------------------------------------------------------------

      itt   = 0
# if defined coarse_grained_parallelism
      taum1disk = mod(itt+2,3) + 1
      taudisk   = mod(itt  ,3) + 1
      taup1disk = mod(itt+1,3) + 1
# else
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk
# endif

      p1 = 0.1
      c100 = 100.
      c1e3 = 1000.

!-----------------------------------------------------------------------
!     update pointers to tau-1, tau, & tau+1 data in the MW
!-----------------------------------------------------------------------

      if (wide_open_mw) then

!       rotate time levels instead of moving data

        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
      endif

      ucksum = 0.0
      vcksum = 0.0
      tcksum = 0.0
      scksum = 0.0

!-----------------------------------------------------------------------
!     initialize every latitude jrow either in the MW (when wide opened)
!     or on disk (when jmw < jmt)
!-----------------------------------------------------------------------

      do jrow=1,jmt

        if (wide_open_mw) then
          j = jrow
        else
          j = jmw
        endif

!-----------------------------------------------------------------------
!       zero out all variables. velocities are internal modes only
!-----------------------------------------------------------------------

        do k=1,km
          do i=1,imt
            u(i,k,j,1,taup1) = c0
            u(i,k,j,2,taup1) = c0
          enddo
        enddo
        do n=1,nt
          do k=1,km
            do i=1,imt
              t(i,k,j,n,taup1) = c0
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       set tracers
!-----------------------------------------------------------------------

        do n=1,nt

          do i=1,imt
            do k=1,kmt(i,jrow)
              if (trim(mapt(n)) .eq. 'temp') then
                t(i,k,j,n,taup1) = theta0 (yt(jrow), zt(k))
              elseif (trim(mapt(n)) .eq. 'salt') then
                t(i,k,j,n,taup1) = 0.03472 - 0.035
              elseif (trim(mapt(n)) .eq. 'dic') then
                t(i,k,j,n,taup1) = 2.315
              elseif (trim(mapt(n)) .eq. 'alk') then
                t(i,k,j,n,taup1) = 2.429
              elseif (trim(mapt(n)) .eq. 'o2') then
                t(i,k,j,n,taup1) = 0.1692
              elseif (trim(mapt(n)) .eq. 'po4') then
                if (k .eq. 1) then
                  t(i,k,j,n,taup1) = 0.543
                else
                  t(i,k,j,n,taup1) = 2.165
                endif
              elseif (trim(mapt(n)) .eq. 'phyt') then
                t(i,k,j,n,taup1) = 0.14*exp((zt(1)-zt(k))/100.e2)
              elseif (trim(mapt(n)) .eq. 'zoop') then
                t(i,k,j,n,taup1) = 0.014*exp((zt(1)-zt(k))/100.e2)
              elseif (trim(mapt(n)) .eq. 'detr') then
                t(i,k,j,n,taup1) = 1.e-4
              elseif (trim(mapt(n)) .eq. 'no3') then
                if (k .eq. 1) then
                  t(i,k,j,n,taup1) = 5.30
                else
                  t(i,k,j,n,taup1) = 30.84
                endif
              elseif (trim(mapt(n)) .eq. 'diaz') then
                t(i,k,j,n,taup1) = 0.014*exp((zt(1)-zt(k))/100.e2)
              elseif (trim(mapt(n)) .eq. 'c14') then
                t(i,k,j,n,taup1) = -150 ! permil (converted later)
              elseif (trim(mapt(n)) .eq. 'cfc11') then
                t(i,k,j,n,taup1) = 0.0
              elseif (trim(mapt(n)) .eq. 'cfc12') then
                t(i,k,j,n,taup1) = 0.0
              else
                if (k .eq. 1) then
                  t(i,k,j,n,taup1) = c1
                else
                  t(i,k,j,n,taup1) = c0
                endif
              endif
            enddo
          enddo

# if !defined idealized_ic && !defined equatorial_thermocline
          ib(:) = 1
          ib(2) = jrow
          ib(3) = 1
          ic(:) = imt
          ic(2) = 1
          ic(3) = km
          if (trim(mapt(n)) .eq. 'temp') then
            C2K = 273.15
!           expecting units of degrees K
            fname = new_file_name ("temperature_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,itemp,taup1)
!             convert to model units (degrees C)
              call getvara ('temperature', iou, imt*km, ib, ic, tmpik
     &,         c1, -C2K)
              t(1:imt,1:km,j,itemp,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          elseif (trim(mapt(n)) .eq. 'salt') then
!           expecting units of psu
            fname = new_file_name ("salinity_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              p001 = 0.001
              p035 = 0.035
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,isalt,taup1)
!             convert to model units ((psu-35)/1000)
              call getvara ('salinity', iou, imt*km, ib, ic, tmpik
     &,         p001, -p035)
              t(1:imt,1:km,j,isalt,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          elseif (trim(mapt(n)) .eq. 'dic') then
!           expecting units of mol m-3 (equals model units, umol cm-3)
            fname = new_file_name ("tco2_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,idic,taup1)
              call getvara ('tco2', iou, imt*km, ib, ic, tmpik, c1, c0)
              t(1:imt,1:km,j,idic,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          elseif (trim(mapt(n)) .eq. 'alk') then
!           expecting units of mol m-3 (equals model units, umol cm-3)
            fname = new_file_name ("alk_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,ialk,taup1)
              call getvara ('alk', iou, imt*km, ib, ic, tmpik, c1, c0)
              t(1:imt,1:km,j,ialk,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          elseif (trim(mapt(n)) .eq. 'o2') then
!           expecting units of mol m-3 (equals model units, umol cm-3)
            fname = new_file_name ("oxygen_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,io2,taup1)
              call getvara ('oxygen', iou, imt*km, ib, ic, tmpik
     &,         c1, c0)
              t(1:imt,1:km,j,io2,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          elseif (trim(mapt(n)) .eq. 'po4') then
!           expecting units of mol m-3
            fname = new_file_name ("phosphate_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,ipo4,taup1)
!             convert to model units (nmol cm-3)
              call getvara ('phosphate', iou, imt*km, ib, ic, tmpik
     &,         c1e3, c0)
              t(1:imt,1:km,j,ipo4,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          elseif (trim(mapt(n)) .eq. 'no3') then
!           expecting units of mol m-3
            fname = new_file_name ("nitrate_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,ino3,taup1)
!             convert to model units (nmol cm-3)
              call getvara ('nitrate', iou, imt*km, ib, ic, tmpik
     &,         c1e3, c0)
              t(1:imt,1:km,j,ino3,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          elseif (trim(mapt(n)) .eq. 'c14') then
!           expecting units of permil
            fname = new_file_name ("bkgc14_ann.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              call openfile (trim(fname), iou)
              tmpik(1:imt,1:km) = t(1:imt,1:km,j,ic14,taup1)
!             convert later
              call getvara ('bkgc14', iou, imt*km, ib, ic, tmpik
     &,         c1, c0)
              t(1:imt,1:km,j,ic14,taup1) = tmpik(1:imt,1:km)
              call closefile (iou)
            else
              if (jrow .eq. 1)
     &          print*, "Warning => Can not find ",trim(fname)
            endif
          endif

          do k=1,km
            do i=1,imt
              if (kmt(i,jrow) .lt. k)  t(i,k,j,n,taup1) = c0
            enddo
          enddo
# endif

        enddo
# if defined uvic_carbon && defined uvic_carbon_14

!       convert c14 from permil to model units (umol cm-3)
        do k=1,km
          do i=1,imt
            t(i,k,j,ic14,taup1) = (t(i,k,j,ic14,taup1)*0.001 + 1)*rstd
     &                            *t(i,k,j,idic,taup1)
          enddo
        enddo
# endif

!-----------------------------------------------------------------------
!       initialize surface boundary
!-----------------------------------------------------------------------
        if (isst .ne. 0) sbc(1:imt,jrow,isst) = t(1:imt,1,j,itemp,taup1)
        if (isss .ne. 0) sbc(1:imt,jrow,isss) = t(1:imt,1,j,isalt,taup1)
        if (issdic .ne. 0)
     &    sbc(1:imt,jrow,issdic) = t(1:imt,1,j,idic,taup1)
        if (isso2 .ne. 0) sbc(1:imt,jrow,isso2) = t(1:imt,1,j,io2,taup1)
        if (issalk .ne. 0)
     &    sbc(1:imt,jrow,issalk) = t(1:imt,1,j,ialk,taup1)
        if (isspo4 .ne. 0)
     &    sbc(1:imt,jrow,isspo4) = t(1:imt,1,j,ipo4,taup1)
        if (issphyt .ne. 0)
     &    sbc(1:imt,jrow,issphyt) = t(1:imt,1,j,iphyt,taup1)
        if (isszoop .ne. 0)
     &    sbc(1:imt,jrow,isszoop) = t(1:imt,1,j,izoop,taup1)
        if (issdetr .ne. 0)
     &    sbc(1:imt,jrow,issdetr) = t(1:imt,1,j,idetr,taup1)
        if (issno3 .ne. 0)
     &    sbc(1:imt,jrow,issno3) = t(1:imt,1,j,ino3,taup1)
        if (issdiaz .ne. 0)
     &    sbc(1:imt,jrow,issdiaz) = t(1:imt,1,j,idiaz,taup1)
        if (issc14 .ne. 0)
     &    sbc(1:imt,jrow,issc14) = t(1:imt,1,j,ic14,taup1)
        if (isscfc11 .ne. 0)
     &    sbc(1:imt,jrow,isscfc11) = t(1:imt,1,j,icfc11,taup1)
        if (isscfc12 .ne. 0)
     &    sbc(1:imt,jrow,isscfc12) = t(1:imt,1,j,icfc12,taup1)

# if defined linearized_advection || defined equatorial_thermocline

!       initialize density profile with tbarz

        do n=1,nt
          do i=1,imt
            do k=1,km
              t(i,k,j,n,taup1) = tbarz(k,n)
            enddo
          enddo
        enddo
        if (jrow .eq. jmt) then
          write (stdout,'(a,a)')
     & '=> Note: initialized T & S with idealized profile T(z), S=35ppt'
        endif
# endif
# if defined obctest

!       construct idealized initial density profile a la Stevens (1990)

        do i=1,imt
          do k=1,km
            if ( jrow .le. (jmt/2) ) then
               t(i,k,j,1,taup1) = tempsouth(k)
               t(i,k,j,2,taup1) = (saltsouth(k)-35.0)*0.001
            else
               t(i,k,j,1,taup1) = tempnorth(k)
               t(i,k,j,2,taup1) = (saltnorth(k)-35.0)*0.001
            endif
          enddo
        enddo
        if (jrow .eq. jmt) then
          write (stdout,'(a,a)')
     & '=> Note: initialized T & S to Stevens (1990) first test case'
     &,                    '                        '
        endif
# endif
# if defined obctest2

!       construct idealized two basins a la Stevens (1990), case two

        do i=1,imt
          do k=1,km
            t(i,k,j,1,taup1) = 6.
            t(i,k,j,2,taup1) = (34.9-35.0)*0.001
          enddo
        enddo

        if (jrow .eq. jmt) then
          write (stdout,'(a,a)')
     & '=> Note: initialized T & S to Stevens (1990) second test case'
     &,                    '                        '
        endif
# endif

!-----------------------------------------------------------------------
!       zero out tracers in land points
!-----------------------------------------------------------------------

        do i=1,imt
          kz = kmt(i,jrow)
          do k=1,km
            if (k .gt. kz) then
              do n=1,nt
                t(i,k,j,n,taup1) = c0
              enddo
            endif
          enddo
        enddo

!-----------------------------------------------------------------------
!       checksum the initial conditions
!-----------------------------------------------------------------------

        ucksum = ucksum + checksum (u(1,1,j,1,taup1), imt, km)
        vcksum = vcksum + checksum (u(1,1,j,2,taup1), imt, km)
        tcksum = tcksum + checksum (t(1,1,j,1,taup1), imt, km)
        scksum = scksum + checksum (t(1,1,j,2,taup1), imt, km)

!-----------------------------------------------------------------------
!       initialize every latitude jrow either on disk (when jmw < jmt)
!       or in the MW (when the last jrow is complete and jmw = jmt)
!-----------------------------------------------------------------------

        if (wide_open_mw) then
          if (jrow .eq. jmt) then
            call copy_all_rows (taup1, tau)
            call copy_all_rows (tau, taum1)
          endif
        else
# if defined coarse_grained_parallelism
          do n=1,num_processors
            if (jrow .ge. jstask(n) .and. jrow .le. jetask(n)) then
# endif
              call putrow (latdisk(taudisk),  nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
              call putrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
# if defined coarse_grained_parallelism
              if (jrow .eq. 1 .or. jrow .eq. jmt) then
                call putrow (latdisk(taum1disk), nslab, jrow
     &,                  u(1,1,j,1,taup1), t(1,1,j,1,taup1))
              endif
# endif
# if defined coarse_grained_parallelism
            endif
          enddo
# endif
        endif

      enddo

!-----------------------------------------------------------------------
!     find inital average surface references
!-----------------------------------------------------------------------
      print*, " "
      print*, "inital average surface references: "
      dmsk(:,:) = 1.
      if (isss .ne. 0) then
        call areaavg (sbc(1,1,isss), dmsk, socn)
        socn = socn + 0.035
        print*, "salinity, socn = ", socn
      endif
      if (issdic .ne. 0) then
        call areaavg (sbc(1,1,issdic), dmsk, dicocn)
        print*, "dissolved inorganic carbon, dicocn = ", dicocn
      endif
      if (isso2 .ne. 0) then
        call areaavg (sbc(1,1,isso2), dmsk, o2ocn)
        print*, "oxygen, o2oc = ", o2ocn
      endif
      if (issalk .ne. 0) then
        call areaavg (sbc(1,1,issalk), dmsk, alkocn)
        print*, "alkalinity, alkocn = ", alkocn
      endif
      if (isspo4 .ne. 0) then
        call areaavg (sbc(1,1,isspo4), dmsk, po4ocn)
        print*, "phosphate, po4ocn = ", po4ocn
      endif
      if (issphyt .ne. 0) then
        call areaavg (sbc(1,1,issphyt), dmsk, phytocn)
        print*, "phytoplankton, phytocn = ", phytocn
      endif
      if (isszoop .ne. 0) then
        call areaavg (sbc(1,1,isszoop), dmsk, zoopocn)
        print*, "zooplankton, zoopocn = ",zoopocn
      endif
      if (issdetr .ne. 0) then
        call areaavg (sbc(1,1,issdetr), dmsk, detrocn)
        print*, "detritus, detocn = ", detocn
      endif
      if (issno3 .ne. 0) then
        call areaavg (sbc(1,1,issno3), dmsk, no3ocn)
        print*, "nitrate, no3ocn = ", no3ocn
      endif
      if (issdiaz .ne. 0) then
        call areaavg (sbc(1,1,issdiaz), dmsk, diazocn)
        print*, "diazotrophs, diazocn = ", diazocn
      endif
      if (issc14 .ne. 0) then
        call areaavg (sbc(1,1,issc14), dmsk, c14ocn)
        print*, "carbon 14, c14ocn = ", c14ocn
      endif
      if (isscfc11 .ne. 0) then
        call areaavg (sbc(1,1,isscfc11), dmsk, cfc11ocn)
        print*, "chlorofluorocarbon 11, cfc11ocn = ", cfc11ocn
      endif
      if (isscfc12 .ne. 0) then
        call areaavg (sbc(1,1,isscfc12), dmsk, cfc12ocn)
        print*, "chlorofluorocarbon 12, cfc12ocn = ", cfc12ocn
      endif
      print*, " "

      write (stdout,*) ' I.C. checksum for t =',tcksum
      write (stdout,*) ' I.C. checksum for s =',scksum
      write (stdout,*) ' I.C. checksum for u =',ucksum
      write (stdout,*) ' I.C. checksum for v =',vcksum
      return
      end

      subroutine sethr (nr, xstart, xend, ystart, yend)

!=======================================================================
!     discretizes the horizontal region to nearest model grid points

!     nr     = the horizontal region number
!     xstart = starting longitude at edge of "t" box region
!     xend   = ending longitude at edge of  "t" box region
!     ystart = starting latitude at edge of  "t" box region
!     yend   = ending latitude at edge of "t" box region
!=======================================================================

# include "param.h"
# include "coord.h"
# include "cregin.h"
# include "levind.h"

!     find the nearest "t" box indices within the region

      jsr = min(indp (ystart, yu, jmt)+1, jmt)
      jer = indp (yend, yu, jmt)

      isr = min(indp (xstart, xu, imt)+1, imt)
      ier = indp (xend, xu, imt)

!     define "edges" of "t" box region

      if (isr .eq. 1) then
        xsrl = xu(1) - dxudeg(1)
      else
        xsrl = xu(isr-1)
      endif
      xerl = xu(ier)
      if (jsr .eq. 1) then
        ysrl = yu(1) - dyudeg(1)
      else
        ysrl = yu(jsr-1)
      endif
      yerl = yu(jer)

      write (hregnm(nr),9000) xsrl, xerl, ysrl, yerl
      write (stdout,*) ' Defining horizontal region # ',nr
     &, ' as "t" cells within ', hregnm(nr)
      if (isr .gt. ier) then
        write (stdout,*) ' Error: isr=',isr,' >  ier=',ier,' in sethr'
        stop '=>sethr'
      endif
      if (jsr .gt. jer) then
        write (stdout,*) ' Error: jsr=',jsr,' >  jer=',jer, 'in sethr'
        stop '=>sethr'
      endif
      do j=jsr,jer
        do i=isr,ier
          if (kmt(i,j) .gt. 0)  then
            mskhr(i,j) = nr
          endif
        enddo
      enddo
9000  format ('lon: ',f5.1,' => ',f5.1,'  lat: ',f5.1,' => ',f5.1)
      return
      end

      subroutine setvr (nr, zstart, zend)

!=======================================================================
!     discretizes the vertical region to nearest model grid points

!     nr     = the vertical region number
!     zstart = starting depth at edge of "t" box region in cm.
!     zend   = ending depth at edge of "t" box region in cm.
!=======================================================================

# include "param.h"
# include "coord.h"
# include "cregin.h"

!     find the nearest "t" box indices within the region

      if (zstart .lt. p5*zw(1)) then
        ksr = 1
        ztopb = 0.0
      else
        ksr = min(indp (zstart, zw, km)+1, km)
        ztopb = zw(ksr-1)*0.01
      endif
      ker = indp (zend, zw, km)
      llvreg(nr,1) = ksr
      llvreg(nr,2) = ker

      write (vregnm(nr),9000) ztopb, zw(ker)*0.01
      write (stdout,*) ' Defining vertical region # ',nr
     &,  ' as "t" cells within ',  vregnm(nr)
      if (ksr .gt. ker) then
        write (stdout,*) ' Error: ksr=',ksr,' >  ker=',ker, 'in setvr'
        stop '=>setvr'
      endif
9000  format (' dpt:',f6.1, '=>',f6.1, 'm')
      return
      end

      function theta0 (ydeg, depth)

!=======================================================================
!     this subroutine returns estimates of global mean potential
!     temperature for model initialization as a function of depth.
!     it is used to produce a reference thermal stratification for the
!     upper 2000m of the MOM`s test case.  below 2000m, the
!     potential temperature returned is 2.0 degrees C.  surface
!     values are set slightly above 18.4 degrees C at the reference
!     latitude "reflat".
!     the estimates are produced from a 7th order polynomial fit to
!     the annual mean world ocean potential temperature observations
!     of Levitus (1982).

!     input [units]:
!       a latitude (ydeg): [degrees]
!       a zt value (depth): [centimetres]
!     output [units]:
!       potential temperature estimate (est): [degrees centigrade]

!     variables:
!       coeft     = coefficients for the polynomial fit of potential
!                   temperature vs. depth
!       reflat    = reference latitude at which observed surface
!                   temperatures approximately equal coeft(1)
!       factor    = the ratio of the cosine of the latitude requested
!                   ("ydeg") to the reference latitude ("reflat")
!                   used to scale the upper 2000 meters of the vertical
!                   temperature profile
!       tmin,tmax = the minimum and maximum potential temperatures
!                   allowed at the time of model initialization

!     reference:
!       Levitus, S., Climatological atlas of the world ocean, NOAA
!     Prof. Paper 13, US Gov`t printing Office, Washington, DC, 1982.

      parameter (ndeg=7)
      dimension coeft(ndeg+1)
      save coeft, tmin, tmax, reflat
      data coeft / 0.184231944E+02,-0.430306621E-01, 0.607121504E-04
     &           ,-0.523806281E-07, 0.272989082E-10,-0.833224666E-14
     &           , 0.136974583E-17,-0.935923382E-22/
      data tmin, tmax, reflat /2.0, 25.0, 34.0/

!=======================================================================

      c0 = 0.0
      pi = atan(1.0) * 4.0
      refcos = abs(cos(pi*reflat/180.))

      coslat = abs(cos(pi*ydeg/180.))
      factor = coslat/refcos
      z = depth * 0.01

      if (z .gt. 2000.) then
        est = 2.0
      else
        est = c0
        bb = 1.0
        do nn=1,ndeg+1
!          if (nn.gt.1) bb = z**(nn-1)
          est = est + coeft(nn)*bb
          bb = bb*z
        enddo
        est = est * factor
      endif

      if (est .gt. tmax) est = tmax
      if (est .lt. tmin) est = tmin

      theta0 = est

      return
      end
#endif
