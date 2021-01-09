      subroutine setembm

!=======================================================================
!     initialize the energy-moisture balance model

!     based on code by: A. Fanning and M. Eby
!=======================================================================
#if defined uvic_embm

      implicit none

# include "param.h"
# include "calendar.h"
# include "solve.h"
# include "switch.h"
# include "coord.h"
# include "grdvar.h"
# include "cembm.h"
# include "atm.h"
# include "insolation.h"
# if defined uvic_ice
#  include "ice.h"
#  if defined uvic_ice_cpts
#   include "thermo.h"
#  endif
#  if defined uvic_ice_evp
#   include "evp.h"
#  endif
# endif
# include "riv.h"
# include "tmngr.h"
# include "levind.h"
# include "csbc.h"
# include "scalar.h"
# include "veg.h"
      real rveg(imt,jmt)
# if defined uvic_embm_annual
      real cosz(imt,jmt)
# endif
# if defined uvic_embm_awind
      real tmpij(2:imt-1,2:jmt-1)
# endif

      character(120) :: fname, new_file_name
      character(3) :: a3

      integer i, ie, ii, iou, is, j, je, jj, js, jz, k, m, n, nsolve
      integer nu, nsum, iice, ib(10), ic(10)

      logical exists, inqvardef

      real dlam, dphi, dtatms, dte, dyz, eccice, grarea, saltmax
      real si, ssh, t1, tair, yz_max, yz_min, wz, calday
      real zrel, C2K, c100
      real dmsk(imt,jmt)

# if defined uvic_ice_cpts
      logical rowflg(ncat) ! flag for computing ridg. matrix row
      real Hi(ncat)        ! ice thickness (m)
      real Hmean(ncat)     ! a dummy variable at setup (m)

# endif
      namelist /tsteps/ dtts, dtuv, dtsf, dtatm, dtatms, namix, segtim
      namelist /embm/   rlapse, rfactor, scatter, rhmax
      namelist /carbon/ co2ccn, co2yri, co2yrf, co2rate, co2for, co2_yr
     &,                 c14_yr, dc14ccn, c14prod
      namelist /paleo/  orbit_yr, pyear, eccen, obliq, mvelp
      namelist /ice/    niats, nivts, dampice, tsno, hsno_max, ice_yr
     &,                  ice_calb, sno_calb
      namelist /veg/    veg_alb, veg_rl, veg_rs, veg_smd, crops_yr
     &,                 icrops, idesert, iice
      namelist /solar/  solarconst, solar_yr

      is = 1
      ie = imt
      js = 1
      je = jmt

      c100 = 100.
      C2K = 273.15

!-----------------------------------------------------------------------
!     define default parameters
!-----------------------------------------------------------------------

!     set defaults namelist embm
      rlapse     = 5.e-5
      rfactor    = 0.2
      scatter    = 0.23
      rhmax      = 0.85

!     set defaults namelist carbon
      co2ccn     = 280
      co2yri     = 1.e20
      co2yrf     = 1.e20
      co2rate    = 0.
      co2for     = 5.35e03
      co2_yr     = 2.e20
      c14_yr     = 2.e20
      dc14ccn    = 0.
      c14prod    = 0.

!     set defaults namelist paleo
      pyear      = 2.e20
      orbit_yr   = 2.e20
!     default user orbit is for 1950
      eccen      = 1.672393E-02
      obliq      = 23.44627
      mvelp      = 102.0391

!     set defaults namelist ice
      niats      = 1
      nivts      = 1
      dampice    = 5.
      tsno       = 0.
      hsno_max   = 1000.
# if defined uvic_embm_snow_transient
!     no limit
      hsno_max   = 1.e30
# endif
      ice_yr     = 2.e20
      ice_calb   = 0.35
      sno_calb   = 0.25

!     set defaults namelist veg
      veg_alb(1:7) = (/0.17, 0.17, 0.22, 0.22, 0.22, 0.30, 0.60/)
      veg_rl(1:7) = (/2.9, 0.91, 0.11, 0.053, 0.039, 0.039, 0.020/)
      veg_rs(1:7) = (/0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
      veg_smd(1:7) = (/2.0, 3.0, 0.1, 0.1, 0.5, 0.01, 0.01/)
      crops_yr = 2.e20
      icrops = 3
      idesert = 6
      iice = 7

!     set defaults namelist solar
      solarconst = 1.368e6
      solar_yr = 2.e20

!     set other stuff
      cdatm     = 1.e-3
      cpatm     = 1.004e7
      sht       = 8.4e5
      shq       = 1.8e5
      shc       = 7.9e5
      rhoatm    = 1.250e-3
      esatm     = 4.6e-05
      pcfactor  = 0.
      cssh      = 3.8011e-3
      cfc11ccnn = 0.
      cfc11ccns = 0.
      cfc12ccnn = 0.
      cfc12ccns = 0.
      dc14ccnn  = 0.
      dc14ccne  = 0.
      dc14ccns  = 0.

      rhoocn    = 1.035
      esocn     = 5.4e-5
      vlocn     = 2.501e10

      cdice     = 5.5e-3
      rhoice    = 0.913
      rhosno    = 0.330
      esice     = 5.347e-5
      slice     = 2.835e10
      flice     = 3.34e9
      condice   = 2.1656e5

      soilmax   = 15.
      eslnd     = 5.347e-5
      cs_alb    = 0.08

      nivc      = 1
      dtatms    = 1800.
      ns        = 30

      dalt_v    = 3.3e-3
      dalt_o    = 1.4e-3
      dalt_i    = 1.4e-3

!-----------------------------------------------------------------------
!     read namelist files
!-----------------------------------------------------------------------
      call getunit (iou, 'control.in', 'f s r')
      read  (iou, tsteps, end=100)
100   continue
      write (stdout, tsteps)
      call relunit (iou)

      call getunit (iou, 'control.in', 'f s r')
      read  (iou, embm, end=101)
101   continue
      write (stdout,embm)
      call relunit (iou)

      call getunit (iou, 'control.in', 'f s r')
      read  (iou, carbon, end=102)
102   continue
      write (stdout,carbon)
      call relunit (iou)
      co2ccni = co2ccn

      call getunit (iou, 'control.in', 'f s r')
      read  (iou, paleo, end=103)
103   continue
      write (stdout,paleo)
      call relunit (iou)

      call getunit (iou, 'control.in', 'f s r')
      read  (iou, ice, end=104)
104   continue
      write (stdout,ice)
      call relunit (iou)

      call getunit (iou, 'control.in', 'f s r')
      read  (iou, veg, end=105)
105   continue
      write (stdout,veg)
      call relunit (iou)

      call getunit (iou, 'control.in', 'f s r')
      read  (iou, solar, end=106)
106   continue
      write (stdout,solar)
      call relunit (iou)

!     ensure pass is between zero and one.
      pass =  min(max((1. - scatter), 0.), 1.)
# if !defined uvic_mom

!     if not coupled set ocean and coupling time to atmosphere
      dtts = dtatm
      dtocn = dtatm
      segtim = dtts*secday
# endif
# if defined uvic_carbon_14_coupled

!     calculate c14ccn from dc14ccn and co2ccn
      c14ccn = (1 + dc14ccn*0.001)*rstd*co2ccn
# endif

!     sort out forcing years
      if (pyear  .gt. 1.e20) pyear = 1850.
      if (co2_yr .gt. 1.e20) co2_yr = pyear
      if (c14_yr .gt. 1.e20) c14_yr = pyear
      if (orbit_yr .gt. 1.e20) orbit_yr = pyear
      if (crops_yr .gt. 1.e20) crops_yr = pyear
      if (ice_yr .gt. 1.e20) ice_yr = pyear
      if (solar_yr .gt. 1.e20) solar_yr = pyear

# if defined uvic_embm_explicit
      if (dtatms .ne. 0.) ns = nint(dtatm/dtatms)

# endif
# if defined time_averages
      if (mod(timavgint, segtim) .gt. 1.e-6 .and. timavgint .gt. 0.)
     &  then
        t1 = nint(timavgint/segtim)*segtim
        if (t1 .lt. segtim) t1 = segtim
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: "timavgint" does not contain an integral number'
     &,   '              of coupling time steps "segtim".              '
        write (stdout,*) '              (changed "timavgint" from '
     &, timavgint,' days to ', t1,' days to insure this condition)'
        timavgint = t1
      endif
      if (timavgint .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging interval "timavgint" = 0. implies no '
     &,   '             averaging when "time_averages" is enabled      '
      endif
      if (timavgint .gt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the interval "timavgint" exceeds the averaging '
     &,   '             period "timavgper" for option "time_averages"  '
      endif
      if (timavgint .lt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging period "timavgper" exceeds interval  '
     &,   '             "timavgint". Setting timavgper = timavgint     '
        timavgper = timavgint
      endif
      if (timavgper .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the averaging period "timavgper" is zero. The  '
     &,   '             average will be over only one time step!       '
      endif
      write (stdout,'(/,1x,a,f10.2,a,/,1x,a,f10.2,a)')
     &  '==> Time averages will be written every ', timavgint, ' days, '
     &, '    with an averaging period of         ', timavgper, ' days. '

# endif
# if !defined uvic_carbon_coupled
!-----------------------------------------------------------------------
!     calculate the relative CO2 forcing term
!-----------------------------------------------------------------------

      call co2forc

#  if !defined uvic_embm_co2_lin && !defined uvic_embm_co2_exp
      write (stdout,*)
      write (stdout,*) 'CO2 concentration (versus today) =',co2ccn/350.
      write (stdout,*) 'Yields radiative forcing (W/m2) = ',anthro*1.e-3

#  endif
# endif
!-----------------------------------------------------------------------
!     calculate the expansion coefficients for Berger's solution for
!     the year of the initial conditions
!-----------------------------------------------------------------------

# if defined uvic_embm_orbit_transient
      orbit_yr = year0 + relyr
# endif
      call orbit (orbit_yr, eccen, obliq, mvelp, lambm0)
      write (stdout,*)
      write (stdout,*) 'Initial Orbital Parameters:'
# if !defined uvic_embm_orbit_user
      write (stdout,*) '  Orbital Year:', orbit_yr
# endif
      write (stdout,*) '  Eccentricity:', eccen
      write (stdout,*) '  Obliquity:   ', obliq
      write (stdout,*) '  Longitude of Perihelion:', mvelp+180.

!-----------------------------------------------------------------------
!     calculate annual average insolation and Coriolis factor
!-----------------------------------------------------------------------

      radian = 360./(2.*pi)
      do j=1,jmt
#  if defined uvic_embm_explicit
        filter(j) = 1. - sin(yt(j)/radian)**ns
#  endif
        do i=1,imt
!         calculate coriolis parameter
          fcor(i,j) = 2.*omega*sin(ulat(i,j)/radian)
        enddo
      enddo

# if defined uvic_embm_annual
      solins(:,:) = 0.
      i = imt*jmt
      do n=1,365
!       subroutine decl is expecting a 365.25 day year
        calday = n*365.25/365.
        call decl (calday, eccen, obliq, mvelp, lambm0, sindec, eccf)
        call zenith (i, c0, daylen, daylen, tlat, tlon, sindec, cosz)
        solins(:,:) = solins(:,:) + solarconst*eccf*cosz(:,:)
      enddo
      solins(:,:) = solins(:,:)/365.
#  endif

!-----------------------------------------------------------------------
!     read diffusion
!-----------------------------------------------------------------------

      dn(:,:,:) = 1.e10
      de(:,:,:) = 1.e10
      fname = new_file_name ("diff.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        call openfile (fname, iou)
        ib(:) = 1
        ic(:) = imt
        ic(2) = jmt
        do n=1,nat
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
          call getvara ('dn_'//trim(a3), iou, imt*jmt, ib, ic, dn(1,1,n)
     &,     c1, c0)
          call getvara ('de_'//trim(a3), iou, imt*jmt, ib, ic, de(1,1,n)
     &,     c1, c0)
        enddo
        call closefile (iou)
      endif

!-----------------------------------------------------------------------
!     set solver parameters
!-----------------------------------------------------------------------

      nsolve = 0
# if !defined uvic_embm_explicit
      itin(1:nat)  = 500        ! max solver iterations
#  if defined uvic_global_sums
      epsin(:) = 5.e-12
#  else
      epsin(1) = 1.e-3
      epsin(2:nat) = 1.e-5
#  endif
# endif
# if defined uvic_embm_essl
      nsolve = nsolve + 1
      iparm(2)= 2               ! solver method
      iparm(3)= 7               ! if iparm(2)=3 then iparm(3)=k,5<k<10
      iparm(4)= 4               ! preconditioning
      iparm(5)= 2               ! stopping criterion
      call einfo(0)
      call errset(2110,256,-1,0,1,2110)
# endif
# if defined uvic_embm_sparskit
      nsolve = nsolve + 1
      ipar(1) = 0               ! always 0 to start an iterative solver
      ipar(2) = 2               ! right preconditioning
      ipar(3) = 1               ! use convergence test scheme 1
      ipar(4) = nwork           ! the 'w' has nwork elements
      ipar(5) = 10              ! use *GMRES(10) (e.g. FGMRES(10))
      ipar(6) = 100             ! use at most 100 matvec's
      fpar(1) = 1.0E-6          ! relative tolerance 1.0E-6
      fpar(2) = 1.0E-10         ! absolute tolerance 1.0E-10
      fpar(11) = 0.0            ! clearing the FLOPS counter
# endif
# if defined uvic_embm_slap
      nsolve = nsolve + 1
# endif
# if defined uvic_embm_mgrid
      nsolve = nsolve + 1
      levelin = 20              ! max coarse grid level
# endif
# if defined uvic_embm_adi
      nsolve = nsolve + 1
      itin(1:nat) = 1
      itout(1:nat) = 1
# endif
# if defined uvic_embm_explicit
      nsolve = nsolve + 1
# endif
      if (nsolve .ne. 1) then
        write(*,*) '==> Error: more or less than one solver defined.'
        write(*,*) '           Use only one of uvic_embm_adi,'
     &,   ' uvic_embm_mgrid, uvic_embm_slap, uvic_embm_essl or'
     &,   ' or uvic_embm_sparskit or uvic_embm_explicit'
        stop '=>setembm'
      endif

# if defined uvic_ice
!-----------------------------------------------------------------------
!     check latent heats will sum to zero
!-----------------------------------------------------------------------

      if (slice .ne. vlocn + flice) write (stdout,'(/,a)')
     &   '==> Warning: changing latent heat of fusion to conserve heat'
        flice = slice - vlocn

# else

          write (stdout,*) '==> Warning: ice model is not defined.'
     &,     ' heat flux may be limited to prevent freezing sst.'

# endif
# if defined uvic_embm_solve2x || defined uvic_embm_solve2y
!-----------------------------------------------------------------------
!     calculate grid ratio for the coarse grid atmospheric solver
!-----------------------------------------------------------------------

#  if defined uvic_embm_solve2y
      do jj=1,jjmtm2
        j = jj*2
#  else
      do j=2,jmtm1
#  endif
#  if defined uvic_embm_solve2x
        do ii=1,iimtm2
          i = ii*2
#  else
        do i=2,imtm1
#  endif
          grarea = dxt(i)*dyt(j)*cst(j)
#  if defined uvic_embm_solve2x
          grarea = grarea + dxt(i+1)*dyt(j)*cst(j)
#  endif
#  if defined uvic_embm_solve2y
          grarea = grarea + dxt(i)*dyt(j+1)*cst(j+1)
#  endif
#  if defined uvic_embm_solve2x && defined uvic_embm_solve2y
          grarea = grarea + dxt(i+1)*dyt(j+1)*cst(j+1)
#  endif
          gr(i,j) = dxt(i)*dyt(j)*cst(j)/grarea
#  if defined uvic_embm_solve2x
          gr(i+1,j) = dxt(i+1)*dyt(j)*cst(j)/grarea
#  endif
#  if defined uvic_embm_solve2y
          gr(i,j+1) = dxt(i)*dyt(j+1)*cst(j+1)/grarea
#  endif
#  if defined uvic_embm_solve2x && defined uvic_embm_solve2y
          gr(i+1,j+1) = dxt(i+1)*dyt(j+1)*cst(j+1)/grarea
#  endif
        enddo
      enddo

# endif
!-----------------------------------------------------------------------
!     calculate grid terms for the atmospheric solver
!-----------------------------------------------------------------------

# if defined uvic_embm_solve2y
      do jj=2,jjmtm2
        j = jj*2
        wtj(j) = 0.5*dyt(j)/(dyt(j-1)+dyt(j)+dyt(j+1)+dyt(j+2))
        wtj(j-1) = 0.5*dyt(j+1)/(dyt(j-1)+dyt(j)+dyt(j+1)+dyt(j+2))
        ygrd(j) = dyt(j)*cst(j)/(cst(j)*dyt(j)+cst(j+1)*dyt(j+1))
        ygrd(j+1) = dyt(j+1)*cst(j+1)/(cst(j+1)*dyt(j+1)+cst(j)*dyt(j))
      enddo

      do j=2,jmtm2
        dsgrd(j) = csu(j-1)/((dyt(j)+dyt(j-1))*csu(j)*
     &             (dyu(j)+dyu(j+1)))
        dngrd(j) = csu(j+1)/((dyt(j+2)+dyt(j+1))*csu(j)*
     &             (dyu(j)+dyu(j+1)))
        asgrd(j) = csu(j-1)/(2.*csu(j)*(dyu(j)+dyu(j+1)))
        angrd(j) = csu(j+1)/(2.*csu(j)*(dyu(j)+dyu(j+1)))
# else
      do j=2,jmtm1
        dsgrd(j) = csu(j-1)/(dyu(j-1)*cst(j)*dyt(j))
        dngrd(j) = csu(j)/(dyu(j)*cst(j)*dyt(j))
        asgrd(j) = csu(j-1)/(2.*cst(j)*dyt(j))
        angrd(j) = csu(j)/(2.*cst(j)*dyt(j))
# endif
      enddo
# if defined uvic_embm_solve2x

      do ii=1,iimtm2
        i = ii*2
        wti(i) = 0.5*dxt(i)/(dxt(i-1)+dxt(i)+dxt(i+1)+dxt(i+2))
        wti(i+1) = 0.5*dxt(i+1)/(dxt(i-1)+dxt(i)+dxt(i+1)+dxt(i+2))
        xgrd(i) = dxt(i)/(dxt(i) + dxt(i+1))
        xgrd(i+1) = dxt(i+1)/(dxt(i) + dxt(i+1))
      enddo

      do i=2,imtm2
        dwgrd(i) = 1./((dxt(i) + dxt(i-1))*(dxu(i) + dxu(i+1)))
        degrd(i) = 1./((dxt(i+2) + dxt(i+1))*(dxu(i) + dxu(i+1)))
        azgrd(i) = 1./(2.*(dxu(i) + dxu(i+1)))
# else
      do i=2,imtm1
        dwgrd(i) = 1./(dxu(i-1)*dxt(i))
        degrd(i) = 1./(dxu(i)*dxt(i))
        azgrd(i) = 1./(2.*dxt(i))
# endif
      enddo

# if defined uvic_ice_cpts
#  if !defined uvic_ice_cpts5 && !defined uvic_ice_cpts10 && defined roth_press
      write(*,*) 'you are strongly discouraged from using roth_press'
      stop ' with fewer than 5 ice categories'
#  endif
!-----------------------------------------------------------------------
!     setup the vectors identifying first and last layer in each bin
!-----------------------------------------------------------------------

#  if defined uvic_ice_cpts1
      nilay(1) = 4
#   if defined uvic_ice_cpts3 || defined uvic_ice_cpts5 || defined uvic_ice_cpts10
!     stop if uvic_ice_cpts1 and uvic_ice_cpts3, 5 or 10 is defined
      stop 'single and multiple ice categories options defined'
#   endif
#  elif defined uvic_ice_cpts3
      nilay(1) = 2
      nilay(2) = 4
      nilay(3) = 8

      hstar(1) = 50.
      hstar(2) = 250.
#  elif defined uvic_ice_cpts5
      nilay(1) = 2
      nilay(2) = 4
      nilay(3) = 4
      nilay(4) = 8
      nilay(5) = 8

      hstar(1) = 40.
      hstar(2) = 90.
      hstar(3) = 200.
      hstar(4) = 350.
#  elif defined uvic_ice_cpts10
      nilay(1) = 2
      nilay(2) = 2
      do n=3,5
        nilay(n) = 4
      enddo
      do n=6,ncat
        nilay(n) = 8
      enddo

      hstar(1) = 25.
      hstar(2) = 50.
      hstar(3) = 75.
      hstar(4) = 100.
      hstar(5) = 140.
      hstar(6) = 190.
      hstar(7) = 330.
      hstar(8) = 500.
      hstar(9) = 700.
#  else
      stop 'setup needed in setembm.F for ',ncat,' ice categories'
#  endif

      hstar(0) = 10.
      hstar(ncat) = 200000. !should not be used, make it real big anyway

      nsum = 0.
      do n=1,ncat
        nsum = nsum + nilay(n)
      enddo
      if (nsum .ne. ntilay) stop 'the sum of nilay must be ntilay'

      print*, 'cpts ice model set up with ncat=',ncat,
     &  ' ice categories and 1 open water category'
      print*, '   category intervals are: ',(hstar(n), n=0, ncat)
      print*, '   layers per category are:',(nilay(n), n=1, ncat)

!     minimum allowable fract
      asmall(0) = a0small
      do n=1,ncat
        asmall(n) = aismall
      enddo

!     matrix used to assist in heat transf. from cat i to j
      do n=1,ncat                     ! if nilay = { 2,4,8 }
        do m=1,n
          ncrel(n,m) = 1              !    ncrel = | 1 2 4 |
        enddo                         !            | 1 1 2 |
        do m=n,ncat                   !            | 1 1 1 |
          ncrel(n,m) = nilay(m)/nilay(n)
        enddo
      enddo
!     vectors identifying first and last layer in each bin
      layer1(1) = 1                   ! if nilay = { 2,4,8 }
      layern(1) = nilay(1)            !   layer1 = { 1,3,7 }
      do n=2,ncat                     !   layern = { 2,6,16}
        layer1(n) = layern(n-1) + 1
        layern(n) = layern(n-1) + nilay(n)
      enddo

!     default ridg matrices,  comp. all rows (rowflg=true)
!     assume ice thickness is mean of range for n < ncat
!     and 1m thicker than lower limit for ncat
      do n=1,ncat
        do k=1,ncat
          M_def(n,k) = 0.
          N_def(n,k) = 0.
          HN_def(n,k) = 0.
        enddo
        rowflg(n) = .true.
        if (n .lt. ncat) then
          Hi(n) = 0.5*(hstar(n-1) + hstar(n))
        else
           Hi(ncat) = hstar(ncat-1) + 1.*centi
        endif
      enddo
      call comp_matrices (rowflg, Hi, Hmean, M_def, N_def, HN_def)

!     setup the salinity profile and the melting temperature
!     for each layer
      salnew = 5.
!      saltmax = 3.2
      saltmax = 5.

      do n=1,ncat
        do k=1,nilay(n)
          zrel = (k-0.5)/nilay(n)
          saltz(k,n) = saltmax*0.5*(1.+sin(3.14159*(
     &                     zrel**(0.40706205/(zrel+0.57265966))-0.5)))
        enddo
        saltz(nilay(n)+1,n) = saltmax
        do k=1,nilay(n)
          tmelz(k,n) = -saltz(k,n)*alpha
        enddo
        print*, 'Category ',n
        write(*,'(A17,10(1x,f8.3))')
     &    '   salt profile:',(saltz(k,n),k=1,nilay(n)+1)
        write(*,'(A17,10(1x,f8.3))')
     &    '   melt temp:   ',(tmelz(k,n),k=1,nilay(n))
      enddo

! it would be nice to do this in a parameter statement
! because these are not variable
      rflice = flice*rhoice    ! specific latent heat of fusion ice
      rflsno = flice*rhosno    ! specific latent heat of fusion snow
      rslice = slice*rhoice    ! specific latent heat of sublim ice
      rslsno = slice*rhosno    ! specific latent heat of sublim snow

      rvlice = vlocn*rhoice    ! specific latent heat of vapour*rhoice
      rvlsno = vlocn*rhosno    ! specific latent heat of vapour*rhosno
      rcpice = cpice*rhoice    ! specific heat capacity of fresh ice
      rcpsno = cpsno*rhosno    ! specific heat capacity of snow

      rcpatm = rhoatm*cpatm
      rvlatm = rhoatm*vlocn
      rslatm = rhoatm*slice

      gamma  = rflice*ALPHA    ! heat capacity C=Cpi+gamma*salinity/T**2

# endif
!-----------------------------------------------------------------------
!     set initial conditions or read a restart
!-----------------------------------------------------------------------

      newcoef(:,:) = .true.

      nats = namix
      dayoyr = 1.
      itt = 0
      irstdy = 0
      msrsdy = 0
# if defined uvic_embm_running_average || defined uvic_embm_awind
      totaltime = 0.
      atbar(:,:) = 0.
      rtbar(:,:) = 0.
# endif
      at(:,:,:,:) = 0.
      tair = 13.
      at(:,:,:,isat) = tair
      ssh = cssh*exp(17.67*tair/(tair + 243.5))
      rh(:,:) = rhmax
      at(:,:,:,ishum) = rhmax*ssh
# if defined uvic_carbon && defined uvic_carbon_co2_2d
      at(:,:,:,ico2) = co2ccn
# endif
      precip(:,:) = 0.
# if !defined uvic_mom
      sbc(:,:,isst) = 0.
      sbc(:,:,isss) = 0.
# endif
# if defined uvic_embm_icedata
      aicel(:,:,:) = 0.
      hicel(:,:,:) = 0.
# endif
# if defined uvic_embm_awind
      awx(:,:) = 0.
      awy(:,:) = 0.
# endif
      soilm(:,:,:) = 0.
      surf(:,:) = 0.
# if defined uvic_ice
      hice(:,:,:) = 0.
      aice(:,:,:) = 0.
      tice(:,:) = 0.
      hsno(:,:,:) = 0.
# endif
# if defined uvic_ice_cpts
      hseff(:,:,:,:) = 0.
      A(:,:,:,:) = 0.
      heff(:,:,:,:) = 0.
      ts(:,:,:,:) = 0.
      E(:,:,:,:) = 0.
# endif
# if defined uvic_ice_cpts_roth_press && defined uvic_ice_cpts
        strength(:,:,:) = 0
# endif
# if defined uvic_ice_evp
      uice(:,:) = 0.
      vice(:,:) = 0.
      sbc(:,:,isu) = 0.
      sbc(:,:,isv) = 0.
      sbc(:,:,igu) = 0.
      sbc(:,:,igv) = 0.
      sig11n(:,:) = 0.
      sig11e(:,:) = 0.
      sig11s(:,:) = 0.
      sig11w(:,:) = 0.
      sig22n(:,:) = 0.
      sig22e(:,:) = 0.
      sig22s(:,:) = 0.
      sig22w(:,:) = 0.
      sig12n(:,:) = 0.
      sig12e(:,:) = 0.
      sig12s(:,:) = 0.
      sig12w(:,:) = 0.
# endif
# if !defined uvic_embm_explicit
      bv(:) = 0.
      xv(:) = 0.
# endif
# if defined uvic_embm_slap
      raux(:) = 0.
      iaux(:) = 0
# endif
# if defined uvic_embm_essl
      aux1(:,:,:) = 0.
      aux2(:) = 0.
      rparm(1) = 0.
      iparm(1) = 0
# endif
# if defined uvic_embm_sparskit
      work(:,:,:) = 0.
# endif
# if defined uvic_embm_icedata

!-----------------------------------------------------------------------
!     set land ice data and tracer grid ocean mask
!-----------------------------------------------------------------------

      call icedata
# else

!-----------------------------------------------------------------------
!     set tracer grid ocean mask
!-----------------------------------------------------------------------

      fname = new_file_name ("tmsk.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        call openfile (fname, iou)
        exists = inqvardef('tmsk', iou)
        ib(:) = 1
        ic(:) = imt
        ic(2) = jmt
        if (exists) call getvara ('tmsk', iou, imt*jmt, ib, ic, tmsk
     &,   c1, c0)
        call closefile (iou)
      endif
      if (.not. exists) then
        tmsk(:,:) = 0.
        do j=1,jmtm1
          do i=1,imtm1
            if (kmt(i,j) .gt. 0.) tmsk(i,j) = 1.
          enddo
        enddo
      endif
# endif

      if (.not. init) then
        fname = new_file_name ("restart_embm.nc")
        call embm_rest_in (fname, is, ie, js, je)
      endif
# if !defined uvic_mom

!-----------------------------------------------------------------------
!     initialize the time manager with specified initial conditions
!     time, user reference time, model time, and how long to integrate.
!-----------------------------------------------------------------------

      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,            ryear, rmonth, rday, rhour, rmin, rsec
     &,            irstdy, msrsdy, runlen, rununits, rundays, dtatm)
# endif
      flux(:,:,:) = 0.0

!-----------------------------------------------------------------------
!     read any atmospheric forcing data
!-----------------------------------------------------------------------

      call setatm

# if defined uvic_embm_awind
!-----------------------------------------------------------------------
!     read average air temperature
!-----------------------------------------------------------------------

      tbar(:,:) = 0.
      fname = new_file_name ("tbar.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'tbar in setembm.f'
      endif
      ib(:) = 1
      ic(:) = 1
      ic(1) = imt-2
      ic(2) = jmt-2
      call openfile (fname, iou)
      exists = inqvardef('slat', iou)
      if (.not. exists) then
        print*, "Error => slat does not exist."
        stop 'slat in setembm.f'
      endif
      call getvara ('slat', iou, ic(1)*ic(2), ib, ic, tmpij, c1, -C2K)
      call closefile (iou)
      tbar(2:imt-1,2:jmt-1) = tmpij(2:imt-1,2:jmt-1)
      call embmbc (tbar)

# endif
!-----------------------------------------------------------------------
!     read land elevations
!-----------------------------------------------------------------------

      elev(:,:) = 0.
      fname = new_file_name ("elev.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        ib(:) = 1
        ic(:) = imt
        ic(2) = jmt
        call openfile (fname, iou)
        call getvara ('elev', iou, imt*jmt, ib, ic, elev, c100, c0)
        call closefile (iou)
      endif
!     check for negative elevations
      where (elev(:,:) .lt. 0.) elev(:,:) = 0.

# if defined uvic_embm_running_average || defined uvic_embm_awind
!-----------------------------------------------------------------------
!     initialize running annual averages
!-----------------------------------------------------------------------

      fname = new_file_name ("restart_embm.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        call openfile (fname, iou)
        exists = inqvardef('rtbar', iou)
        call closefile (iou)
      endif
      if (.not. exists .or. init) then
        totaltime = 0.
        atbar(:,:) = 0.
        rtbar(:,:) = tbar(:,:)
      endif

# endif
!-----------------------------------------------------------------------
!     set velocity grid ocean mask
!-----------------------------------------------------------------------

      umsk(:,:) = 0.
      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (tmsk(i,j), tmsk(i+1,j), tmsk(i,j+1)
     &,                    tmsk(i+1,j+1))
        enddo
      enddo
      call embmbc (umsk)
!     remove isolated bays
      do j=2,jmtm1
        do i=2,imtm1
          tmsk(i,j) = max (umsk(i,j), umsk(i-1,j), umsk(i,j-1)
     &,                    umsk(i-1,j-1))
        enddo
      enddo
      call embmbc (tmsk)
      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (tmsk(i,j), tmsk(i+1,j), tmsk(i,j+1)
     &,                    tmsk(i+1,j+1))
        enddo
      enddo
      call embmbc (umsk)

#  if uvic_ice && !defined uvic_ice_cpts
      do j=1,jmt
        do i=1,imt
          if (tmsk(i,j) .ge. 0.5) then
            if (hice(i,j,1) .le. 0.) aice(i,j,1) = 0.
            if (hice(i,j,2) .le. 0.) aice(i,j,2) = 0.
          endif
        enddo
      enddo
      call embmbc (aice(1,1,1))
      call embmbc (aice(1,1,2))

#  endif
!-----------------------------------------------------------------------
!     set the river model
!-----------------------------------------------------------------------

      call rivinit

!-----------------------------------------------------------------------
!     set ocean coalbedo
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          if (kmt(i,j) .gt. 0) then
!           varies from 0.895 at the equator to 0.815 at the pole
            sbc(i,j,isca) = 0.855 + 0.04*cos(abs(tlat(i,j))*2./radian)
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     read vegetation class
!-----------------------------------------------------------------------

      rveg(:,:) = 0.
# if defined uvic_embm_cropdata
      fname = new_file_name ("veg_type.nc")
# else
      fname = new_file_name ("veg_class.nc")
# endif
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        call openfile (fname, iou)
        call getvara ('veg', iou, imt*jmt, ib, ic, rveg, c1, c0)
        call closefile (iou)
      endif
      do j=1,jmt
        do i=1,imt
          if (kmt(i,j) .gt. 0) then

            iveg(i,j) = iice
          elseif (rveg(i,j) .gt. 0.6 .and. rveg(i,j) .lt. 7.4) then
            iveg(i,j) = nint(rveg(i,j))
          else
            print*, "==> Warning: invalid vegetation type, ", rveg(i,j)
     &,             " at point (i,j), ", i, j
            print*, "    Setting point to vegetation type ice"
            iveg(i,j) = iice
          endif
        enddo
      enddo

      call gvsbc

# if defined uvic_ice_evp
!----------------------------------------------------------------------
!     initialize elastic viscous plastic variables
!-----------------------------------------------------------------------

      dlam = dxu(int(imt/2))/100.
      dphi = dyu(int(jmt/2))/100.
      diff1 = 0.004
      diff1 = diff1*dlam
      diff2 = diff1*dlam**2
      eccice = 2.
      ecc2 = 1./(eccice**2)
      ecc2m = 2.*(1.-ecc2)
      ecc2p = (1.+ecc2)
      zetamin = 4.e11
      eyc = 0.25
      dte = dtatm/float(ndte)
      dtei = 1./dte
      floor = 1.e-11
      do j=2,jmtm1
        do i=2,imtm1
           xyminevp=(min(cst(j)*dxt(i),dyt(j)))**2
        enddo
      enddo

# endif
# if defined uvic_embm_even_fluxes && defined uvic_mom
!-----------------------------------------------------------------------
!     check for even flux averaging
!-----------------------------------------------------------------------

      if (mod(namix,2) .ne. mod(nint(segtim*daylen/dtatm),2)) then
        write(*,*) '==> Error: time steps between mixing and coupling'
        write(*,*) '           in the atmosphere must both be even to'
        write(*,*) '           use uvic_embm_even_fluxes.'
        stop '=>setembm'
      endif

      if (mod(nmix,2) .ne. mod(nint(segtim*daylen/dtocn),2)) then
        write(*,*) '==> Error: time steps between mixing and coupling'
        write(*,*) '           in the ocean must both be even to use'
        write(*,*) '           uvic_embm_even_fluxes.'
        stop '=>setembm'
      endif

      if (mod(namix,2) .ne. mod(nats,2)) then
        write(*,*) '==> Warning: restart was not saved with even flux'
        write(*,*) '             averaging. Starting with a mixing time'
        write(*,*) '             step.'
        nats = namix
      endif

# endif
!-----------------------------------------------------------------------
!     check ice velocity calculation
!-----------------------------------------------------------------------

      if (nivts .gt. nint(segtim*daylen/dtatm)) then
        write(*,*) '==> Warning: ice velocities will be calculated'
        write(*,*) '             every coupling time.'
        nivts =  nint(segtim*daylen/dtatm)
      endif
# if defined uvic_embm_solve2y
!-----------------------------------------------------------------------
!     check atmosphere is even size in y
!-----------------------------------------------------------------------

      if (mod(jmtm2,2) .ne. 0) then
        write(*,*) '==> Error: atmosphere must be even sized in'
        write(*,*) '           latitude to use uvic_embm_solve2y.'
        stop '=>setembm'
      endif
# endif
# if defined uvic_embm_solve2x
!-----------------------------------------------------------------------
!     check atmosphere is even size in x
!-----------------------------------------------------------------------

      if (mod(imtm2,2) .ne. 0) then
        write(*,*) '==> Error: atmosphere must be even sized in'
        write(*,*) '           latitude to use uvic_embm_solve2x.'
        stop '=>setembm'
      endif
# endif
# if defined time_averages
!-----------------------------------------------------------------------
!     zero time average accumulators
!-----------------------------------------------------------------------

      call ta_embm_snap (is, ie, js, je, 0)

# endif
# if defined time_step_monitor
!-----------------------------------------------------------------------
!     zero integrated time average accumulators
!-----------------------------------------------------------------------

      call ta_embm_tsi (0)

# endif
#endif
      return
      end
