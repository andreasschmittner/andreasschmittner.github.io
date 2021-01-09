      program rot_nc
!-----------------------------------------------------------------------
!     Interpolates fields from one netcdf file into another.

!     Variables must be defined in both the input and output files (fi
!     and fo) for at least one time slice. Fields are interpolated from
!     the input file to the output file so the original contents of the
!     output file are lost. If times between the files differ new time
!     slices will be created in the output file. Any variables not 
!     defined in both the input and output files are ignored. Both 
!     components of vector fields are required. If the grids are rotated
!     relative to each other, the rotation must be defined with the 
!     variables psi, theta and phi.

!     If a variable is not being interpolated, make sure it is defined
!     in one of s2d, v2d, s3d or v3d. New or missing variables will have
!     to be added

!     Land or ocean masking may also be set for all variables (see ms2d,
!     mv2d, ms3d and mv3d). Masks may be read from a file (fmski, fmsko) 
!     or if the masking variable (vmski, vmsko) is defined, taken from 
!     the input or output file. The output mask may also be interpolated
!     from the input mask (see intrp_mask). If the input or output mask
!     is undefined, no masking is done.
!-----------------------------------------------------------------------

      implicit none

      integer maxs2d, maxv2d, maxs3d, maxv3d, max3d, kmax
      parameter (maxs2d=1000, maxv2d=100, maxs3d=500, maxv3d=50)
      parameter (max3d=5)
      character (120) :: s2d(maxs2d), v2d(maxv2d,2), s3d(maxs3d,max3d)
      character (120) :: v3d(maxv3d,2,max3d), fi, fo, fmski, fmsko
      character (120) :: vmski, vmsko, name3d(max3d)
      character(3) :: a3

      logical verbose, exists, masking, intrp_mask, inqvardef

      integer i, iou, ii, im, io, is, j, ji, jm, jo, js, k, ka, kb
      integer l, m, mlev, n, ntrec, ns2d, nv2d, ns3d, nv3d 
      integer ki(max3d), ko(max3d), i3d(max3d)
      integer :: is2d(maxs2d), iv2d(maxv2d), is3d(maxs3d,max3d)
      integer :: iv3d(maxv3d,max3d), ms2d(maxs2d), mv2d(maxv2d)
      integer :: ms3d(maxs3d,max3d),  mv3d(maxs3d,max3d)

      real :: fs2d(maxs2d), fv2d(maxv2d), fs3d(maxs3d,max3d)
      real :: fv3d(maxv3d,max3d), es2d(maxs2d),  ev2d(maxv2d)
      real :: es3d(maxs3d,max3d), ev3d(maxv3d,max3d)
      real time, psi, theta, phi, valmask, wtb
      real, allocatable :: dsi(:,:), dvi(:,:,:), dso(:,:), dvo(:,:,:)
      real, allocatable :: tmps1(:,:), tmps2(:,:)
      real, allocatable :: tmpv1(:,:,:), tmpv2(:,:,:)
      real, allocatable :: xsi(:), xvi(:), xso(:), xvo(:)
      real, allocatable :: ysi(:), yvi(:), yso(:), yvo(:)
      real, allocatable :: zsi(:), zso(:)
      real, allocatable :: xm(:), ym(:)
      integer, allocatable :: msi(:,:), mvi(:,:), mso(:,:), mvo(:,:)
      
      include 'netcdf.inc'

!-----------------------------------------------------------------------
!     initialise some stuff
!-----------------------------------------------------------------------

      psi   = 0.
      theta = 0.
      phi   = 0.
      
!     define rotation in radians
!     rotate into "rotated" model
!      phi   = -2.268928028
!      theta = -0.226892803
!      psi   = 3.141592654

!     rotate out of "rotated" model
!      psi   = 2.268928028
!      theta = 0.226892803
!      phi   = -3.141592654

!     set verbose to true to see field names as they are processed
      verbose = .true.
!     set intrp_mask to true to interpolate the mask from the input file
      intrp_mask = .false.
!     set masking to false if masking is not required
      masking = .true.

!     set file names
      fi = "in.nc"
      fo = "out.nc"
      fmski = "maski.nc"
      fmsko = "masko.nc"

!     set initial 2d and 3d vector and scalar fields to blank
      s2d(:) = " "
      v2d(:,:) = " "
      s3d(:,:) = " "
      v3d(:,:,:) = " "
!     set default masking (0=no, 1=land, -1=ocean)
      ms2d(:) = 0
      mv2d(:) = 0
      ms3d(:,:) = 0
      ms3d(:,1:3) = 1
      ms3d(:,4:5) = -1
      mv3d(:,:) = 0
      mv3d(:,1:3) = 1
      mv3d(:,4:5) = -1
!     set default interpolation (0=bilinear, 1=nearest)
      is2d(:) = 0
      iv2d(:) = 0
      is3d(:,:) = 0
      iv3d(:,:) = 0
!     set interpolation for 3rd dimension (0=no, 1=yes)
      i3d(1:2) = 1 !set only for ocean depths
      i3d(3:max3d) = 0
!     set default fill values (used in masking)
      fs2d(:) = nf_fill_float
      fv2d(:) = nf_fill_float
      fs3d(:,:) = nf_fill_float
      fv3d(:,:) = nf_fill_float
!     set default extapolation value (used in extrapolating values)
      es2d(:) = nf_fill_float
      ev2d(:) = nf_fill_float
      es3d(:,:) = nf_fill_float
      ev3d(:,:) = nf_fill_float
!     set default mask variable
      vmski = "kmt"
      vmsko = "kmt"
      time = 0.
!     mask anything over valmask
      valmask = 1.e15

!-----------------------------------------------------------------------
!     set the names of 2d and 3d scalars and vectors to be interpolated.
!     these may be commented out if not desired. if not found in the 
!     input and output file they are ignored
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     2d scalar definitions
!-----------------------------------------------------------------------
      m = 0
      call def_s2d (s2d, maxs2d, m, "kmt")
        is2d(m) = 1
        ms2d(m) = 1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "mskhr")
        is2d(m) = 1
        ms2d(m) = 1
        fs2d(m) = 0
      call def_s2d (s2d, maxs2d, m, "tlat")
      call def_s2d (s2d, maxs2d, m, "tlon")
      call def_s2d (s2d, maxs2d, m, "ulat")
      call def_s2d (s2d, maxs2d, m, "ulon")
      call def_s2d (s2d, maxs2d, m, "elev")
        ms2d(m) = -1
        fs2d(m) = 0
      call def_s2d (s2d, maxs2d, m, "sat")
      call def_s2d (s2d, maxs2d, m, "shum")
      call def_s2d (s2d, maxs2d, m, "slat")
      call def_s2d (s2d, maxs2d, m, "precip")
      call def_s2d (s2d, maxs2d, m, "evap")
      call def_s2d (s2d, maxs2d, m, "outlwr")
      call def_s2d (s2d, maxs2d, m, "uplwr")
      call def_s2d (s2d, maxs2d, m, "upsens")
      call def_s2d (s2d, maxs2d, m, "dnswr")
      call def_s2d (s2d, maxs2d, m, "upltnt")
      call def_s2d (s2d, maxs2d, m, "p_alb")
      call def_s2d (s2d, maxs2d, m, "a_alb")
      call def_s2d (s2d, maxs2d, m, "sf_alb")
      call def_s2d (s2d, maxs2d, m, "a_calb")
      call def_s2d (s2d, maxs2d, m, "s_alb")
      call def_s2d (s2d, maxs2d, m, "tice")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "hice")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "aice")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "hsno")
      call def_s2d (s2d, maxs2d, m, "soilm")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "runoff")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "surf")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "flxadj_t")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flxadj_s")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "psno")
      call def_s2d (s2d, maxs2d, m, "solins")
      call def_s2d (s2d, maxs2d, m, "ws")
      call def_s2d (s2d, maxs2d, m, "tmsk")
        is2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "aicel")
      call def_s2d (s2d, maxs2d, m, "hicel")      
      call def_s2d (s2d, maxs2d, m, "nriv")
        is2d(m) = 1
        ms2d(m) = -1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "rivers")
        is2d(m) = 1
        ms2d(m) = -1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "discharge")
        is2d(m) = 1
        ms2d(m) = 1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "weights")
        is2d(m) = 1
        ms2d(m) = 1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "veg")
        is2d(m) = 1
        ms2d(m) = -1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "vegtype")
        is2d(m) = 1
        ms2d(m) = -1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "disch")
        ms2d(m) = 1
        fs2d(m) = 0.
      call def_s2d (s2d, maxs2d, m, "co2")
      call def_s2d (s2d, maxs2d, m, "rtbar")
      call def_s2d (s2d, maxs2d, m, "apress")
      call def_s2d (s2d, maxs2d, m, "ws")
      call def_s2d (s2d, maxs2d, m, "wa")
      call def_s2d (s2d, maxs2d, m, "rh")
      call def_s2d (s2d, maxs2d, m, "rhum")
      call def_s2d (s2d, maxs2d, m, "sst")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sss")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "hflx")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sflx")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ssc")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sso2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ssa")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "atbar")
      call def_s2d (s2d, maxs2d, m, "avgp")
      call def_s2d (s2d, maxs2d, m, "accp")
      call def_s2d (s2d, maxs2d, m, "soilm1")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "soilm2")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "hice2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "aice1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "aice2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "hsno1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "hsno2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "strength1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "strength2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig11n")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig11s")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig11e")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig11w")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig22n")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig22s")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig22e")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig22w")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig12n")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig12s")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig12e")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "sig12w")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "TS1")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "TSTAR_GB")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "ALBLAND")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "ET")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "CS")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "RESP_S")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "land_map")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "TSOIL")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "LYING_SNOW")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "CV")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "VEG_FRAC")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "FRAC_VS")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "M")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "FSMC")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "RESP_S_DR")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "MNEG")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "sbc_sca")
      call def_s2d (s2d, maxs2d, m, "sbc_lwr")
      call def_s2d (s2d, maxs2d, m, "sbc_sens")
      call def_s2d (s2d, maxs2d, m, "sbc_evap")
      call def_s2d (s2d, maxs2d, m, "sbc_npp")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "sbc_sr")
        ms2d(m) = -1
      call def_s2d (s2d, maxs2d, m, "dissipation")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flux_heat")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flux_salt")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ps")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ps1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ps2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "psi")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "psi1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "psi2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ubarm1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ubarm2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ubar1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ubar2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ptd1")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "ptd2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flux_carbon")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flux_o2")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flux_c14")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flux_cfc11")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "flux_cfc12")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "totalk")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "vdepth")
        ms2d(m) = 1
      call def_s2d (s2d, maxs2d, m, "pe")
        ms2d(m) = 1
!     assume number of unnamed tracers, pfts, etc. =< 20
      do n=1,20
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        call def_s2d (s2d, maxs2d, m, "at_"//trim(a3))
        call def_s2d (s2d, maxs2d, m, "dn_"//trim(a3))
        call def_s2d (s2d, maxs2d, m, "de_"//trim(a3))
        call def_s2d (s2d, maxs2d, m, "at1_"//trim(a3))
        call def_s2d (s2d, maxs2d, m, "at2_"//trim(a3))
        call def_s2d (s2d, maxs2d, m, "TSTAR_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "ALBSNF_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "ALBSNC_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "HT_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "LAI_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "CVEG_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "G_LEAF_PHEN_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "G_LEAF_DR_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "NPP_DR_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "RESP_W_DR_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "CATCH_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "Z0_")
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "FRAC_"//trim(a3))
        ms2d(m) = -1
        call def_s2d (s2d, maxs2d, m, "hseff1_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "hseff2_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "heff1_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "heff2_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "A1_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "A2_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "ts1_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "ts2_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "E1_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "E2_"//trim(a3))
        ms2d(m) = 1
        call def_s2d (s2d, maxs2d, m, "flux_"//trim(a3))
        ms2d(m) = 1
      enddo

!-----------------------------------------------------------------------
!     2d vector definitions
!-----------------------------------------------------------------------
      m = 0
      call def_v2d (v2d, maxv2d, m, "taux","tauy")
        mv2d(m) = 1
      call def_v2d (v2d, maxv2d, m, "wx","wy")
      call def_v2d (v2d, maxv2d, m, "wx_q","wy_q")
      call def_v2d (v2d, maxv2d, m, "wx_t","wy_t")
      call def_v2d (v2d, maxv2d, m, "awx","awy")
      call def_v2d (v2d, maxv2d, m, "uice","vice")
        mv2d(m) = 1
      call def_v2d (v2d, maxv2d, m, "xint","yint")
        mv2d(m) = 1
      call def_v2d (v2d, maxv2d, m, "su","sv")
        mv2d(m) = 1
      call def_v2d (v2d, maxv2d, m, "gu","gv")
        mv2d(m) = 1

!     3d scalar in ocean depth
      name3d(1) = "ocean depth"
      m = 0
      call def_s3d (s3d, maxs3d, max3d, m, 1, "temp")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "temperature")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "salinity")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "dic")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "anthco2")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "tco2")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "alk")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "palk")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "silicate")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "phosphate")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "o2")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "oxygen")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "aou")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "o2sat")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "nutr")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "nitrate")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "chlorophyll")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "no3")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "c14")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "bkgc14")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "bombc14")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "cfc11")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "cfc12")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "pcfc11")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "pcfc12")
      
!     assume number of unnamed tracers, pfts, etc. =< 20
      do n=1,20
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        call def_s3d (s3d, maxs3d, max3d, m, 1, "tracer_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "tracer1_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "tracer2_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "tracer_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "d_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "-xadv_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "-yadv_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "-zadv_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "xdiff_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "ydiff_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "zdiff_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "source_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "convect_"//trim(a3))
        call def_s3d (s3d, maxs3d, max3d, m, 1, "filter_"//trim(a3))
      enddo
      call def_s3d (s3d, maxs3d, max3d, m, 1, "w")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "adv_vbtiso")
      call def_s3d (s3d, maxs3d, max3d, m, 1, "kv")

!-----------------------------------------------------------------------
!     3d scalar definitions in limited ocean depth
!-----------------------------------------------------------------------
      name3d(2) = "limited ocean depth"
      m = 0
      call def_s3d (s3d, maxs3d, max3d, m, 2, "phyt")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "zoop")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "detr")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "diaz")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "onpp")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "graz")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "morp")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "morpt")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "morz")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "remi")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "excr")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "expo")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "npp_D")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "graz_D")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "morp_D")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "phyt1")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "zoop1")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "detr1")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "phyt2")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "zoop2")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "detr2")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "diaz1")
      call def_s3d (s3d, maxs3d, max3d, m, 2, "diaz2")

!-----------------------------------------------------------------------
!     3d scalar definitions in ice category
!-----------------------------------------------------------------------
      name3d(3) = "ice category"
      m = 0
      call def_s3d (s3d, maxs3d, max3d, m, 3, "ticen")
      call def_s3d (s3d, maxs3d, max3d, m, 3, "hicen")
      call def_s3d (s3d, maxs3d, max3d, m, 3, "aicen")
      call def_s3d (s3d, maxs3d, max3d, m, 3, "hsnon")

!-----------------------------------------------------------------------
!     3d scalar definitions in plant funcional type
!-----------------------------------------------------------------------
      name3d(4) = "plant funcional type"
      m = 0
      call def_s3d (s3d, maxs3d, max3d, m, 4, "GPP")
      call def_s3d (s3d, maxs3d, max3d, m, 4, "NPP")
      call def_s3d (s3d, maxs3d, max3d, m, 4, "HT")
      call def_s3d (s3d, maxs3d, max3d, m, 4, "LAI")
      call def_s3d (s3d, maxs3d, max3d, m, 4, "C_VEG")

!-----------------------------------------------------------------------
!     3d scalar definitions in land type
!-----------------------------------------------------------------------
      name3d(5) = "land type"
      m = 0
      call def_s3d (s3d, maxs3d, max3d, m, 5, "FRAC")

!-----------------------------------------------------------------------
!     3d vector definitions in ocean depth
!-----------------------------------------------------------------------
      m = 0
      call def_v3d (v3d, maxv3d, max3d, m, 1, "u","v")
      call def_v3d (v3d, maxv3d, max3d, m, 1, "u1","v1")
      call def_v3d (v3d, maxv3d, max3d, m, 1, "u2","v2")
      call def_v3d (v3d, maxv3d, max3d, m, 1, "adv_vetiso","adv_vntiso")

!-----------------------------------------------------------------------
!     get input grid and allocate arrays
!-----------------------------------------------------------------------
      inquire (file=trim(fo), exist=exists)
      if (.not. exists) then
        print*, 'can not open input file: ',fi
        stop
      endif
      call openfile (fi, iou)
      ntrec = 1
      ii = 1
      ji = 1
      ki(:) = 0
      if (inqvardef ('time', iou)) call getdimlen ('time', iou, ntrec)
      if (inqvardef ('xu', iou)) call getdimlen ('xu', iou, ii)
      if (inqvardef ('yu', iou)) call getdimlen ('yu', iou, ji)
      if (inqvardef ('zw', iou)) call getdimlen ('zw', iou, ki(1))
      if (inqvardef ('xt', iou)) call getdimlen ('xt', iou, ii)
      if (inqvardef ('yt', iou)) call getdimlen ('yt', iou, ji)
      if (inqvardef ('zt', iou)) call getdimlen ('zt', iou, ki(1))
      ki(2) = ki(1)
      if (inqvardef ('zl', iou)) call getdimlen ('zl', iou, ki(2))
      if (inqvardef ('cat', iou)) call getdimlen ('cat', iou, ki(3))
      if (inqvardef ('pft', iou)) call getdimlen ('pft', iou, ki(4))
      if (inqvardef ('type', iou)) call getdimlen ('type', iou, ki(5))
      allocate ( xsi(ii) )
      allocate ( ysi(ji) )
      if (ki(1) .gt. 0) allocate ( zsi(ki(1)) )
      allocate ( dsi(ii,ji) )
      allocate ( msi(ii,ji) )
      allocate ( xvi(ii) )
      allocate ( yvi(ji) )
      allocate ( dvi(ii,ji,2) )
      allocate ( mvi(ii,ji) )
      if (inqvardef ('xt', iou))
     &  call getvara ('xt', iou, ii, (/1/), (/ii/), xsi, 1., 0.)
      if (inqvardef ('xu', iou))
     &  call getvara ('xu', iou, ii, (/1/), (/ii/), xvi, 1., 0.)
      if (inqvardef ('yt', iou))
     &  call getvara ('yt', iou, ji, (/1/), (/ji/), ysi, 1., 0.)
      if (inqvardef ('yu', iou))
     &  call getvara ('yu', iou, ji, (/1/), (/ji/), yvi, 1., 0.)
      if (inqvardef ('zt', iou))
     &  call getvara ('zt', iou, ii, (/1/), (/ki(1)/), zsi, 1., 0.)
      call closefile (iou)

!-----------------------------------------------------------------------
!     get output grid and allocate arrays
!-----------------------------------------------------------------------
      inquire (file=trim(fo), exist=exists)
      if (.not. exists) then
        print*, 'can not open ouput file: ',fo
        stop
      endif
      call openfile (fo, iou)
      io = 1
      jo = 1
      ko(:) = 0
      if (inqvardef ('xu', iou)) call getdimlen ('xu', iou, io)
      if (inqvardef ('yu', iou)) call getdimlen ('yu', iou, jo)
      if (inqvardef ('zw', iou)) call getdimlen ('zw', iou, ko(1))
      if (inqvardef ('xt', iou)) call getdimlen ('xt', iou, io)
      if (inqvardef ('yt', iou)) call getdimlen ('yt', iou, jo)
      if (inqvardef ('zt', iou)) call getdimlen ('zt', iou, ko(1))
      ki(2) = ki(1)
      if (inqvardef ('zl', iou)) call getdimlen ('zl', iou, ko(2))
      if (inqvardef ('cat', iou)) call getdimlen ('cat', iou, ko(3))
      if (inqvardef ('pft', iou)) call getdimlen ('pft', iou, ko(4))
      if (inqvardef ('type', iou)) call getdimlen ('type', iou, ko(5))
      allocate ( xso(io) )
      allocate ( yso(jo) )
      if (ko(1) .gt. 0) allocate ( zso(ko(1)) )
      allocate ( dso(io,jo) )
      allocate ( tmps1(io,jo) )
      allocate ( tmps2(io,jo) )
      allocate ( mso(io,jo) )
      allocate ( xvo(io) )
      allocate ( yvo(jo) )
      allocate ( dvo(io,jo,2) )
      allocate ( tmpv1(io,jo,2) )
      allocate ( tmpv2(io,jo,2) )
      allocate ( mvo(io,jo) )
      if (inqvardef ('xt', iou))
     &  call getvara ('xt', iou, io, (/1/), (/io/), xso, 1., 0.)
      if (inqvardef ('xu', iou))
     &  call getvara ('xu', iou, io, (/1/), (/io/), xvo, 1., 0.)
      if (inqvardef ('yt', iou))
     &  call getvara ('yt', iou, jo, (/1/), (/jo/), yso, 1., 0.)
      if (inqvardef ('yu', iou))
     &  call getvara ('yu', iou, jo, (/1/), (/jo/), yvo, 1., 0.)
      if (inqvardef ('zt', iou))
     &  call getvara ('zt', iou, io, (/1/), (/ko(1)/), zso, 1., 0.)
      call closefile (iou)

!-----------------------------------------------------------------------
!     get masks
!-----------------------------------------------------------------------

      if (masking) then
!       get input ocean mask from the mask or input file
        inquire (file=trim(fmski), exist=exists)
        if (exists) then
!         get input ocean mask from a mask file
          call openfile (fmski, iou)
        else
!         get input ocean mask from input file
          call openfile (fi, iou)
        endif        
        if (inqvardef (vmski, iou)) then
          im = 0
          jm = 0
          if (inqvardef ('xt', iou)) call getdimlen ('xt', iou, im)
          if (inqvardef ('yt', iou)) call getdimlen ('yt', iou, jm)
          if (im .gt. 0 .and. jm .gt. 0) then
            allocate ( xm(im) )
            allocate ( ym(jm) )
            call getvara ('xt', iou, im, (/1/), (/im/), xm, 1., 0.)
            call getvara ('yt', iou, jm, (/1/), (/jm/), ym, 1., 0.)
            is = 0
            do i=1,min(im, ii)
              if (abs(xm(i) - xsi(1)) .lt. 1.e-3) is = i
            enddo
            deallocate ( xm )
            js = 0
            do j=1,min(jm, ji)
              if (abs(ym(j) - ysi(1)) .lt. 1.e-3) js = j
            enddo
            deallocate ( ym )
            if (is .gt. 0 .and. js .gt. 0) then
              call getvara (vmski, iou, ii*ji, (/is,js/), (/ii,ji/)
     &,         dsi, 1., 0.)
            else
              masking = .false.
            endif 
          else
            masking = .false.
          endif
        else
          masking = .false.
        endif
        call closefile (iou)
      endif

      if (masking) then
!       get output ocean mask from the mask or output file
!       or interpolate it from input file mask
        inquire (file=trim(fmsko), exist=exists)
        if (exists .or. .not. intrp_mask) then
          if (exists) then
!           get output ocean mask from a mask file
            call openfile (fmsko, iou)
          else
!           get output ocean mask from output file
            call openfile (fo, iou)
          endif
          if (inqvardef (vmsko, iou)) then
            im = 0
            jm = 0
            if (inqvardef ('xt', iou)) call getdimlen ('xt', iou, im)
            if (inqvardef ('yt', iou)) call getdimlen ('yt', iou, jm)
            if (im .gt. 0 .and. jm .gt. 0) then
              allocate ( xm(im) )
              allocate ( ym(jm) )
              call getvara ('xt', iou, im, (/1/), (/im/), xm, 1., 0.)
              call getvara ('yt', iou, jm, (/1/), (/jm/), ym, 1., 0.)
              is = 0
              do i=1,min(im, io)
                if (abs(xm(i) - xso(1)) .lt. 1.e-3) is = i
              enddo
              deallocate ( xm )
              js = 0
              do j=1,min(jm, jo)
                if (abs(ym(j) - yso(1)) .lt. 1.e-3) js = j
              enddo
              deallocate ( ym )
              if (is .gt. 0 .and. js .gt. 0) then
                call getvara (vmsko, iou, io*jo, (/is,js/), (/io,jo/)
     &,           dso, 1., 0.)
              else
                masking = .false.
              endif 
            else
              masking = .false.
            endif
          else
            masking = .false.
          endif
          call closefile (iou)
        else
          call rot_intrp_sclr (dsi, xsi, ysi, ii, ji, dso, xso
     &,     yso, io, jo, phi, theta, psi, -abs(nf_fill_float), 1)
        endif
      endif

      if (masking) then
!       write mask to the output file
        call openfile (fo, iou)
        if (inqvardef (vmski, iou)) then
          call putvara (vmski, iou, io*jo, (/1,1/), (/io,jo/), dso
     &,     1., 0.)
        endif
        call closefile (iou)
!       calculate input velocity mask
        msi = dsi
        do i=1,ii-1
          do j=1,ji-1
            mvi(i,j) = min(msi(i,j),msi(i+1,j),msi(i,j+1),msi(i+1,j+1))
          enddo
        enddo
        mvi(ii,:) = mvi(2,:)
        mvi(:,ji) = 0 
!       calculate a output velocity mask
        mso = dso
        do i=1,io-1
          do j=1,jo-1
            mvo(i,j) = min(mso(i,j),mso(i+1,j),mso(i,j+1),mso(i+1,j+1))
          enddo
        enddo
        mvo(io,:) = mvo(2,:)
        mvo(:,jo) = 0
      else
!       if no masking set all masks to 1
        msi(:,:) = 1
        mvi(:,:) = 1
        mso(:,:) = 1
        mvo(:,:) = 1
      endif

!-----------------------------------------------------------------------
!     loop through all time records
!-----------------------------------------------------------------------
      do m=1,ntrec
        call openfile (fi, iou)
        if (inqvardef ('time', iou)) then
          call getvars ('time', iou, m, time, 1., 0.)    
          call closefile (iou)
          call openfile (fo, iou)
          if (inqvardef ('time', iou)) then
            call putvars ('time', iou, m, time, 1., 0.)    
          endif
        endif
        call closefile (iou)
        if (verbose) print*, 'time record: ', m, ' time: ', time

!-----------------------------------------------------------------------
!       rotate and interpolate 2d scalar data
!-----------------------------------------------------------------------
        do n=1,maxs2d
          if (trim(s2d(n)) .ne. vmski) then
            call openfile (fi, iou)
            if (inqvardef (s2d(n), iou)) then
              if (verbose) print*, '2d scalar: ',trim(s2d(n))
              call getvara (s2d(n), iou, ii*ji, (/1,1,m/)          
     &,         (/ii,ji,1/), dsi, 1., 0.)             
              if (ms2d(n) .eq. 1) then
                call msk_lnd (dsi(:,:), msi, es2d(n), ii, ji, 1)
              endif
              if (ms2d(n) .eq. -1) then
                call msk_ocn (dsi(:,:), msi, es2d(n), ii, ji, 1)
              endif
              call rot_intrp_sclr (dsi, xsi, ysi, ii, ji, dso, xso
     &,         yso, io, jo, phi, theta, psi, -abs(es2d(n)), is2d(n))
              if (ms2d(n) .ne. 0) then
                call extrap2 (dso(:,:), es2d(n), xso, io, jo)
              endif
              if (ms2d(n) .eq. 1) then
                call msk_lnd (dso(:,:), mso, fs2d(n), io, jo, 1)
              endif
              if (ms2d(n) .eq. -1) then
                call msk_ocn (dso(:,:), mso, fs2d(n), io, jo, 1)
              endif
            else
              dso(:,:) = fs2d(n)
            endif
            call closefile (iou)

            call openfile (fo, iou)
            if (inqvardef (s2d(n), iou)) then
              where (dso(:,:) .gt. valmask) dso(:,:) = fs2d(n)
              call putvara (s2d(n), iou, io*jo, (/1,1,m/)          
     &,         (/io,jo,1/), dso, 1., 0.)
            endif
            call closefile (iou)

          endif
        enddo

!-----------------------------------------------------------------------
!       rotate and interpolate 2d vector data
!-----------------------------------------------------------------------
        do n=1,maxv2d
          if (trim(v2d(n,1)) .ne. vmski) then
            call openfile (fi, iou)
            if (inqvardef (v2d(n,1), iou) .and. 
     &        inqvardef (v2d(n,2), iou)) then
              if (verbose) print*, '2d vector: ',trim(v2d(n,1)), ' '
     &,         trim(v2d(n,2))
              call getvara (v2d(n,1), iou, ii*ji, (/1,1,m/)          
     &,         (/ii,ji,1/), dvi(1,1,1), 1., 0.)    
              call getvara (v2d(n,2), iou, ii*ji, (/1,1,m/)          
     &,         (/ii,ji,1/), dvi(1,1,2), 1., 0.)    
              if (mv2d(n) .eq. 1) then
                call msk_lnd (dvi(:,:,1), mvi, ev2d(n), ii, ji, 1)
                call msk_lnd (dvi(:,:,2), mvi, ev2d(n), ii, ji, 1)
              endif
              if (mv2d(n) .eq. -1) then
                call msk_ocn (dvi(:,:,1), mvi, ev2d(n), ii, ji, 1)
                call msk_ocn (dvi(:,:,2), mvi, ev2d(n), ii, ji, 1)
              endif
              call rot_intrp_vctr (dvi, xvi, yvi, ii, ji, dvo, xvo
     &,         yvo, io, jo, phi, theta, psi, -abs(ev2d(n)), iv2d(n))
              if (mv2d(n) .ne. 0) then 
                call extrap2 (dvo(:,:,1), ev2d(n), xvo, io, jo)
                call extrap2 (dvo(:,:,2), ev2d(n), xvo, io, jo)
              endif
              if (mv2d(n) .eq. 1) then
                call msk_lnd (dvo(:,:,1), mvo, fv2d(n), io, jo, 1)
                call msk_lnd (dvo(:,:,2), mvo, fv2d(n), io, jo, 1)
              endif
              if (mv2d(n) .eq. -1) then
                call msk_ocn (dvo(:,:,1), mvo, fv2d(n), io, jo, 1)
                call msk_ocn (dvo(:,:,2), mvo, fv2d(n), io, jo, 1)
              endif
            else
              dvo(:,:,:) = fv2d(n)
            endif
            call closefile (iou)

            call openfile (fo, iou)
            if (inqvardef (v2d(n,1), iou)) then
              where (dvo(:,:,1) .gt. valmask) dvo(:,:,1) = fv2d(n)
              call putvara (v2d(n,1), iou, io*jo, (/1,1,m/)          
     &,         (/io,jo,1/), dvo(1,1,1), 1., 0.)
            endif   
            if (inqvardef (v2d(n,2), iou)) then
              where (dvo(:,:,2) .gt. valmask) dvo(:,:,2) = fv2d(n)
              call putvara (v2d(n,2), iou, io*jo, (/1,1,m/)          
     &,         (/io,jo,1/), dvo(1,1,2), 1., 0.)    
            endif   
            call closefile (iou)
          endif
        enddo

!-----------------------------------------------------------------------
!       loop through each type of 3d dimension
!-----------------------------------------------------------------------
        do i=1,max3d
!-----------------------------------------------------------------------
!         rotate and interpolate 3d scalar data
!-----------------------------------------------------------------------
          do n=1,maxs3d
          
            do k=1,ko(i)
              if (k .lt. 1000) write(a3, '(i3)') k
              if (k .lt. 100) write(a3, '(i2)') k
              if (k .lt. 10) write(a3, '(i1)') k
              ka = k
              kb = k
              mlev = 1
              if (i3d(i) .ne. 0) then
                ka = 1
                do j=2,ki(i)
                  if (zsi(j) .lt. zso(k)) ka = j
                enddo
                kb = min(ka+1,ki(i))
                ka = max(kb-1,1)
                wtb = min(1.0, max(0.0, (zso(k) - zsi(ka))
     &                /(zsi(kb) - zsi(ka))))
                mlev = ka
              endif

              if (trim(s3d(n,i)) .ne. vmski) then
                call openfile (fi, iou)
                if (inqvardef (s3d(n,i), iou)) then
                  if (verbose) print*, '3rd dim: ', trim(name3d(i))
     &,             '  level: ', a3, '  3d scalar: ', trim(s3d(n,i))
                  call getvara (s3d(n,i), iou, ii*ji, (/1,1,ka,m/)          
     &,            (/ii,ji,1,1/), dsi, 1., 0.)    
                  if (ms3d(n,i) .eq. 1) then
                    call msk_lnd (dsi(:,:), msi, es3d(n,i), ii, ji
     &,               mlev)
                  endif
                  if (ms3d(n,i) .eq. -1) then
                    call msk_ocn (dsi(:,:), msi, es3d(n,i), ii, ji
     &,               mlev)
                  endif
                  call rot_intrp_sclr (dsi, xsi, ysi, ii, ji, dso
     &,             xso, yso, io, jo, phi, theta, psi, -abs(es3d(n,i))
     &,             is3d(n,i))
     
                  if (i3d(i) .ne. 0) then
                    mlev = kb
                    tmps1(:,:) = dso(:,:)
                    call getvara (s3d(n,i), iou, ii*ji, (/1,1,kb,m/)          
     &,              (/ii,ji,1,1/), dsi, 1., 0.)    
                    if (ms3d(n,i) .eq. 1) then
                      call msk_lnd (dsi(:,:), msi, es3d(n,i), ii, ji
     &,                 mlev)
                    endif
                    if (ms3d(n,i) .eq. -1) then
                      call msk_ocn (dsi(:,:), msi, es3d(n,i), ii, ji
     &,                 mlev)
                    endif
                    call rot_intrp_sclr (dsi, xsi, ysi, ii, ji, dso
     &,               xso, yso, io, jo, phi, theta, psi
     &,               -abs(es3d(n,i)), is3d(n,i))
                    tmps2(:,:) = dso(:,:)
                    call intrp_vert (tmps1, tmps2, wtb, 
     &                -abs(es3d(n,i)), dso, io, jo)
                    mlev = k
                  endif
                  
                  if (ms3d(n,i) .ne. 0) then 
                    call extrap2 (dso(:,:), es3d(n,i), xso, io, jo)
                  endif
                  if (ms3d(n,i) .eq. 1) then
                    call msk_lnd (dso(:,:), mso, fs3d(n,i), io, jo
     &,               mlev)
                  endif
                  if (ms3d(n,i) .eq. -1) then
                    call msk_ocn (dso(:,:), mso, fs3d(n,i), io, jo
     &,               mlev)
                  endif
                else
                  dso(:,:) = fs3d(n,i)
                endif
                call closefile (iou)

                call openfile (fo, iou)
                if (inqvardef (s3d(n,i), iou)) then
                  where (dso(:,:) .gt. valmask) dso(:,:) = fs3d(n,i)
                  call putvara (s3d(n,i), iou, io*jo, (/1,1,k,m/)          
     &,             (/io,jo,1,1/), dso, 1., 0.)    
                endif
                call closefile (iou)
              endif
            enddo

          enddo

!-----------------------------------------------------------------------
!         rotate and interpolate 3d vector data
!-----------------------------------------------------------------------
          do n=1,maxv3d

            do k=1,ko(i)
              if (k .lt. 1000) write(a3, '(i3)') k
              if (k .lt. 100) write(a3, '(i2)') k
              if (k .lt. 10) write(a3, '(i1)') k
              ka = k
              kb = k
              mlev = 1
              if (i3d(i) .ne. 0) then
                 ka = 1
                do j=2,ki(i)
                  if (zsi(j) .lt. zso(k)) ka = j
                enddo
                kb = min(ka+1,ki(i))
                ka = max(kb-1,1)
                wtb = min(1.0, max(0.0, (zso(k) - zsi(ka))
     &                /(zsi(kb) - zsi(ka))))
                mlev = ka
              endif

              if (trim(v3d(n,1,i)) .ne. vmski) then
                call openfile (fi, iou)
                if (inqvardef (v3d(n,1,i), iou) .and.
     &            inqvardef (v3d(n,2,i), iou)) then
                  if (verbose) print*, '3rd dim: ', trim(name3d(i))
     &,             '  level: ', a3, '  3d vector: '
     &,             trim(v3d(n,1,i)), ' ', trim(v3d(n,2,i))
                  call getvara (v3d(n,1,i), iou, ii*ji, (/1,1,ka,m/)  
     &,             (/ii,ji,1,1/), dvi(1,1,1), 1., 0.)
                  call getvara (v3d(n,2,i), iou, ii*ji, (/1,1,ka,m/)          
     &,             (/ii,ji,1,1/), dvi(1,1,2), 1., 0.)    
                  if (mv3d(n,i) .eq. 1) then
                    call msk_lnd (dvi(:,:,1), mvi, ev3d(n,i), ii, ji
     &,               mlev)
                    call msk_lnd (dvi(:,:,2), mvi, ev3d(n,i), ii, ji
     &,               mlev)
                  endif
                  if (mv3d(n,i) .eq. -1) then
                    call msk_ocn (dvi(:,:,1), mvi, ev3d(n,i), ii, ji
     &,               mlev)
                    call msk_ocn (dvi(:,:,2), mvi, ev3d(n,i), ii, ji
     &,               mlev)
                  endif
                  call rot_intrp_vctr (dvi, xvi, yvi, ii, ji, dvo, xvo
     &,             yvo, io, jo, phi, theta, psi, -abs(ev3d(n,i))
     &,             iv3d(n,i))
     
                  if (i3d(i) .ne. 0) then
                    mlev = kb
                    tmpv1(:,:,:) = dvo(:,:,:)
                    call getvara (v3d(n,1,i), iou, ii*ji, (/1,1,kb,m/)  
     &,               (/ii,ji,1,1/), dvi(1,1,1), 1., 0.)
                    call getvara (v3d(n,2,i), iou, ii*ji, (/1,1,kb,m/)          
     &,               (/ii,ji,1,1/), dvi(1,1,2), 1., 0.)    
                    if (mv3d(n,i) .eq. 1) then
                      call msk_lnd (dvi(:,:,1), mvi, ev3d(n,i), ii, ji
     &,                 kb)
                      call msk_lnd (dvi(:,:,2), mvi, ev3d(n,i), ii, ji
     &,                 kb)
                    endif
                    if (mv3d(n,i) .eq. -1) then
                      call msk_ocn (dvi(:,:,1), mvi, ev3d(n,i), ii, ji
     &,                 kb)
                      call msk_ocn (dvi(:,:,2), mvi, ev3d(n,i), ii, ji
     &,                 kb)
                    endif
                    call rot_intrp_vctr (dvi, xvi, yvi, ii, ji, dvo
     &,               xvo, yvo, io, jo, phi, theta, psi
     &,               -abs(ev3d(n,i)), iv3d(n,i))
                    tmpv2(:,:,:) = dvo(:,:,:)
                    call intrp_vert (tmpv1(:,:,1), tmpv2(:,:,1), wtb
     &,              -abs(ev3d(n,i)), dvo(:,:,1), io, jo)
                    call intrp_vert (tmpv1(:,:,2), tmpv2(:,:,2), wtb
     &,              -abs(ev3d(n,i)), dvo(:,:,2), io, jo)
                    mlev = k
                  endif
                  
                  if (mv3d(n,i) .ne. 0) then 
                    call extrap2 (dvo(:,:,1), ev3d(n,i), xvo, io, jo)
                    call extrap2 (dvo(:,:,2), ev3d(n,i), xvo, io, jo)
                  endif
                  if (mv3d(n,i) .eq. 1) then
                    call msk_lnd (dvo(:,:,1), mvo, fv3d(n,i), io, jo
     &,               k)
                    call msk_lnd (dvo(:,:,2), mvo, fv3d(n,i), io, jo
     &,               k)
                  endif
                  if (mv3d(n,i) .eq. -1) then
                    call msk_ocn (dvo(:,:,1), mvo, fv3d(n,i), io, jo
     &,               k)
                    call msk_ocn (dvo(:,:,2), mvo, fv3d(n,i), io, jo
     &,               k)
                  endif
                else
                  dvo(:,:,:) = fv3d(n,i)
                endif
                call closefile (iou)              

                call openfile (fo, iou)
                if (inqvardef (v3d(n,1,i), iou)) then
                  where (dvo(:,:,1) .gt. valmask) dvo(:,:,1) = fv3d(n,i)
                  call putvara (v3d(n,1,i), iou, io*jo, (/1,1,k,m/)          
     &,             (/io,jo,1,1/), dvo(1,1,1), 1., 0.)    
                endif
                if (inqvardef (v3d(n,2,i), iou)) then
                  where (dvo(:,:,2) .gt. valmask) dvo(:,:,2) = fv3d(n,i)
                  call putvara (v3d(n,2,i), iou, io*jo, (/1,1,k,m/)          
     &,             (/io,jo,1,1/), dvo(1,1,2), 1., 0.)    
                endif
                call closefile (iou)              
              endif
            enddo

          enddo

        enddo

      enddo

      end
      
!-----------------------------------------------------------------------
      subroutine def_s2d (s2d, maxs2d, m, name)
!-----------------------------------------------------------------------
      
      integer maxs2d, m
      character(*) :: name
      character(120) :: s2d(maxs2d)
      
      m = m + 1
      if (m .gt. maxs2d) stop 'increase maxs2d'
      s2d(m) = trim(name)

      return
      end
      
!-----------------------------------------------------------------------
      subroutine def_v2d (v2d, maxv2d, m, name1, name2)
!-----------------------------------------------------------------------
      
      integer maxv2d, m
      character(*) :: name1, name2
      character(120) :: v2d(maxv2d,2)
      
      m = m + 1
      if (m .gt. maxv2d) stop 'increase maxv2d'
      v2d(m,1:2) = (/trim(name1),trim(name2)/)

      return
      end

!-----------------------------------------------------------------------
      subroutine def_s3d (s3d, maxs3d, max3d, m, n, name)
!-----------------------------------------------------------------------
      
      integer maxs3d, max3d, m, n
      character(*) :: name
      character(120) :: s3d(maxs3d,max3d)
      
      m = m + 1
      if (m .gt. maxs3d) stop 'increase maxs3d'
      s3d(m,n) = trim(name)

      return
      end
      
!-----------------------------------------------------------------------
      subroutine def_v3d (v3d, maxv3d, max3d, m, n, name1, name2)
!-----------------------------------------------------------------------
      
      integer maxv3d, max3d, m, n
      character(*) :: name1, name2
      character(120) :: v3d(maxv3d,2,max3d)
      
      m = m + 1
      if (m .gt. maxv3d) stop 'increase maxv3d'
      v3d(m,1:2,n) = (/trim(name1),trim(name2)/)

      return
      end

!-----------------------------------------------------------------------
      subroutine msk_ocn (data, mask, val, id, jd, lev)
!-----------------------------------------------------------------------
      
      integer id, jd, lev
      real data(id,jd), val
      integer mask(id,jd)
      
      where (mask(:,:) .ge. lev) data(:,:) = val

      return
      end

!-----------------------------------------------------------------------
      subroutine msk_lnd (data, mask, val, id, jd, lev)
!-----------------------------------------------------------------------
      
      integer id, jd, lev
      real data(id,jd), val
      integer mask(id,jd)
      
      where (mask(:,:) .lt. lev) data(:,:) = val

      return
      end

      
