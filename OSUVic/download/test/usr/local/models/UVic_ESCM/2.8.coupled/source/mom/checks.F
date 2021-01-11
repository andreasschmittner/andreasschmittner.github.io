      subroutine checks (errorc, vmixset, hmixset)
      logical errorc, vmixset, hmixset
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "hmixc.h"
#include "iounit.h"
#include "levind.h"
#if defined isopycmix || defined isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#include "vmixc.h"
#if defined xbts
# include "xbt.h"
#endif

!-----------------------------------------------------------------------
!     do consistency checks before allowing MOM to continue
!-----------------------------------------------------------------------

      write (stdout,'(/,20x,a,/)')
     &         'G E N E R A L    C O N S I S T E N C Y    C H E C K S'

      if (imt .lt. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: parameter "imt" less than 3 is not allowed         '
        errorc = .true.
      endif

      if (jmt .lt. 4) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: parameter "jmt" less than 4 is not allowed          '
        errorc = .true.
      endif
      if (num_processors .gt. 1 .and. num_processors .le. jmt-2) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "num_processors" is set for multitasking using    '
#if defined coarse_grained_parallelism
     &,'             option coarse_grained_parallelism (microtasking)  '
#else
     &,'             fine grained parallelism (autotasking)            '
#endif
      elseif (num_processors .eq. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "num_processors" = 1 for unitasking               '
      elseif (num_processors .gt. jmt-2 .or. num_processors .lt. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "num_processors" cannot be > jmt-2 or < 1           '
        write (stdout,*)' num_processors = ',num_processors
        write (stdout,*)
     &  'set num_processors using environment variable NCPUS in script '
        errorc = .true.
      endif

#if defined coarse_grained_parallelism
# if defined fourth_order_window && !defined pressure_gradient_average
      if (jmw .ne. 5) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: jmw must = 5 when "coarse_grained_parallelism"      '
     &,'           is enabled but not "pressure_gradient_average"      '
        errorc = .true.
      endif
# endif
# define other_fourth_order =  defined biharmonic || defined fourth_order_tracer_advection || defined fct
# if defined other_fourth_order && defined pressure_gradient_average
      if (jmw .ne. 7) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: jmw must = 7 when "coarse_grained_parallelism"      '
     &,'           is enabled but not "pressure_gradient_average"      '
        errorc = .true.
      endif
# endif
# if !defined fourth_order_window
      if (jmw .ne. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: jmw must = 3 when "coarse_grained_parallelism"      '
     &,'           is enabled.                                         '
        errorc = .true.
      endif
# endif
      if (ncrows .ne. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "ncrows" must = 1 when "coarse_grained_parallelism" '
     &,'            is enabled.                                        '
        errorc = .true.
      endif
# if defined equivalence_mw
      write (stdout,'(/,(1x,a))')
     & '==> Error: "equivalence_mw" is not allowed with option         '
     &,'           "coarse_grained_parallelism"                        '
      errorc = .true.
# endif
# if defined trajectories
      write (stdout,'(/,(1x,a))')
     & '==> Error: "trajectories" is not allowed with option           '
     &,'           "coarse_grained_parallelism"                        '
      errorc = .true.
# endif
#endif
#if defined pressure_gradient_average
      if (.not. eb) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: option "pressure_gradient_average" may need       '
     &,'             euler backward mixing timesteps to suppress the   '
     &,'             computational modes.                              '
      endif
#endif

#if defined linearized_advection
# if defined fct
      write (stdout,'(/,(1x,a))')
     & '==> Error: fct" is incompatable with "linearized_advection".   '
        errorc = .true.
# endif
# if defined fourth_order_tracer_advection
      write (stdout,'(/,(1x,a))')
     & '==> Error: "fourth_order_tracer_advection" is incompatable     '
     &,'            with "linearized_advection".                       '
        errorc = .true.
# endif
# if defined quicker
      write (stdout,'(/,(1x,a))')
     & '==> Error: "quicker" is incompatable with                      '
     &,'           "linearized_advection".                             '
        errorc = .true.
# endif
# if !defined linearized_density
      write (stdout,'(/,(1x,a))')
     & '==> Error: "linearized_density" must be enabled when using     '
     &,'           "linearized_advection".                             '
        errorc = .true.
# endif
# if defined minimize_sbc_memory && defined time_mean_sbc_data
      write (stdout,'(/,(1x,a))')
     & '==> Error: "minimize_sbc_memory" and "time_mean_sbc_data"      '
     &,'           are incompatible. "minimize_sbc_memory" is          '
     &,'           intended for option "time_varying_sbc_data".        '
        errorc = .true.
# endif
#endif

#if defined consthmix
      if (hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "consthmix"  cannot be enabled because another      '
     &,'            horizontal mixing scheme has been enabled          '
        errorc = .true.
      else
        hmixset = .true.
      endif
#endif
      if (.not.hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal mixing scheme has been enabled        '
        errorc = .true.
      endif

#if defined constvmix
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "constvmix"  cannot be enabled because another      '
     &,'            vertical mixing scheme has been enabled            '
        errorc = .true.
      else

!       set vmixset = true for enabeling "constvmix"

        vmixset = .true.
      endif
#endif
#if defined bryan_lewis_vertical
      if (.not. vmixset) then
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a/1x,a/1x,a))')
     & '==> Error: "bryan_lewis_vertical" tracer diffusion requires    '
     &,'           "constvmix", "ppvmix" or "tcvmix" to be enabled.    '
     &,'           If "constvmix" is enabled, "bryan_lewis_vertical"   '
     &,'           coefficients replace kappa_h. If "ppvmix" or        '
     &,'           "tcvmix" is enabled, "bryan_lewis_vertical"         '
     &,'           coefficients are added to the predicted coeffs.     '
        errorc = .true.
      endif
#endif
#if defined bryan_lewis_horizontal
# if !defined consthmix
        write (stdout,'(/,(1x,a/1x,a))')
     & '==> Error: "bryan_lewis_horizontal" tracer diffusion requires  '
     &,'           "consthmix" to be enabled.                          '
        errorc = .true.
# endif
#endif

      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No vertical mixing scheme has been enabled          '
        errorc = .true.
      endif
#if defined tcvmix
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tcvmix" cannot be enabled because another vertical '
     &,'           mixing scheme has been enabled                      '
        errorc = .true.
      else
        vmixset = .true.
      endif
#endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: there is no vertical mixing scheme enabled          '
        errorc = .true.
      endif
#if defined rigid_lid_surface_pressure
        write (stdout,9000)
     & '==> Warning: when using the "rigid_lid_surface_pressure" method'
     &,'             a checkerboard may develop in "ps". Using option  '
     &,'             "remove_ps_checkerboard" may help to remove it.   '
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
# if defined hypergrid || defined oldrelax
        write (stdout,9000)
     & '==> Error: "hypergrid" or "oldrelax" is not allowed when       '
     &,'           using "implicit_free_surface" or                    '
     &,'           "rigid_lid_surface_pressure" options.               '
        errorc = .true.
# endif
# if defined sf_5_point
        write (stdout,9000)
     & '==> Warning: "sf_5_point" is ignored when using the            '
     &,'             "implicit_free_surface" or                        '
     &,'             "rigid_lid_surface_pressure" option.              '
     &,'             conjugate_gradient with 9 point numerics is used. '
# endif
#endif

#if !defined stream_function && !defined implicit_free_surface
# if !defined rigid_lid_surface_pressure
        write (stdout,9000)
     & '==> Error: choose "stream_function", "implicit_free_surface"   '
     &,'           or the "rigid_lid_surface_pressure"  method         '
        errorc = .true.
# endif
#endif

#if defined stream_function
# if !defined sf_5_point && !defined sf_9_point
        write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using      '
     &,'           the stream_function method                          '
        errorc = .true.
# endif
# if defined sf_5_point && defined sf_9_point
        write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using      '
     &,'           the stream_function method (not both)               '
        errorc = .true.
# endif
#endif
#if defined coarse_grained_parallelism
        write (stdout,9000)
     & '==> Warning: "coarse_grained_parallelism" is enabled. If un-   '
     &,'             familiar with this option, use script "rum_mwsim" '
     &,'             to execute the memory window  simulator "mwsim.F" '
     &,'             to see how the latitude rows are tasked.          '
        write (stdout,*) '              Note: '
     &, num_processors, ' processors were requested.'
#else

!     look for potential load imbalance

      n = mod((jmt-2),ncrows)
      if (n .ne. 0) then
        num_mw = (jmt-2)/ncrows + (jmt-3)/(ncrows*((jmt-2)/ncrows))
        write (stdout
     &,'(1x,a,i4,/14x,a,i4,a,i4/14x,a,i4,/14x,a,i4,a,/14x,a,i4,a)')
     &'==> Warning: the number of rows (jmw) in the MW (memory window)='
     &, jmw, 'with calculations starting on MW row =', jsmw
     &,'and ending on MW row =', jemw
# if defined fourth_order_window
     &                               -1
# endif
     &,'disk latitude rows to be calculated are from 2 through', jmt-1
     &,'to compute these rows, the MW is moved northward '
     &, num_mw,' times.','however, the last time there are only ',n
     &,' row(s) to calculate.'
        write (stdout,'(/14x,a,/,14x,a,/,14x,a)')
     &  'If you are using multiprocessors, this will result in a load'
     &, 'imbalance on the last MW. The following is a list of "jmw"'
     &, 'which will not cause a load imbalance:'
# if defined fourth_order_window
        jb = 3
        jjs  = 4
# else
        jb = 2
        jjs  = 3
# endif
        do jj=jjs,jmt-2
          ncrow = jj-jb
          if (mod(jmt-2,ncrow) .eq. 0) then
            write (stdout,'(14x,a,i4,a,i4)')
     &      'jmw = ', jj, ',  number of MW moves = ',(jmt-2)/ncrow
          endif
        enddo
      endif
#endif
#if defined fourth_order_window
      if (jmw .lt. 4) then
        write (stdout,9000)
     & '==> Error: the MW can not have fewer than 4 rows when using any'
     &,'           fourth order options                                '
        write (stdout,*)'          you have set jmw=',jmw
        errorc = .true.
      endif
      if (jmw .gt. 4) then
        write (stdout,9000)
     & '==> Warning: "jmw" > 4 ("jmw"=4 will use the minimum memory)   '
        write (stdout,*)'          you have set jmw=',jmw
      endif
#else
      if (jmw .lt. 3) then
        write (stdout,9000)
     & '==> Error: the MW can not have fewer than 3 rows (jmw >2)      '
        write (stdout,*)'          you have set jmw=',jmw
        errorc = .true.
      endif
      if (jmw .gt. 3) then
        write (stdout,9000)
     & '==> Warning: "jmw" > 3 ("jmw"=3 will use the minimum memory)   '
      endif
#endif
      if (jmw .gt. jmt) then
        write (stdout,9000)
     & '==> Error: the MW can not have more rows than "jmt"            '
        write (stdout,*)'          you have set jmw=',jmw, ', jmt=',jmt
        errorc = .true.
      endif
      if (jmw .eq. jmt) then
        write (stdout,9000)
     & '==> Warning: The MW is open all the way ("jmw" = "jmt") which  '
     &,'             is the maximum memory configuration. Note that    '
     &,'             latitude rows are kept in the MW and not on disk! '
      endif

      if (nkflds .lt. 2) then
        write (stdout,9000)
     & '==> Error: "nkflds" must be at least 2                         '
        write (stdout,*)'          nkflds is set = ',nkflds
        errorc = .true.
      endif

#if defined equivalence_mw
        write (stdout,9000)
     & '==> Warning: The MW is equivalenced to the poisson solver      '
     &,'             workspace to reduce memory requirements. Verify   '
     &,'             that results are correct by comparing to a case   '
     &,'             with the "equivalence_mw" option removed.         '
      if (jmw .eq. jmt) then
        write (stdout,9000)
     & '==> Error: The MW cannot be equivalenced to the poisson solver '
     &,'           workspace but when "jmw" = "jmt"                    '
     &,'           remove the "equivalence_mw" option!                 '
        errorc = .true.
      endif
#endif

#if defined restorst
      if (dampts(1) .eq. c0 .or. dampts(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping time scale "dampts" must be set > 0     '
     &,'           when enabeling option "restorst"                    '
        errorc = .true.
      endif
      if (dampts(1) .ne. dampts(2)) then
        write (stdout,9000)
     & '==> Warning: are you aware that the damping time scale "dampts"'
     &,'             is set differently for T and S in option          '
     &,'             "restorst"? (they can be set the same in the      '
     &,'             namelist)                                         '
      endif
      if (dampdz(1) .eq. c0 .or. dampdz(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping thickness "dampdz" must be set nonzero  '
     &,'           when enabeling the "restorst" otpion                '
        errorc = .true.
      endif
      if (dampdz(1) .ne. dampdz(2)) then
        write (stdout,9000)
     & '==> Warning: are you aware that the damping thickness "dampdz" '
     &,'             is set differently for T and S in option          '
     &,'             "restorst"?                                       '
      endif
#else
      if (dampts(1) .ne. c0 .or. dampts(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping time scale "dampts" is > zero but     '
     &,'             the "restorst" otpion is not enabled              '
      endif
      if (dampdz(1) .ne. c0 .or. dampdz(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping thickness "dampdz" is > zero but      '
     &,'             the "restorst" otpion is not enabled              '
      endif
#endif

# if defined fourth_order_tracer_advection
        write (stdout,9000)
     & '==> Warning: "fourth_order_tracer_advection" has a more        '
     &,'             stringent stability requirement than the standard '
     &,'             second order advection. If problems arise, the    '
     &,'             time step may have to be reduced. Refer to the    '
     &,'             manual                                            '
# endif

#if defined isopycmix || defined isoneutralmix
# if defined consthmix && !defined biharmonic
        write (stdout,9000)
     & '==> Note: consthmix will only affect mixing of momentum        '
     &,'          since isopycmix was specified for tracer diffusion.  '
     &,'          kappa_h and Ah will be used as background mixing     '
     &,'          coefficients                                         '
# endif
# if defined biharmonic
        write (stdout,9000)
     & '==> Error: biharmonic is incompatable with isopycmix           '
     &,'           (use "smagnlmix" instead)                           '
        errorc = .true.
# endif
#else
# if defined gent_mcwilliams
        write (stdout,9000)
     & '==> Error: "isopycmix" must be enabled for "gent_mcwilliams"   '
     &,'           to work.                                            '
        errorc = .true.
# endif
#endif

#if defined obc
# if defined symmetry
        write (stdout,9000)
     &  '==> Error: "obc" and "symmetry" together is not compatible    '
     &, '           choose one or the other option.                    '
        errorc = .true.
# endif
# if !defined conjugate_gradient
        write (stdout,9000)
     &  '==> Error: "obc" works only with "conjugate gradient".        '
     &, '           It is not tested with other poisson solvers.       '
        errorc = .true.
# endif
# if !defined consthmix || !defined constvmix
        write (stdout,9000)
     &  '==> Warning: At the open boundaries laplacian formulation is  '
     &, '             used.                                            '
# endif
# if !defined stream_function
        write (stdout,9000)
     &  '==> Error: "obc" works only with "stream_function". It is     '
     &, '           not tested with other external mode solvers.       '
        errorc = .true.
# endif
# if defined cyclic
#  if defined obc_west || defined obc_east
        write (stdout,9000)
     &  '==> Error: "obc_west" or "obc_east" and "cyclic" together is  '
     &, '           not compatible. Choose one or the other option.    '
        errorc = .true.
#  endif
# endif
#endif

#if defined orlanski && !defined obc
        write (stdout,9000)
     &  '==> Error: For Orlanski radiation condition you have to       '
     &, '           choose an open boundary.                           '
        errorc = .true.
#endif

#if defined obcparameter
# if !defined obc
        write (stdout,9000)
     &  '==> Error: For writing obcparameter into snapshots you have   '
     &, '           to choose open boundaries.                         '
        errorc = .true.
# endif
# if !defined snapshots
        write (stdout,9000)
     &  '==> Error: For writing obcparameter into snapshots you have   '
     &, '           to choose "snapshots".                             '
        errorc = .true.
# endif
#endif

#if defined implicit_free_surface && defined stream_function
        write (stdout,9000)
     & '==> Error: "implicit_free_surface" and "stream_function" are   '
     &,'           incompatible. choose one or the other option.       '
        errorc = .true.
#endif
#if defined rigid_lid_surface_pressure && defined stream_function
        write (stdout,9000)
     & '==> Error: "rigid_lid_surface_pressure" and "stream_function"  '
     &,'           are incompatible. choose one or the other option.   '
        errorc = .true.
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
#if defined diagnostic_surf_height
        write (stdout,9000)
     & '==> Warning: "diagnostic_surf_height" does not apply under     '
     &,'             the "rigid_lid_surface_pressure" or               '
     &,'             "implicit_free_surface" options.                  '
#endif
#endif
#if defined stream_function
# if !defined conjugate_gradient && !defined oldrelax
#  if !defined hypergrid
        write (stdout,9000)
     & '==> Error: either "conjugate_gradient", "oldrelax", or         '
     &,'           "hypergrid" must be used with stream_function.      '
        errorc = .true.
#  endif
# endif
# if defined oldrelax && defined hypergrid
        write (stdout,9000)
     & '==> Error: "oldrelax" and "hypergrid" are incompatible         '
     &,'           choose one or the other option.                     '
        errorc = .true.
# endif
# if defined oldrelax && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error: "oldrelax" and "conjugate_gradient" are             '
     &,'           incompatible. choose one or the other option.       '
        errorc = .true.
# endif
# if defined oldrelax && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error: "oldrelax" and "conjugate_gradient" are             '
     &,'           incompatible. choose one or the other option.       '
        errorc = .true.
# endif
# if defined hypergrid && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error: "hypergrid" and "conjugate_gradient" are            '
     &,'           incompatible. choose one or the other option.       '
        errorc = .true.
# endif
# if defined hypergrid && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error: "hypergrid" and "conjugate_gradient" are            '
     &,'           incompatible. choose one or the other option.       '
        errorc = .true.
# endif
# if defined sf_9_point
        write (stdout,9000)
     & '==> Warning: although "sf_9_point" conserves total energy      '
     &,'             it has a 2 grid point null space.                 '
# endif
# if defined sf_5_point
        write (stdout,9000)
     & '==> Warning: although "sf_5_point" has no null space, it does  '
     &,'             not conserve total energy.                        '
# endif
# if defined sf_9_point && defined oldrelax
        write (stdout,9000)
     & '==> Error: "oldrelax" cannot be used with 9 point numerics     '
     &,'           use sf_5_point with oldrelax.                       '
        errorc = .true.
# endif
#endif

#if defined free_slip
        write (stdout,9000)
     & '==> Warning: "free-slip" only accounts for no transfer of      '
     &,'              tangential viscous stress to lateral boundaries. '
     &,'              Non-linear effects have been neglected for now,  '
     &,'              but a complete "free_slip" is being researched   '
     &,'              Note: do not use yet... still needs work         '
#endif

#if defined biharmonic && defined free_slip
        write (stdout,9000)
     & '==> Error: "biharmonic" does not account for "free_slip"       '
        errorc = .true.
#endif
#if defined biharmonic && !defined consthmix
        write (stdout,9000)
     & '==> Error: "consthmix" must be enabled when using "biharmonic" '
        errorc = .true.
#endif
#if defined isopycmix || defined isoneutralmix
# if defined consthmix
        if ((ah+ahisop) .gt. 1.e11) then
          write (stdout,9000)
     & '==> Error: "ahisop"+"ah" is too large for the                  '
     &,'           "isopycmix" mixing option                           '
          errorc = .true.
        endif
# endif
#endif

        if (dtsf .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the external mode time step "dtsf"      '
          errorc = .true.
        endif

        if (dtuv .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the internal mode time step "dtuv"      '
          errorc = .true.
        endif

        if (dtts .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the density time step "dtts"            '
          errorc = .true.
        endif

        critv = 1.e-6
        if (mod(rundays,dtts*secday) .gt. critv) then
          t1 = nint(rundays/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: there must be an integral number of density time  '
     &,'             steps within "rundays" (the integration time).    '
          write (stdout,*) '               (changed "rundays" from     '
     &,   rundays,' days to ', t1,' days to insure this condition)     '
          rundays = t1
        endif

        if (itmb) then
          write (stdout,9000)
     & '==> Warning: "itmb" is set to "true". set it to "false" in     '
     &,'             subsequent runs to prevent the time independent   '
     &,'             basin mask from being written more than once.     '
     &,'             This reduces the size of the diagnostic file.     '
        endif

        if (itrmb) then
          write (stdout,9000)
     & '==> Warning: "itrmb" is set to "true". set it to "false" in    '
     &,'             subsequent runs to prevent the time independent   '
     &,'             region masks from being written more than once.   '
     &,'             This reduces the size of the diagnostic file.     '
        endif

        if (itavg) then
          write (stdout,9000)
     & '==> Warning: "itavg" is set to "true". set it to "false" in    '
     &,'             subsequent runs to prevent the time independent   '
     &,'             region masks from being written more than once.   '
     &,'             This reduces the size of the diagnostic file.     '
        endif
#if defined trajectories
        if (trajint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: "trajectories" is set but "trajint"=0. this means '
     &,'             particle positions will be written every time step'
        endif
        if (trajint .lt. c0) then
          write (stdout,9000)
     & '==> Error: "trajectories" is set but "trajint" < 0. this means '
     &,'            particle positions will not be written.            '
          errorc = .true.
        endif
# if !defined lyapunov
          write (stdout,9000)
     & '==> Warning: "trajectories" is enabled but "lyapunov" is not.  '
# endif
#else
# if defined lyapunov
          write (stdout,9000)
     & '==> Warning: "trajectories" is not enabled but "lyapunov" is.  '
# endif
#endif
#if defined meridional_tracer_budget
        if (mod(tmbint,dtts*secday) .gt. critv .and. tmbint .gt. c0)
     &    then
          t1 = nint(tmbint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "tmbint" does not contain an integral number of   '
     &,'             density time steps "dtts".                        '
          write (stdout,*) '               (changed "tmbint" from '
     &,   tmbint,' days to ', t1,' days to insure this condition)'
          tmbint = t1
        endif
        if (tmbint .eq. c0) then
          write (stdout,9000)
     & '==> Error: averaging interval "tmbint"  must be specified when '
     &,'           the "meridional_tracer_budget" option is enabled    '
          errorc = .true.
        endif
        if (tmbint .gt. tmbper) then
          write (stdout,9000)
     & '==> Warning: the interval "tmbint" exceeds the averaging       '
     &,'             period "tmbper" for option                        '
     &,'             "meridional_tracer_budget"                        '
        endif
        if (tmbint .lt. tmbper) then
          write (stdout,9000)
     & '==> Warning: averaging period "tmbper" exceeds interval        '
     &,'             "tmbint". Setting tmbper = tmbint                 '
          tmbper = tmbint
        endif
        if (tmbper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "tmbper" is zero. The        '
     &,'             average will be over only one time step!          '
        endif
#else
      if (tmbint .gt. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "tmbint" is > zero but the '
     &,'             the "meridional_tracer_budget" option is not on.  '
      endif
#endif
#if defined time_averages
        if (mod(timavgint,dtts*secday).gt.critv .and. timavgint .gt. c0)
     &    then
          t1 = nint(timavgint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "timavgint" does not contain an integral number of'
     &,'             density time steps "dtts".                        '
          write (stdout,*) '               (changed "timavgint" from '
     &,   timavgint,' days to ', t1,' days to insure this condition)'
          timavgint = t1
        endif
        if (timavgint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: averaging interval "timavgint" = 0.0 implies no   '
     &,'             averaging when option "time_averages" is enabled  '
        endif
        if (timavgint .gt. timavgper) then
          write (stdout,9000)
     & '==> Warning: the interval "timavgint" exceeds the averaging    '
     &,'             period "timavgper" for option "time_averages"     '
        endif
        if (timavgint .lt. timavgper) then
          write (stdout,9000)
     & '==> Warning: averaging period "timavgper" exceeds interval     '
     &,'             "timavgint". Setting timavgper = timavgint        '
          timavgper = timavgint
        endif
        if (timavgper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "timavgper" is zero. The     '
     &,'             average will be over only one time step!          '
        endif
#else
      if (timavgint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "timavgint"  is > zero but   '
     &,'             the "time_averages" option is not enabled         '
      endif
#endif
#if defined xbts
        if (mod(xbtint,dtts*secday) .gt. critv .and. xbtint .gt. c0)
     &    then
          t1 = nint(xbtint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "xbtint" does not contain an integral number of   '
     &,'    density time steps "dtts".                                 '
          write (stdout,*) '               (changed "xbtint" from '
     &,   xbtint,' days to ', t1,' days to insure this condition)'
          xbtint = t1
        endif
        if (xbtint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: averaging interval "xbtint" = 0.0 implies no      '
     &,'             averaging when option "xbts" is enabled           '
        endif
        if (numxbt .eq. 0) then
          write (stdout,9000)
     &   '==> Error: the number of XBT stations detected is zero       '
     &,  '           at least one station must be specified            '
          errorc = .true.
        endif
        if (xbtint .gt. xbtper) then
          write (stdout,9000)
     & '==> Warning: the interval "xbtint" exceeds the  averaging      '
     &,'             period "xbtper" for option "xbts"                 '
        endif
        if (xbtint .lt. xbtper) then
          write (stdout,9000)
     & '==> Warning: averaging period "xbtper" exceeds interval        '
     &,'             "xbtint". Setting xbtper = xbtint                 '
          xbtper = xbtint
        endif
        if (xbtper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "xbtper" is zero. The        '
     &,'             average will be over only one time step!          '
        endif
#else
      if (xbtint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "xbtint"  is > zero but    '
     &,'             the "xbts" option is not enabled                  '
      endif
#endif
#if defined diagnostic_surf_height
        if (mod(dspint,dtts*secday) .gt. critv .and. dspint .gt. c0)
     &    then
          t1 = nint(dspint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "dspint" does not contain an integral number of   '
     &,'             density time steps "dtts".                        '
          write (stdout,*) '               (changed "dspint" from '
     &,   dspint,' days to ', t1,' days to insure this condition)'
          dspint = t1
        endif
        if (dspint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: averaging interval "dspint" = 0.0 implies no      '
     &,'             averaging when option "diagnostic_surf_height" is '
     &,'             enabled                                           '
        endif
        if (dspint .gt. dspper) then
          write (stdout,9000)
     & '==> Warning: the interval "dspint" exceeds the averaging       '
     &,'             period "dspper" for option                        '
     &,'             "diagnostic_surf_height"                          '
        endif
        if (dspint .lt. dspper) then
          write (stdout,9000)
     & '==> Warning: averaging period "dspper" exceeds interval        '
     &,'             "dspint". Setting dspper = dspint                 '
          dspper = dspint
        endif
        if (dspper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "dspper" is zero. The        '
     &,'             average will be over only one time step!          '
        endif
#else
      if (dspint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "dspint"  is > zero but    '
     &,'             option "diagnostic_surf_height" is not enabled    '
      endif
#endif

      if ((dtuv .ne. dtsf) .or. (dtuv .ne. dtts)) then
        write (stdout,9000)
     & '==> Warning: use of unequal time steps implies the transient   '
     &,'             response is unimportant and multiple equilibria   '
     &,'             do not exist.                                     '
      endif

!     check for mixing coefficients larger than stability permits

      dymin  = dyt(2)
      dxmin  = dxt(2)
      do jrow=2,jmtm1
        dymin  = min(dymin,dyt(jrow))
      enddo
      do i=2,imtm1
        dxmin  = min(dxmin,dxt(i))
      enddo
      dzmin  = dzt(1)
      xlmax  = dtxcel(1)
      do k=2,km
        xlmax  = max(xlmax,dtxcel(k))
        dzmin  = min(dzmin,dzt(k))
      enddo

      if (xlmax .gt. c1) then
        write (stdout,9000)
     & '==> Warning: use of accelerated time steps implies the         '
     &,'             transient response is unimportant and multiple    '
     &,'             equilibria do not exist. stability tests will     '
     &,'             use "dtts" multiplied by the maximum "dtxcel"     '
      endif

#if defined consthmix
      dtxl = dtts*xlmax
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
# if defined isopycmix || defined isoneutralmix
        fimax = 0.
        do k=1,km
          do i=2,imtm1
            fimax = max(fimax,fisop(i,j,k))
          enddo
        enddo
        if ((dtxl*(ah+ahisop*fimax))/dxdymn .ge. p25) then
# else
        if ((dtxl*ah)/dxdymn .ge. p25) then
# endif
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
# if defined isopycmix || defined isoneutralmix
     & '==> Warning: lateral diffusive criteria exceeded for "ah" +    '
     &,'             "ahisop". use a smaller "dtts", "dtxcel", and/or  '
     &,'             "ah" + "ahisop"                                   '
# else
     & '==> Warning: lateral diffusive criteria exceeded for "ah".     '
     &,'              use a smaller "dtts", "dtxcel", and/or "ah"      '
# endif
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
        if ((dtuv*am)/dxdymn .ge. p25) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: lateral diffusive criteria exceeded for "am".     '
     &,'             use a smaller "dtuv" and/or "am"                  '
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
#endif
#if defined sponges
# if !defined source_term
      write (stdout,9000)
     & '==> Error: ifdef "sponges" is enabled but "source_term"        '
     &,'           is not. "source_term" needs to be enabled           '
      errorc = .true.
# endif
#endif
#if defined shortwave
# if !defined source_term
      write (stdout,9000)
     & '==> Error: ifdef "shortwave" is enabled but "source_term"      '
     &,'           is not. "source_term" needs to be enabled           '
      errorc = .true.
# endif
# if !defined uvic_embm
      write (stdout,9000)
     & '==> Warning: ifdef "shortwave" is enabled but "coupled" is not.'
     &,'             Is the shortwave component of radiation defined ? '
# endif
#else
      if (dzt(1) .lt. 20.0e2) then
        write (stdout,9000)
     & '==> Warning: if shallow mixed layers develop, then enabling    '
     &,'             ifdef "shortwave" may help to deepen them. note   '
     &,'             that either you or the atmosphere must provide    '
     &,'             the solar short wave as a boundary condition.     '
      endif
#endif
#if defined constvmix
      do k=1,km
        if ((dtts*dtxcel(k)*kappa_h)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_h" '
     &,'             use a smaller "dtts", "dtxcel", and/or "kappa_h"  '
         write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*kappa_m)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_m" '
     &,'             use a smaller "dtuv" and/or "kappa_m"             '
      endif
#endif
#if defined implicitvmix
      write (stdout,9000)
     & '==> Warning: "implicitvmix" turns convective adjustment off.   '
     &,'             It will be handled by the vertical mixing         '
# if defined fullconvect
      write (stdout,9000)
     & '==> Error: the full convective scheme is incompatible with     '
     &,'           option "implicitvmix"                               '
        errorc = .true.
# endif
#else
# if defined fullconvect
      write (stdout,9000)
     & '==> Warning: the full convective scheme is enabled.            '
     &,'             it will ignore "ncon" and remove all instability  '
# else
      write (stdout,9000)
     & '==> Warning: the standard convective scheme is enabled. Verify '
     &,'             that "ncon" is set to what is wanted:             '
      write (stdout,'(32x,a,i3)') '"ncon" is currently set to ',ncon
# endif
#endif
#if defined tcvmix && !defined implicitvmix
      do k=1,km
        if ((dtts*dtxcel(k)*diff_cbt_limit)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
!-----------------------------------------------------------------------
     & '==> Warning: vertical diffusive criteria exceeded for          '
     &,'             "diff_cbt_limit".  use a smaller "dtts", "dtxcel" '
     &,'             ,and/or "diff_cbt_limit"                          '
          write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*visc_cbu_limit)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded for          '
     &,'             "visc_cbu_limit". use a smaller "dtuv" or         '
     &,'             "visc_cbu_limit"                                  '
      endif
#endif

!     check range of implicit factors

#if !defined implicitvmix && !defined isopycmix && !defined isoneutralmix
      if (aidif .ne. c0) then
        write (stdout,9000)
     & '==> Warning: "aidif" was not = 0.0  but "implicivmix" or       '
     &,'             "isopycmix" is not defined, so vertical mixing    '
     &,'              will be fully explicit.                          '
      endif
#endif
#if defined damp_inertial_oscillation
      if (acor .lt. p5) then
        write (stdout,9000)
     & '==> Error: "acor" must be between 0.5 & 1.0 when option        '
     &,'           damp_inertial_oscillation is enabled.               '
          errorc = .true.
      endif
# if defined fct
      write (stdout,9000)
     & '==> Error: non-linear instability will result if "fct" is used '
     &,'           "damp_inertial_oscillation"                         '
          errorc = .true.
# endif
#else
      if (acor .ne. 0) then
        write (stdout,9000)
     & '==> Error: "acor" must=0 when option damp_inertial_oscillation '
     &,'           is not enabled.                                     '
          errorc = .true.
      else

!       check for marginally resolved inertial oscillation

        fmax = epsln
        do jrow=2,jmtm1
          do i=2,imtm1
            fmax = max(fmax,abs(cori(i,jrow,1)))
          enddo
        enddo
        if (dtuv .gt. (1.0/6.0)*(c2*pi)/fmax) then
          write (stdout,9000)
     & '==> Error: the inertial oscillation is not resolved. reduce    '
     &,'           "dtuv" or use option "damp_inertial_oscillation"    '
          errorc = .true.
        endif
      endif
#endif
#if defined consthmix

!-----------------------------------------------------------------------
!     search for topographic instabilities (based  on the  work of
!     Peter Killworth  ...  eqn 11 from ocean modeling nov 1987)
!-----------------------------------------------------------------------

      num   = 50
      do j=2,jmtm1
        dysq = dyt(j)**2
        do i=2,imtm1
          if (kmu(i+1,j-1) .ne. 0 .and. kmu(i+1,j) .ne. 0) then
            dxsq = (dxt(i)*cst(j))**2
            clix = am*dtuv/dxsq
            h1   = zw(kmu(i+1,j-1))
            h2   = zw(kmu(i+1,j))
            hx   = (8.0*h1*h2/(h1+h2)**2 + dxsq/dysq)/(4.0 + dxsq/dysq)
            if (clix .ge. hx .and. num .ge. 0) then
              num = num - 1
              write(stdout,*)
              write (stdout,'(a,a,i4,a,i4,a)')
     &        '==> Warning: Killworth topographic roughness condition'
     &,       ' exceeded at location (i,j) = (',i+1,',',j,')'
              if (num .eq. 0) then
                write (stdout,9000)
     &         '==> Warning: msgs terminated after 50 cases were found '
              endif
            endif
          endif
        enddo
      enddo
#endif

!     verify that the domain boundary is valid

#if defined cyclic && defined solid_walls
        write (stdout,9000)
     & '==> Error: Both "cyclic" and "solid_walls cannot be enabled    '
        errorc = .true.
#endif
#if !defined cyclic && !defined solid_walls && !defined obc
        write (stdout,9000)
     & '==> Error: Either "cyclic", "solid_walls
     &  or "obc" must be enabled'
        errorc = .true.
#endif
#if defined solid_walls
# if !defined obc_west && !defined obc_east
      jw = 0
      je = 0
      do j=1,jmt
        if (kmt(1,j) .ne. 0) jw = j
        if (kmt(imt,j) .ne. 0) je = j
      enddo
      if (je .ne. 0 .or. jw .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the eastern or western boundary. This is not        '
     &,'           correct when enabeling option "solid_walls"         '
        if (jw .ne. 0) write (stdout,*) '           at i=1 and j=',jw
        if (je .ne. 0) write (stdout,*) '           at i=imt and j=',je
      endif
# endif
#endif
#if !defined obc_south && !defined obc_north
      in = 0
      is = 0
      do i=1,imt
        if (kmt(i,1) .ne. 0) is = i
        if (kmt(i,jmt) .ne. 0) in = i
      enddo
      if (is .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the southern boundary.                              '
        write (stdout,*) '           at j=1 and i=',is
      endif
#if !defined symmetry
      if (in .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the northern boundary.                              '
        write (stdout,*) '           at j=jmt and i=',in
      endif
#endif
#endif

!     verify that each ocean point is at least 2 levels deep

      numk = 0
      do jrow=1,jmt
        do i=1,imt
          if (kmt(i,jrow) .eq. 1) then
            numk = numk + 1
            errorc = .true.
            write (stdout,*)
     &       ' Error: kmt(',i,',',jrow,') = 1 is not allowed    '
          endif
        enddo
      enddo
      if (numk .ne. 0) then
        write (stdout,9000)
     & '==> Error: "kmt" must be at least 2 levels deep at all ocean   '
     &,'           points.                                             '
      endif

#if !defined quicker && defined ncar_upwind3
      write (stdout,9000)
     & '==> Error: "ncar_upwind3" only works when option               '
     &,'           option "quicker"  is enabled                        '
      errorc = .true.
#endif

#if defined  fct
# if  defined fourth_order_tracer_advection
        write (stdout,9000)
     & '==> Error: "fct" is not compatable with option                 '
     &,'           "fourth_order_tracer_advection"                     '
        errorc = .true.
# endif
# if  defined quicker
        write (stdout,9000)
     & '==> Error: "fct" is not compatable with option                 '
     &,'           "quicker"                                           '
        errorc = .true.
# endif
# if defined fct_dlm1 && defined fct_dlm2
      write (stdout,9000)
     & '==> Error: more than one fct delimiter type is specified       '
     &,'           specify either fct_dlm1 or fct_dlm2                 '
     &,'           Default: using fct_dlm1                             '
        errorc = .true.
# endif
# if !defined fct_dlm1 && !defined fct_dlm2
      write (stdout,9000)
     & '==> Warning: fct delimiter type is not specified               '
     &,'             specify either fct_dlm1 or fct_dlm2               '
     &,'    Default: using fct_dlm1                                    '
# endif
#else
# if defined fct_dlm1 || defined fct_dlm2 || defined fct_3d
      write (stdout,9000)
     & '==> Warning: main fct option is not specified                  '
     &,'             minor fct options will be ignored                 '
# endif
#endif
#if defined pressure_gradient_average && !defined fourth_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option '
     &,'           "pressure_gradient_average". in "size.h" but it is  '
     &,'            not.                                               '
        errorc = .true.
#endif
#if defined biharmonic && !defined fourth_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option '
     &,'           "biharmonic". in "size.h" but it is not.            '
        errorc = .true.
#endif
#if defined fct && !defined fourth_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option '
     &,'           "fct". in "size.h" but it is not.                   '
        errorc = .true.
#endif
#if defined quicker && !defined fourth_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option '
     &,'           "quicker". in "size.h" but it is not.               '
        errorc = .true.
#endif
#if defined fourth_order_tracer_advection && !defined fourth_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option '
     &,'           "fourth_order_tracer_advection". in "size.h" but it '
     &,'            is not.                                            '
        errorc = .true.
#endif
#if defined held_larichev
      write (stdout,9000)
     & '==> Error: "held_larichev" is still under development and is   '
     &,'           unsuitable for use until it is completed.           '
        errorc = .true.
#endif
      write (stdout,'(/,20x,a,/)')
     &         ' E N D    C O N S I S T E N C Y    C H E C K S'
      if (errorc) stop '=>checks'

      return
9000  format (/,(1x,a))
      end
