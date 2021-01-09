#if defined drive_mwsim
      program mwsim

!=======================================================================
!     Memory Window Simulator shows how the memory window handles
!     various combinations of "jmw", "jmt", and options

!     based on code by: R. C. Pacanowski
!=======================================================================

      parameter (jmt=8)

!------------------------ parameters from size.h ----------------------

#if defined pressure_gradient_average
# define fourth_order_window
#endif

#if defined biharmonic
# define fourth_order_window
#endif

#if defined fourth_order_tracer_advection
# define fourth_order_window
#endif

#if defined fct
# define fourth_order_window
#endif

#if defined coarse_grained_parallelism

!     for coarse grained parallelism: use the following setting for jmw
!       (don`t monkey around with this unless you know what to do)

# if defined fourth_order_window
#  if defined pressure_gradient_average
#   if defined biharmonic  || defined fourth_order_tracer_advection || fct
      parameter (jmw=7)
#   else
      parameter (jmw=5)
#   endif
#  else
      parameter (jmw=5)
#  endif
# else
      parameter (jmw=3)
# endif
#else

!     for uni-tasking: "jmw" is set to minimum for each option class
!     "jmw" may be increased up to "jmt" for fine grained parallelism

# if defined fourth_order_window
#  if defined pressure_gradient_average
#   if defined biharmonic  || defined fourth_order_tracer_advection || fct
      parameter (jmw=5)
#   else
      parameter (jmw=4)
#   endif
#  else
      parameter (jmw=4)
#  endif
# else
      parameter (jmw=3)
# endif
#endif
      parameter (jsmw=2, jemw=jmw-1)
#if defined coarse_grained_parallelism
      parameter (jextra = 0
# if defined fourth_order_window
#  if defined pressure_gradient_average
     &                    + 1
#  endif
#  if defined biharmonic || defined fourth_order_tracer_advection || fct
     &                    + 1
#  endif
# endif
     &                    )
#else
          parameter (jextra = 0)
#endif

!---------------------  parameters from param.h  ----------------------

!     ncrows = number of calculated rows starting with row 2 in the MW.
!              (the remaining rows are buffer rows).

#if defined coarse_grained_parallelism
      parameter (ncrows = 1)
#else
# if defined fourth_order_window
#  if defined pressure_gradient_average
#   if defined biharmonic || defined fourth_order_tracer_advection || fct
      parameter (ncrows = jmw - 4 + 2*(jmw/jmt))
#   else
      parameter (ncrows = jmw - 3 + jmw/jmt)
#   endif
#  else
      parameter (ncrows = jmw - 3 + jmw/jmt)
#  endif
# else
      parameter (ncrows = jmw - 2)
# endif
#endif

!------------------------- start of  mwsim --------------------------

      logical first_mw
        last_clinic_row(joff) = min(jsmw+ncrows-1+jextra,jmt-1-joff)
#if defined pressure_gradient_average
        last_tracer_row(joff) = min(jsmw+ncrows+jextra,jmt-1-joff)
#endif

      print '(/20x,a,i4,/)', 'Memory Window Simulator for jmt = ',jmt

!-----------------------------------------------------------------------
!       Solve equations for rows within each MW
!-----------------------------------------------------------------------

#if defined coarse_grained_parallelism

!     generate results for 1..jmt-2 processors

      do num_processors=1,jmt-2

!       ntasks = number of tasks needed to solve rows 2 through jmt-1.
!       set one task per processor

        ntasks = num_processors

        call checks (jmw, ntasks, ncrows, jmt)

        do n = 1,ntasks

!         define starting and ending rows for each task

          trows = float(jmt-2)/ntasks
          jstask = max(int((n-1)*trows - jextra + 1.0001),  1 - jextra)
          jetask = min(int(n*trows + 1 + jextra + 1.0001),jmt + jextra)

          print '(3x,a,i2,a,i2,a,i2)'
     &,  'TASK # ',n,',   jstask=',jstask,', jetask=',jetask

!         num_mwpt = number of memory windows per task

          num_mwpt = jetask - jstask + 1 - (jmw - ncrows)

          do mw = 1,num_mwpt
            first_mw = (mw .eq. 1)
            if (first_mw) then
              js   = 1
            else
              js = jmw - ncrows + 1
            endif
            joff = jstask + mw - 2
            if (js+joff .lt. 1) js = 1 - joff
            je = min(jmw,jmt-joff)
#else

!       num_mw  = number of memory windows needed to solve latitude
!                 rows 2 through jmt-1

        num_mw = (jmt-2)/ncrows + (jmt-3)/(ncrows*((jmt-2)/ncrows))
        do mw = 1,num_mw
          first_mw = (mw .eq. 1)
          if (first_mw) then
            js = 1
          else
            js = jmw - ncrows + 1
          endif

          joff = (mw-1)*ncrows
          je = min(jmw,jmt-joff)

# if defined fourth_order_window
          print '(3x,a,i2,a,i2,a,i2)'
     &,  'TESTING 4th order MW # ',mw,', with jmw = ',jmw
     &,  ', and jmt = ',jmt
# else
          print '(3x,a,i2,a,i2,a,i2)'
     &,  'TESTING 2nd order MW # ',mw,', with jmw = ',jmw
     &,  ', and jmt = ',jmt
# endif
#endif

!-----------------------------------------------------------------------
!           load latitude rows into memory window
!-----------------------------------------------------------------------

            call loadmw (joff, js, je, first_mw, jmw, ncrows)

!-----------------------------------------------------------------------
!           integrate
!-----------------------------------------------------------------------

            jscalc = 2 + jextra
            jecalc = last_clinic_row(joff)
#if defined pressure_gradient_average

!           compute one extra row of tracers

            jstrac = jscalc + 1 - 1/mw
            jetrac = last_tracer_row(joff)
#else
            jstrac = jscalc
            jetrac = jecalc
#endif
            call tracer (joff, jstrac, jetrac)
            call clinic (joff, jscalc, jecalc)

!-----------------------------------------------------------------------
!           save updated rows
!-----------------------------------------------------------------------

            call putmw  (joff, jscalc, jecalc)

          enddo
#if defined coarse_grained_parallelism
        enddo
      enddo
#endif
      print '(/a/)', 'Parameters can be changed in "mwsim.F"'
      stop
      end

      subroutine tracer (joff, jstrac, jetrac)
      if (jstrac .le. jetrac) then
        print '(5x,4(a,i2))'
     &,'Tracer: computing rows ',jstrac,' to ', jetrac
     &,', jrows = ',jstrac+joff, ' to ', jetrac+joff
      else
        print '(5x,4(a,i2))'
     &,'Tracer: no computation for rows',jstrac,' to ', jetrac
     &,', jrows = ',jstrac+joff, ' to ', jetrac+joff
      endif
      return
      end

      subroutine clinic (joff, jscalc, jecalc)
      if (jscalc .le. jecalc) then
        print '(5x,4(a,i2))'
     &,'Clinic: computing rows ',jscalc,' to ', jecalc
     &,', jrows = ',jscalc+joff, ' to ', jecalc+joff
      else
        print '(5x,4(a,i2))'
     &,'Clinic: no computation for rows ',jscalc,' to ', jecalc
     &,', jrows = ',jscalc+joff, ' to ', jecalc+joff
      endif
      return
      end

      subroutine putmw (joff, jscalc, jecalc)
      print '(5x,4(a,i2))'
     &,'Writing rows ',jscalc,' to ', jecalc
     &,', jrows = ',jscalc+joff, ' to ', jecalc+joff
      return
      end

      subroutine loadmw (joff, js, je, first_mw, jmw, ncrows)
      logical first_mw

      Print '(5x,a)','Loadmw...'
      if (.not. first_mw) then
        nrows = jmw - ncrows
        do move=1,nrows
          jfrom = jmw - (nrows - move)
          jto   = move
          print '(7x,a,i2,a,i2)'
     &,  'Copying data from row ',jfrom,' to ', jto
        enddo
      else
        print '(7x,a)', 'No copying of data'
      endif

      if (js .gt. je) then
        print '(7x,a,i2,a,i2,a,i2,a,i2,a,i2,a,l2)'
     &,   'Skip reading rows js=',js,', je=',je, ' joff=',joff
     &,   ' jrows =',js+joff,':',je+joff, ' first_mw= ',first_mw
      else
        print '(7x,a,i2,a,i2,a,i2,a,i2,a,i2,a,l2)'
     &,   'Reading rows js=',js,', je=',je, ' joff=',joff
     &,   ' jrows =',js+joff,':',je+joff, ' first_mw= ',first_mw
      endif
      return
      end

      subroutine checks (jmw, ntasks, ncrows, jmt)
      integer stdout
      parameter (stdout = 6)
      logical error
      error = .false.
      print *,' '
      print *,' '
      if (jmw .eq. 3) print '(a,i2,a,i2,a,i2)'
     &,'TESTING => 2nd order window with ',ntasks,' tasks. jmw =',jmw
     &,', jmt=',jmt
      if (jmw .eq. 5) print '(a,i2,a,i2,a,i2)'
     &,'TESTING => 4th order window with ',ntasks,' tasks. jmw =',jmw
     &,', jmt=',jmt
      if (jmw .eq. 7) print '(a,i2,a,i2,a,i2)'
     &,'TESTING => 4th order window + press_grad_aver with ',ntasks
     &,' tasks. jmw =',jmw
     &,', jmt=',jmt
#if defined coarse_grained_parallelism
      if (ntasks .gt. jmt-2) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  parameter "ntasks" cannot be > jmt-2'
        print *,'=>Error: ntasks=',ntasks,' > jmt-2'
        error = .true.
      endif
# if defined fourth_order_window && !defined pressure_gradient_average
      if (jmw .ne. 5) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  jmw must = 5 when "coarse_grained_parallelism"     '
     &,'            is enabled but not "pressure_gradient_average"     '
        error = .true.
      endif
# endif
# define other =  defined biharmonic || defined fourth_order_tracer_advection || fct
# if defined other_fourth_order && defined pressure_gradient_average
      if (jmw .ne. 7) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  jmw must = 7 when "coarse_grained_parallelism"     '
     &,'            is enabled but not "pressure_gradient_average"     '
        error = .true.
      endif
# endif
# if !defined fourth_order_window
      if (jmw .ne. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  jmw must = 3 when "coarse_grained_parallelism"     '
     &,'            is enabled.                                        '
        error = .true.
      endif
# endif
      if (ncrows .ne. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  "ncrows" must = 3 when "coarse_grained_parallelism"'
     &,'            is enabled.                                        '
        error = .true.
      endif
#endif

#if defined coarse_grained_parallelism
      print *,' ENABLED OPTION = coarse_grained_parallelism'
#endif
#if defined biharmonic
      print *,' ENABLED OPTION = biharmonic'
#endif
#if defined pressure_gradient_average
      print *,' ENABLED OPTION = pressure_gradient_average'
#endif

      return
      end
#else
      subroutine mwsim
      return
      end
#endif
