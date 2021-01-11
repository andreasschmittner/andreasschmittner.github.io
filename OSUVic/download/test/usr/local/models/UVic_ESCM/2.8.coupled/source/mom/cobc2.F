#if defined obc_west || defined obc_east
      subroutine cobc2 (joff, js, je, i)

!=======================================================================
!     compute tracers at "tau+1" for i= 2 and/or imt-1
!     using open boundary conditions similar to Stevens (1990).

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       i    = boundary longitude (2 or imt-1)

!     similar to "cobc.F" but at meridional boundaries

!     based on code by: A. Biastoch
!=======================================================================

# include "param.h"
# include "accel.h"
# include "cobc.h"
# include "coord.h"
# include "csbc.h"
# include "emode.h"
# include "grdvar.h"
# include "hmixc.h"
# if defined isopycmix
#  include "isopyc.h"
# endif
# include "levind.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt(km)
      dimension ADV_Tx_obc(km,jmw)
# include "fdift.h"

!     at open boundaries always laplacian formulation is used (at tau)

      DIFF_Ty_obc(i,k,j) = ahc_north(jrow)*tmask(i,k,j+1)*
     &                     (t(i,k,j+1,n,tau) - t(i,k,j,n,tau))
     &                   - ahc_south(jrow)*tmask(i,k,j-1)*
     &                     (t(i,k,j,n,tau) - t(i,k,j-1,n,tau))

!-----------------------------------------------------------------------

# if defined obc_west

      if (i .eq. 2) then

!-----------------------------------------------------------------------
!       1) compute the advective velocity "uad"
!          at the east face of the "t" grid box

!       2) compute phase velocity at the western boundary: c1w
!-----------------------------------------------------------------------

        var = -dxu(i+1)/dtts

        do k=1,km
          do j=js,je
            jrow = j + joff

            uad(k,j) = (u(i,k,j,1,tau)*dyt(j)+u(i,k,j-1,1,tau)
     &                 *dyt(j-1))/(dyt(j)+dyt(j-1))
            if (uad(k,j) .gt. c0) uad(k,j) = c0

            do m=1,nt
              var2(k,j,m) = t(i+2,k,j,m,taum1)-t(i+1,k,j,m,taum1)
              if (var2(k,j,m) .eq. c0) then
                c1w(k,j,m) = var
              else
                c1w(k,j,m) = var*(t(i+1,k,j,m,tau)-t(i+1,k,j,m,taum1))
     &                          /var2(k,j,m)
                if (c1w(k,j,m) .gt. c0)  c1w(k,j,m) = c0
                if (c1w(k,j,m) .lt. var) c1w(k,j,m) = var
              endif
            enddo

            if (jrow .eq. 2)     uad(k,j)=c0
            if (jrow .eq. jmtm1) uad(k,j)=c0
          enddo
        enddo
      endif
# endif

# if defined obc_east

      if (i .eq. imtm1) then

!-----------------------------------------------------------------------
!       1) compute the advective velocity "uad"
!          at the west face of the "t" grid box

!       2) compute phase velocity at the western boundary: c1e
!-----------------------------------------------------------------------

        var = dxu(i-2)/dtts

        do k=1,km
          do j=js,je
            jrow = j + joff

            uad(k,j) = (u(i-1,k,j,1,tau)*dyt(j)+u(i-1,k,j-1,1,tau)
     &                 *dyt(j-1))/(dyt(j)+dyt(j-1))
            if (uad(k,j) .lt. c0) uad(k,j) = c0

            do m=1,nt
              var2(k,j,m) = t(i-1,k,j,m,taum1)-t(i-2,k,j,m,taum1)
              if (var2(k,j,m) .eq. c0) then
                c1e(k,j,m) = var
              else
                c1e(k,j,m) = -var*(t(i-1,k,j,m,tau)
     &                             -t(i-1,k,j,m,taum1))/var2(k,j,m)
                if (c1e(k,j,m) .lt. c0)  c1e(k,j,m) = c0
                if (c1e(k,j,m) .gt. var) c1e(k,j,m) = var
              endif
            enddo

            if (jrow .eq. 2)     uad(k,j) = c0
            if (jrow .eq. jmtm1) uad(k,j) = c0
          enddo
        enddo
      endif
# endif

      do n=1,nt

!-----------------------------------------------------------------------
!       diffusive flux across northern face of "T" cells
!       (use constant horizontal diffusion)
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km
            do ii=i-1,i

!             only del**2 diffusion

              jrow = j + joff
              ah_cstdxur(ii,j) = ah*cstr(jrow)*dxur(ii)
              diff_fe(ii,k,j)  = ah_cstdxur(ii,j)*
     &                           (t(ii+1,k,j,n,tau) - t(ii,k,j,n,tau))
              cstdxtr(ii,j)    = cstr(jrow)*dxtr(ii)
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       diffusive flux across bottom face of "T" cells
!       use constant vertical diffusion
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,kmm1
            diff_fb(i,k,j) = diff_cbt(i,k,j)*dzwr(k)*
     &                          (t(i,k,j,n,tau) - t(i,k+1,j,n,tau))
          enddo
        enddo

!       for vertical diffusion of tracers:
!       set k=0 elements to reflect surface tracer flux.
!       set bottom level "kz" to reflect insulation condition.
!       for vertical advection of tracers:
!       set k=0 elements of vertical advection array to the rigid
!       lid condition (w(0) = 0). also assume w(km)=0. If "kz" <
!       "km" and there is a bottom slope, w(kz) may be non zero.

        do j=js,je
          jrow            = j + joff
          kb              = kmt(i,jrow)
          diff_fb(i,0,j)  = stf(i,j,n)
          diff_fb(i,kb,j) = btf(i,j,n)
        enddo

!-----------------------------------------------------------------------
!       calculate lateral diffusion and advection of tracers
!-----------------------------------------------------------------------

# if defined obc_west

        if (i .eq. 2) then

!         radiation condition at the western wall
!         and diffusion (del**2) of tracer

          do j=js,je
            do k=1,km
              jrow = j + joff
              ADV_Tx_obc(k,j) = -(c1w(k,j,n)+uad(k,j))
     &                  *(t(i+1,k,j,n,tau)-t(i,k,j,n,tau))/dxu(i)
            enddo
          enddo
        endif
# endif
# if defined obc_east

        if (i .eq. imtm1) then

!         radiation condition at the eastern wall
!         and diffusion (del**2) of tracer

          do j=js,je
            do k=1,km
              jrow = j + joff
              ADV_Tx_obc(k,j)=-(c1e(k,j,n)+uad(k,j))
     &                  *(t(i,k,j,n,tau)-t(i-1,k,j,n,tau))/dxu(i-1)
            enddo
          enddo
        endif
# endif

!-----------------------------------------------------------------------
!       construct tracer source terms here
!       force n/s wall to observed values
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km
            source(i,k,j) = c0
          enddo
        enddo

# if !defined orlanski
        call obcsponge2 (joff, js, je, i, n, t(1,1,1,1,tau),
     &                  uad, source)
# endif

!-----------------------------------------------------------------------
!       calculate the new tracer quantities allowing for implicit
!       treatment of vertical diffusion
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km
            twodt(k) = dtts*dtxcel(k)
            t(i,k,j,n,taup1) = t(i,k,j,n,tau) + twodt(k)*
     &                        (DIFF_Tx(i,k,j) + DIFF_Tz(i,k,j)
     &                       + DIFF_Ty_obc(i,k,j)
     &                       + ADV_Tx_obc(k,j)
     &                       + source(i,k,j)
     &                        )*tmask(i,k,j)
          enddo
        enddo

# if defined trace_indices
        if (n.eq.1) write (stdout,'(2x,6(a,i4))')
     &    "=> In  cobc2: js=",js," je=",je," joff=",joff
     &,   " jrows=",js+joff," to ",je+joff, " at i= ",i
# endif

#if defined implicitvmix || defined isopycmix

!-----------------------------------------------------------------------
!       add dT/dt component due to implicit vertical diffusion
!-----------------------------------------------------------------------

        call ivdift (joff, js, je, i, i, n, twodt)
#endif

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, js, je, i, i, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

!-----------------------------------------------------------------------
!     construct diagnostics (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, js, je, i, i, idiag)

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------

      call asbct (joff, j, j, istrt, iend, isst, 1)
      call asbct (joff, j, j, istrt, iend, isss, 2)

      return
      end

# if !defined orlanski
      subroutine obcsponge2 ( joff, js, je, i, n, tm1, uad, source )

!=======================================================================
!     newtonian damping variables for obc regions. damp to prescribed
!     values if there are inflow conditions. data must be prepared
!     using the "mkobc"routines included in the programs for working
!     with the MOM dataset.

!     input:
!       i    = (2,imt-1) open boundary row
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       n    = (1,2) = (T,S) tracer component
!       tm1  = tracer at "tau"
!       uad  = advective vel. for tracers at open walls
!            ==> restore to boundary values for inflow conditions

!     output:
!       source = newtonian damping term

!     based on code by: A. Biastoch
!=======================================================================

# include "param.h"
# include "iounit.h"
# include "obc_data.h"
# include "switch.h"
# include "tmngr.h"

      dimension source(imt,km,jsmw:jemw)
      dimension tm1(imt,km,jmw,nt)
      dimension uad(km,jmt)

      if ((joff + js .eq. 2) .and. n .eq. 1) then

!-----------------------------------------------------------------------
!       decide whether to read sponge data or not
!-----------------------------------------------------------------------

        begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
        methodobc = 3
        call timeinterp (begtim, indxob2, tobc2, obc2dpm, 12, .true.
     &,           methodobc, inextdobc2, iprevdobc2, wprevobc2,
     &            readob2, inextobc2, iprevobc2)

!       read in the next data record from disk when needed

        if (readob2) then
          call getunit (ionew6, 'obc2', opt_obc2)
          read (ionew6, rec=inextdobc2) obctnext, spdpmn, im, kk, jm
# if defined obc_west
     &,  obcw,(((obbuf_west (j,k,nn,inextobc2),j=1,jmt),k=1,km),nn=1,2)
# endif
# if defined obc_east
     &,  obce,(((obbuf_east (j,k,nn,inextobc2),j=1,jmt),k=1,km),nn=1,2)
# endif
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obc2 record =',inextdobc2,' into buffer =',inextobc2
     &,   ' method #',methodobc,' at day =', begtim, stamp

        call relunit (ionew6)
        endif
      endif
      if (n .le. 2) then

!-----------------------------------------------------------------------
!       construct newtonian damping term using obc data
!-----------------------------------------------------------------------

# if defined obc_west
          if (i .eq. 2 .and. obcw .ne. c0) then
            tnext = c1-wprevobc2
            do j=js,je
              jrow = joff + j
              do k=1,km
                if (uad(k,j) .eq. 0.) then
                   data = tnext    *obbuf_west(jrow,k,n,inextobc2)
     &                  + wprevobc2*obbuf_west(jrow,k,n,iprevobc2)
                   source(i,k,j) = source(i,k,j) -
     &                               obcw*(tm1(i,k,j,n) - data)
                endif
              enddo
            enddo
          endif
# endif

# if defined obc_east
          if (i .eq. imtm1 .and. obce .ne. c0) then
            tnext = c1-wprevobc2
            do j=js,je
              jrow = joff + j
              do k=1,km
                if (uad(k,j) .eq. 0.) then
                   data = tnext    *obbuf_east(jrow,k,n,inextobc2)
     &                  + wprevobc2*obbuf_east(jrow,k,n,iprevobc2)
                   source(i,k,j) = source(i,k,j) -
     &                               obce*(tm1(i,k,j,n) - data)
                endif
              enddo
            enddo
          endif
# endif
      endif

      return
      end

# endif

#else
      subroutine cobc2 (joff, js, je, i)
      return
      end
#endif
