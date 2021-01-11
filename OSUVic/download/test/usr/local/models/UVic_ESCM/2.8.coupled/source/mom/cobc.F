#if defined obc_north || defined obc_south
      subroutine cobc (joff, j, is, ie)

!=======================================================================
!     compute tracers at "tau+1" for rows 2 and/or jmt-1
!     using open boundary conditions as in Stevens (1990).

!     input:
!       j    = boundary row in the MW
!       joff = offset relating "j" in the MW to latitude "jrow"
!            => jrow=j+joff  = boundary latitude (2 or jmt-1)
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

!     Reference: Stevens, D.P., 1990: On open boundary condidtions for
!                three dimensional primitive equation ocean circulation
!                models, Geophys. Astrophys. Fluid Dynamics, 51, 103-133

!     based on code by: R. Redler and A. Biastoch
!=======================================================================

# include "param.h"
      parameter (istrt=2, iend=imt-1)
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
      dimension ADV_Ty_obc(imt,km)
# include "fdift.h"

!     at open boundaries, laplacian formulation is used (at tau)

      DIFF_Ty_obc(i,k,j) = ahc_north(jrow)*tmask(i,k,j+1)*
     &                     (t(i,k,j+1,n,tau) - t(i,k,j,n,tau))
     &                   - ahc_south(jrow)*tmask(i,k,j-1)*
     &                     (t(i,k,j,n,tau) - t(i,k,j-1,n,tau))

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------

      jrow = j + joff

# if defined obc_south

      if (jrow .eq. 2) then

!-----------------------------------------------------------------------
!       1) compute the advective velocity "vad"
!          at the north face of the "t" grid cell

!       2) compute phase velocity at the southern boundary: c1s
!-----------------------------------------------------------------------

        var = -dyu(jrow+1)/dtts

        do k=1,km
          do i=2,imtm1

            vad(i,k) = (u(i,k,j,2,tau)*dxt(i)+u(i-1,k,j,2,tau)
     &                  *dxt(i-1))/(dxt(i)+dxt(i-1))
            if (vad(i,k) .gt. c0) vad(i,k) = c0

            do m=1,nt
              var1(i,k,m) = t(i,k,j+2,m,taum1)-t(i,k,j+1,m,taum1)
              if (var1(i,k,m) .eq. c0) then
                c1s(i,k,m) = var
              else
                c1s(i,k,m) = var*(t(i,k,j+1,m,tau)-t(i,k,j+1,m,taum1))
     &                           /var1(i,k,m)
                if (c1s(i,k,m) .gt. c0)  c1s(i,k,m) = c0
                if (c1s(i,k,m) .lt. var) c1s(i,k,m) = var
              endif
            enddo

          enddo
#  if defined cyclic
          vad(  1,k) = vad(imtm1,k)
          vad(imt,k) = vad(    2,k)
          do m=1,nt
              c1s(  1,k,m) = c1s(imtm1,k,m)
              c1s(imt,k,m) = c1s(    2,k,m)
          enddo
#  else
          vad(imt,k) = c0
          vad(  1,k) = c0
#  endif
        enddo
      endif
# endif

# if defined obc_north

      if (jrow .eq. jmtm1) then

!-----------------------------------------------------------------------
!       1) compute the advective velocity "vad"
!          at the south face of the "t" grid box

!       2) compute phase velocity at the northern boundary: c1n
!          (var1 is calculated in subroutine tracer because of
!           running out of bounds)
!-----------------------------------------------------------------------

        var = dyu(jrow-1)/dtts

        do k=1,km
          do i=2,imtm1

            vad(i,k) = (u(i,k,j-1,2,tau)*dxt(i)+u(i-1,k,j-1,2,tau)
     &                 *dxt(i-1))/(dxt(i)+dxt(i-1))
            if (vad(i,k) .lt. c0) vad(i,k) = c0

            do m=1,nt
              if (var1(i,k,m) .eq. c0) then
                 c1n(i,k,m) = var
              else
                 c1n(i,k,m) = -var*(t(i,k,j-1,m,tau)-t(i,k,j-1,m,taum1))
     &                           /var1(i,k,m)
                if (c1n(i,k,m) .lt. c0)  c1n(i,k,m) = c0
                if (c1n(i,k,m) .gt. var) c1n(i,k,m) = var
              endif
            enddo

          enddo
#  if defined cyclic
          vad(  1,k) = vad(imtm1,k)
          vad(imt,k) = vad(    2,k)
          do m=1,nt
             c1n(  1,k,m) = c1n(imtm1,k,m)
             c1n(imt,k,m) = c1n(    2,k,m)
          enddo
#  else
          vad(imt,k) = c0
          vad(  1,k) = c0
#  endif
        enddo

      endif
# endif

      do n=1,nt

!-----------------------------------------------------------------------
!       diffusive flux across eastern face of "T" cells
!       (use constant horizontal diffusion)
!-----------------------------------------------------------------------

        do k=1,km
          do i=istrt-1,iend

!          only del**2 diffusion

           ah_cstdxur(i,j) = ah*cstr(jrow)*dxur(i)
           diff_fe(i,k,j)  = ah_cstdxur(i,j)*
     &                       (t(i+1,k,j,n,tau) - t(i,k,j,n,tau))
           cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          enddo
        enddo

!-----------------------------------------------------------------------
!       diffusive flux across bottom face of "T" cells
!       use constant vertical diffusion
!-----------------------------------------------------------------------

        do k=1,kmm1
          do i=istrt,iend
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

        do i=istrt,iend
           kb              = kmt(i,jrow)
           diff_fb(i,0,j)  = stf(i,j,n)
           diff_fb(i,kb,j) = btf(i,j,n)
        enddo

!-----------------------------------------------------------------------
!       calculate meridional diffusion and advection of tracers
!-----------------------------------------------------------------------

# if defined obc_south

        if (jrow .eq. 2) then

!         radiation condition at the southern wall
!         and meridional diffusion (del**2) of tracer

          do k=1,km
            do i=istrt,iend
              ADV_Ty_obc(i,k) = -(c1s(i,k,n)+vad(i,k))
     &                  *(t(i,k,j+1,n,tau)-t(i,k,j,n,tau))/dyu(jrow)
            enddo
          enddo
        endif
# endif
# if defined obc_north

        if (jrow .eq. jmtm1) then

!         radiation condition at the northern wall
!         and meridional diffusion (del**2) of tracer

          do k=1,km
            do i=istrt,iend
              ADV_Ty_obc(i,k) = -(c1n(i,k,n)+vad(i,k))
     &                  *(t(i,k,j,n,tau)-t(i,k,j-1,n,tau))/dyu(jrow-1)
            enddo
          enddo
        endif
# endif

!-----------------------------------------------------------------------
!       construct tracer source terms here
!       force n/s wall to observed values
!-----------------------------------------------------------------------

        do k=1,km
          do i=istrt,iend
            source(i,k,j) = c0
          enddo
        enddo

# if !defined orlanski
        call obcsponge1 (j, jrow, istrt, iend, n, t(1,1,1,1,tau),
     &                   vad, source)
# endif

!-----------------------------------------------------------------------
!       calculate the new tracer quantities allowing for implicit
!       treatment of vertical diffusion
!-----------------------------------------------------------------------

        do k=1,km
          twodt(k) = dtts*dtxcel(k)
          do i=istrt,iend
            t(i,k,j,n,taup1) = t(i,k,j,n,tau) + twodt(k)*
     &                        (DIFF_Tx(i,k,j) + DIFF_Tz(i,k,j)
     &                       + DIFF_Ty_obc(i,k,j)
     &                       + ADV_Ty_obc(i,k)
     &                       + source(i,k,j)
     &                        )*tmask(i,k,j)
          enddo
        enddo

#if defined implicitvmix || defined isopycmix

!-----------------------------------------------------------------------
!       add dT/dt component due to implicit vertical diffusion
!-----------------------------------------------------------------------

        call ivdift (joff, j, j, istrt, iend, n, twodt)
#endif

        call setbcx (t(1,1,j,n,taup1), imt, km)

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, j, j, istrt, iend, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

!-----------------------------------------------------------------------
!     construct diagnostics (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, j, j, istrt, iend, idiag)

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------

      call asbct (joff, j, j, istrt, iend, isst, 1)
      call asbct (joff, j, j, istrt, iend, isss, 2)

#if defined trace_indices
      write (stdout,'(2x,3(a,i4))')
     & "=> In   cobc: j=",j," joff=",joff," jrow=",j+joff
# endif

      return
      end

# if !defined orlanski
      subroutine obcsponge1 ( j, jrow, is, ie, n, tm1, vad, source)

!=======================================================================
!     newtonian damping variables for obc regions. damp to prescribed
!     values if there are inflow conditions. data must be prepared
!     using the "mkobc" routines included in the programs for working
!     with the MOM dataset.

!     input:
!       j    = open boundary row in the MW
!       jrow = (2,jmtm1) open boundary row
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       n    = (1,2) = (T,S) tracer component
!       tm1  = tracer at "tau"
!       vad  = advective vel. for tracers at open walls
!            ==> restore to boundary values for inflow conditions

!     output:
!       source = newtonian damping term

!     based on code by: A. Biastoch
!=======================================================================

#  include "param.h"
#  include "iounit.h"
#  include "obc_data.h"
#  include "switch.h"
#  include "tmngr.h"

      dimension source(imt,km,jsmw:jemw)
      dimension tm1(imt,km,jmw,nt)
      dimension vad(imt,km)

      if (n .eq. 1 .and. is .eq. 2) then

!-----------------------------------------------------------------------
!       decide whether to read sponge data or not
!-----------------------------------------------------------------------

        begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
        methodobc = 3
        call timeinterp (begtim, indxob1, tobc1, obc1dpm, 12, .true.
     &,           methodobc, inextdobc1, iprevdobc1, wprevobc1,
     &            readob1, inextobc1, iprevobc1)

!       read in the next data record from disk when needed

        if (readob1) then
          call getunit (ionew4, 'obc1', opt_obc1)
          read (ionew4, rec=inextdobc1) obctnext, spdpmn, im, kk, jm
#  if defined obc_south
     &,  obcs,(((obbuf_south(i,k,nn,inextobc1),i=1,imt),k=1,km),nn=1,2)
#  endif
#  if defined obc_north
     &,  obcn,(((obbuf_north(i,k,nn,inextobc1),i=1,imt),k=1,km),nn=1,2)
#  endif
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obc1 record =',inextdobc1,' into buffer =',inextobc1
     &,   ' method #',methodobc,' at day =', begtim, stamp

        call relunit (ionew4)
        endif
      endif
      if (n .le. 2) then

!-----------------------------------------------------------------------
!       construct newtonian damping term using obc data
!-----------------------------------------------------------------------

#  if defined obc_south
          if (jrow .eq. 2 .and. obcs .ne. c0) then
            tnext = c1-wprevobc1
            do k=1,km
              do i=is,ie
                if (vad(i,k) .eq. 0.) then
                   data = tnext    *obbuf_south(i,k,n,inextobc1)
     &                  + wprevobc1*obbuf_south(i,k,n,iprevobc1)
                   source(i,k,j) = source(i,k,j) -
     &                               obcs*(tm1(i,k,j,n) - data)
                endif
              enddo
            enddo
          endif
#  endif

#  if defined obc_north
           if (jrow .eq. jmtm1 .and. obcn .ne. c0) then
            tnext = c1-wprevobc1
            do k=1,km
              do i=is,ie
                if (vad(i,k) .eq. 0.) then
                   data = tnext    *obbuf_north(i,k,n,inextobc1)
     &                  + wprevobc1*obbuf_north(i,k,n,iprevobc1)
                   source(i,k,j) = source(i,k,j) -
     &                               obcn*(tm1(i,k,j,n) - data)
                endif
              enddo
            enddo
          endif
#  endif
      endif

      return
      end

# endif

#else
      subroutine cobc (joff, j, is, ie)
      return
      end
#endif
