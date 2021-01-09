#if defined trajectories
      subroutine ptraji

!=======================================================================
!     initialize particle positions

!     based on code by: R. C. Pacanowski
!=======================================================================

#include "param.h"
#include "coord.h"
#include "ptraj.h"
#include "iounit.h"
#include "switch.h"
#include "tmngr.h"

!-----------------------------------------------------------------------
!     distribute particles within volume defined using starting and
!     ending longitudes, latitudes, and depths given by
!     ptslon, ptelon, ptslat, ptelat. ptsdpt, and ptedpt.
!     for example: if nptraj = 500 particles then

!     pacific distribution: drop 100 particles between 180E and 190E
!                           5S to 5N, and down to 200 meters.

!      ptslon = 180.0
!      ptelon = 190.0
!      ptslat = -5.0
!      ptelat = 5.0
!      ptsdpt = 5.1e2
!      ptedpt = 200.e2
!      call pdist (1, 100)

!     deep central atlantic distribution: drop 400 particles between
!     330E and 340E, 5S to 5N and between 200 and 3000 meters depth.

!      ptslon = 330.0
!      ptelon = 340.0
!      ptslat = -5.0
!      ptelat =  5.0
!      ptsdpt = 200.1e2
!      ptedpt = 3000.e2
!      call pdist (101, nptraj)
!-----------------------------------------------------------------------

!     distribute all particles within volume defined by the following:

      ptslon = 180.0
      ptelon = 190.0
      ptslat = -5.0
      ptelat = 5.0
      ptsdpt = 5.1e2
      ptedpt = 200.e2
      call pdist (1, nptraj)
# if defined lyapunov

!-----------------------------------------------------------------------
!     initialize deformation rate matrix
!-----------------------------------------------------------------------

      do n=1,nptraj
        em(1,1,n) = c1
        em(1,2,n) = c0
        em(2,1,n) = c0
        em(2,2,n) = c1
      enddo
# endif

!-----------------------------------------------------------------------
!     write initial particle positions to file
!-----------------------------------------------------------------------

      reltim = relyr

      call getunit (io, 'particles.dta'
     &,             'unformatted sequential append ieee')

      iotext = 'read (iotraj) reltim'
      write (io) stamp, iotext, expnam
      write (io) reltim

# if defined lyapunov
      iotext ='read (iotraj) nptraj, pxyz, pijk, em'
      write (io) stamp, iotext, expnam
      write (io) nptraj, pxyz, pijk, em
      write (stdout,*)
     &' => Initial particle positions (lyapunov) written'
     &,    ' unformatted to file particles.dta on ts = ',itt,' ', stamp
# else
      iotext ='read (iotraj) nptraj, pxyz, pijk'
      write (io) stamp, iotext, expnam
      write (io) nptraj, pxyz, pijk
      write (stdout,*) ' => Initial particle positions written'
     &,    ' unformatted to file particles.dta on ts = ',itt,' ', stamp
# endif
      call relunit (io)

      return
      end

      subroutine pdist (ns, ne)

!=======================================================================
!     distribute n particles (ne-ns+1) within volume given by
!     ptslon, ptelon, ptslat, ptelat, ptsdpt, ptedpt by uniformly
!     placing approximately n**(1/3) particles along each dimension

!     based on code by: R. C. Pacanowski
!=======================================================================

# include "param.h"
# include "coord.h"
# include "ptraj.h"

      if (ns .gt. nptraj .or. ne .gt. nptraj) then
        print *, ' => Error: ns=',ns,', ne=',ne
        stop "=>pdist"
      endif

!-----------------------------------------------------------------------
!     constrain the volume (containing starting positions of
!     particles) to lie within the model domain.
!-----------------------------------------------------------------------

      if (ptslon .lt. xu(2))     ptslon = xu(2)
      if (ptslon .gt. xu(imtm1)) ptslon = xu(imtm1)
      if (ptslat .lt. yu(1))     ptslat = yu(1)
      if (ptslat .gt. yu(jmtm1)) ptslat = yu(jmtm1)
      if (ptsdpt .lt. zt(1))     ptsdpt = zt(1)
      if (ptsdpt .gt. zt(km))    ptsdpt = zt(km)

      if (ptelon .lt. xu(2))     ptelon = xu(2)
      if (ptelon .gt. xu(imtm1)) ptelon = xu(imtm1)
      if (ptelat .lt. yu(1))     ptelat = yu(1)
      if (ptelat .gt. yu(jmtm1)) ptelat = yu(jmtm1)
      if (ptedpt .lt. zt(1))     ptedpt = zt(1)
      if (ptedpt .gt. zt(km))    ptedpt = zt(km)

      if (ptslon .gt. ptelon) then
        t      = ptslon
        ptslon = ptelon
        ptelon = t
      endif
      if (ptslat .gt. ptelat) then
        t      = ptslat
        ptslat = ptelat
        ptelat = t
      endif
      if (ptsdpt .gt. ptedpt) then
        t      = ptsdpt
        ptsdpt = ptedpt
        ptedpt = t
      endif

!-----------------------------------------------------------------------
!     distribute the particles throughout the volume
!-----------------------------------------------------------------------

      cubr       = (float(ne-ns+1))**0.333333
      distx      = (ptelon - ptslon)/cubr
      disty      = (ptelat - ptslat)/cubr
      distz      = (ptedpt - ptsdpt)/cubr
      pxyz(1,ns) = p5*distx + ptslon
      pxyz(2,ns) = p5*disty + ptslat
      pxyz(3,ns) = p5*distz + ptsdpt
      do n=ns,ne
        if (n .gt. ns) then
          pxyz(1,n) = pxyz(1,n-1)
          pxyz(2,n) = pxyz(2,n-1)
          pxyz(3,n) = pxyz(3,n-1)
        endif
        pxyz(1,n)  = pxyz(1,n) + distx
        if (pxyz(1,n) .gt. ptelon) then
          pxyz(1,n) = ptslon + (pxyz(1,n)-ptelon)
          pxyz(2,n) = pxyz(2,n) + disty
          if (pxyz(2,n) .gt. ptelat) then
            pxyz(2,n) = ptslat + (pxyz(2,n)-ptelat)
            pxyz(3,n) = pxyz(3,n) + distz
          endif
        endif

        i = indp (pxyz(1,n), xu, imt)
        if (xu(i) .gt. pxyz(1,n)) then
          pijk(1,n) = i
        else
          pijk(1,n) = i+1
        endif

        jrow = indp (pxyz(2,n), yu, jmt)
        if (yu(jrow) .gt. pxyz(2,n)) then
          pijk(2,n) = jrow
        else
          pijk(2,n) = jrow+1
        endif

        k = indp (pxyz(3,n), zt, km)
        if (zt(k) .gt. pxyz(3,n)) then
          pijk(3,n) = k
        else
          pijk(3,n) = k+1
        endif
      enddo
      write (stdout,9000)
      num_part = ne-ns+1
      write (stdout,*) num_part,' particles were initialized'
     &,   ' to lie within the volume described by:'
      write (stdout,'(1x,f8.2," < lon <",f8.2)') ptslon, ptelon
      write (stdout,'(1x,f8.2," < lat <",f8.2)') ptslat, ptelat
      write (stdout,'(1x,e10.3," < dpt < ",e10.3)') ptsdpt, ptedpt
      write (stdout,*) ' '
      return
9000  format (/20x,'P A R T I C L E    I N I T I A L I Z A T I O N'/)
      end

      subroutine ptraj (j, jrow)

!=======================================================================
!     integrate particle trajectories

!     note:

!     all indicies refer to the "xu","yv" and "zt" grids.

!     it may be useful to interpolate other quantities (eg: potential
!     temperature ...) to the particle positions and save them as well.
!     this can help in understanding where and when such quantities
!     are conserved.

!     based on code by: R. C. Pacanowski
!=======================================================================

# include "param.h"
# include "coord.h"
# include "ptraj.h"
# include "grdvar.h"
# include "iounit.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
      character(30) :: fname
      save fname

!-----------------------------------------------------------------------
!     initialize so that every particle needs to be considered
!-----------------------------------------------------------------------

      if (jrow .eq. 2) then
        do n=1,nptraj
          ptdone(n) = .false.
        enddo
      endif

!-----------------------------------------------------------------------
!    calculate trajectory for all particles between jrow and jrow-1
!-----------------------------------------------------------------------

      rrad   = c1/radian
      cmdeg  = 8.982799e-8

      do n=1,nptraj
        if (.not. ptdone(n) .and. pijk(2,n) .eq. jrow) then
          ptdone(n) = .true.

!-----------------------------------------------------------------------
!         the particle is bounded by the volume with verticies given by
!         the eight nearest surrounding model grid points on the "xu",
!         "yu", and "zt" grids. (i,j,k) is the index of the deepest
!         northeast corner of this bounding volume.
!-----------------------------------------------------------------------

          i  = pijk(1,n)
!         j in the MW corresponds to  jrow = pijk(2,n)
          k  = pijk(3,n)

!-----------------------------------------------------------------------
!         compute volume weights for linear interpolation of velocity
!         at verticies of bounding volume to the particle position.

!         distances between particle and bounding volume faces

!         xe = distance to the east face
!         xw = distance to the west face
!         yn = distance to the north face
!         ys = distance to the south face
!         za = distance above to the top face
!         zb = distance below to the bottom face
!-----------------------------------------------------------------------

          xe = (xu(i) - pxyz(1,n))
          xw = (pxyz(1,n) - xu(i-1))
          yn = (yu(jrow) - pxyz(2,n))
          ys = (pxyz(2,n) - yu(jrow-1))
          za = (pxyz(3,n) - zt(k-1))
          zb = (zt(k) - pxyz(3,n))
          dv = c1/((xt(i)-xt(i-1))*(yt(jrow)-yt(jrow-1))*(zt(k)-zt(k-1)))

!-----------------------------------------------------------------------
!         construct velocity at position of particle by 3-d linear
!         interpolation.
!-----------------------------------------------------------------------

          xeyszb = xe*ys*zb*dv
          xwyszb = xw*ys*zb*dv
          xeysza = xe*ys*za*dv
          xwysza = xw*ys*za*dv
          xeynzb = xe*yn*zb*dv
          xwynzb = xw*yn*zb*dv
          xeynza = xe*yn*za*dv
          xwynza = xw*yn*za*dv

          uu = u(i-1,k-1,j,1,tau)*xeyszb   + u(i,k-1,j,1,tau)*xwyszb
     &        +u(i-1,k  ,j,1,tau)*xeysza   + u(i,k  ,j,1,tau)*xwysza
     &        +u(i-1,k-1,j-1,1,tau)*xeynzb + u(i,k-1,j-1,1,tau)*xwynzb
     &        +u(i-1,k  ,j-1,1,tau)*xeynza + u(i,k  ,j-1,1,tau)*xwynza

          vv = u(i-1,k-1,j,2,tau)*xeyszb   + u(i,k-1,j,2,tau)*xwyszb
     &        +u(i-1,k  ,j,2,tau)*xeysza   + u(i,k  ,j,2,tau)*xwysza
     &        +u(i-1,k-1,j-1,2,tau)*xeynzb + u(i,k-1,j-1,2,tau)*xwynzb
     &        +u(i-1,k  ,j-1,2,tau)*xeynza + u(i,k  ,j-1,2,tau)*xwynza

!         interpolate vertical velocities at the bases of
!         the "u" cells.

          if (pxyz(3,n) .gt. zw(k-1)) then
            za = pxyz(3,n) - zw(k-1)
            zb = zw(k) - pxyz(3,n)
            dv = c1/((xt(i)-xt(i-1))*(yt(jrow)-yt(jrow-1))
     &              *(zw(k)-zw(k-1)))
            xeyszb = xe*ys*zb*dv
            xwyszb = xw*ys*zb*dv
            xeysza = xe*ys*za*dv
            xwysza = xw*ys*za*dv
            xeynzb = xe*yn*zb*dv
            xwynzb = xw*yn*zb*dv
            xeynza = xe*yn*za*dv
            xwynza = xw*yn*za*dv
            ww = adv_vbu(i-1,k-1,j)*xeyszb   + adv_vbu(i,k-1,j)*xwyszb
     &          +adv_vbu(i-1,k  ,j)*xeysza   + adv_vbu(i,k  ,j)*xwysza
     &          +adv_vbu(i-1,k-1,j-1)*xeynzb + adv_vbu(i,k-1,j-1)*xwynzb
     &          +adv_vbu(i-1,k  ,j-1)*xeynza + adv_vbu(i,k  ,j-1)*xwynza
          else
            if (k-2 .eq. 0) then
              za = pxyz(3,n)
            else
              za = pxyz(3,n) - zw(k-2)
            endif
            zb = zw(k-1) - pxyz(3,n)
            dv = c1/((xt(i)-xt(i-1))*(yt(jrow)-yt(jrow-1))
     &              *(zw(k)-zw(k-1)))
            xeyszb = xe*ys*zb*dv
            xwyszb = xw*ys*zb*dv
            xeysza = xe*ys*za*dv
            xwysza = xw*ys*za*dv
            xeynzb = xe*yn*zb*dv
            xwynzb = xw*yn*zb*dv
            xeynza = xe*yn*za*dv
            xwynza = xw*yn*za*dv
            ww = adv_vbu(i-1,k-2,j)*xeyszb   + adv_vbu(i,k-2,j)*xwyszb
     &          +adv_vbu(i-1,k-1,j)*xeysza   + adv_vbu(i,k-1,j)*xwysza
     &          +adv_vbu(i-1,k-2,j-1)*xeynzb + adv_vbu(i,k-2,j-1)*xwynzb
     &          +adv_vbu(i-1,k-1,j-1)*xeynza + adv_vbu(i,k-1,j-1)*xwynza
          endif
# if defined lyapunov

!-----------------------------------------------------------------------
!     construct the shear and compression above and below the particle
!-----------------------------------------------------------------------

!     for du/dx and du/dy

      uxn  = (u(i,k-1,j,1,tau) - u(i-1,k-1,j,1,tau))*dxtr(i)*csu(jrow)
      uxs  = (u(i,k-1,j-1,1,tau) - u(i-1,k-1,j-1,1,tau))*dxtr(i)
     &                                                  *csu(jrow-1)
      uxn2 = (u(i,k,j,1,tau) - u(i-1,k,j,1,tau))*dxtr(i)*csu(jrow)
      uxs2 = (u(i,k,j-1,1,tau) - u(i-1,k,j-1,1,tau))*dxtr(i)*csu(jrow-1)

      uye   = (u(i,k-1,j,1,tau) - u(i,k-1,j-1,1,tau))*dytr(jrow)
      uyw   = (u(i-1,k-1,j,1,tau) - u(i-1,k-1,j-1,1,tau))*dytr(jrow)
      uye2  = (u(i,k,j,1,tau) - u(i,k,j-1,1,tau))*dytr(jrow)
      uyw2  = (u(i-1,k,j,1,tau) - u(i-1,k,j-1,1,tau))*dytr(jrow)

!     for dv/dx and dv/dy

      vxn  = (u(i,k-1,j,2,tau) - u(i-1,k-1,j,2,tau))*dxtr(i)*csu(jrow)
      vxs  = (u(i,k-1,j-1,2,tau) - u(i-1,k-1,j-1,2,tau))*dxtr(i)
     &                                                  *csu(jrow-1)
      vxn2 = (u(i,k,j,2,tau) - u(i-1,k,j,2,tau))*dxtr(i)*csu(jrow)
      vxs2 = (u(i,k,j-1,2,tau) - u(i-1,k,j-1,2,tau))*dxtr(i)*csu(jrow-1)

      vye   = (u(i,k-1,j,2,tau) - u(i,k-1,j-1,2,tau))*dytr(jrow)
      vyw   = (u(i-1,k-1,j,2,tau) - u(i-1,k-1,j-1,2,tau))*dytr(jrow)
      vye2  = (u(i,k,j,2,tau) - u(i,k,j-1,2,tau))*dytr(jrow)
      vyw2  = (u(i-1,k,j,2,tau) - u(i-1,k,j-1,2,tau))*dytr(jrow)

!     interpolate du/dx, du/dy, dv/dx, and dv/dy to particle position

      dxr   = c1/((xw+xe)*(zb+za))
      dyr   = c1/((ys+yn)*(zb+za))

      yszb = ys*zb*dyr
      ynzb = yn*zb*dyr
      ysza = ys*za*dyr
      ynza = yn*za*dyr

      xwzb = xw*zb*dxr
      xezb = xe*zb*dxr
      xwza = xw*za*dxr
      xeza = xe*za*dxr

      dudx = uxn*yszb + uxs*ynzb + uxn2*ysza + uxs2*ynza
      dudy = uye*xwzb + uyw*xezb + uye2*xwza + uyw2*xeza

      dvdx = vxn*yszb + vxs*ynzb + vxn2*ysza + vxs2*ynza
      dvdy = vye*xwzb + vyw*xezb + vye2*xwza + vyw2*xeza

!-----------------------------------------------------------------------
!     integrate the deformation matrix. note that this is not quite
!     correct when a particle encounters a boundary and "slips" along
!     the boundary.

!     the Lyapunov exponents can be computed from the eigenvalues of
!     the deformation rate matrix "em". the exponents are given by:
!     log(abs(eigen(1..2)))/T where T is the integration time.
!     let c = (em(2,2)-em(1,1))**2 + 4*(em(1,2)*em(2,1))
!     if (c >= 0.0) then
!       abs(eigen(1..2)) = abs(((em(2,2)+em(1,1))**2 + or - sqrt(c))/2)
!     else
!       abs(eigen(1)) = abs(eigen(2)) =
!       sqrt((em(2,2)+em(1,1))**2 + abs(c)))/2
!     endif
!-----------------------------------------------------------------------

      em11 = em(1,1,n) + (dudx*em(1,1,n) + dudy*em(2,1,n))*dtuv
      em12 = em(1,2,n) + (dudx*em(1,2,n) + dudy*em(2,2,n))*dtuv
      em21 = em(2,1,n) + (dvdx*em(1,1,n) + dvdy*em(2,1,n))*dtuv
      em22 = em(2,2,n) + (dvdx*em(1,2,n) + dvdy*em(2,2,n))*dtuv

      em(1,1,n) = em11
      em(1,2,n) = em12
      em(2,1,n) = em21
      em(2,2,n) = em22
# endif

!-----------------------------------------------------------------------
!         remember where the particle was
!-----------------------------------------------------------------------

          xold = pxyz(1,n)
          yold = pxyz(2,n)
          zold = pxyz(3,n)

!-----------------------------------------------------------------------
!         integrate the particle trajectory forward for one time step
!         taking convergence of meridians into account.
!-----------------------------------------------------------------------

          rcos      = cmdeg/cos(pxyz(2,n)*rrad)
          pxyz(1,n) = pxyz(1,n) + dtuv*uu*rcos
          pxyz(2,n) = pxyz(2,n) + dtuv*vv*cmdeg
          pxyz(3,n) = pxyz(3,n) - dtuv*ww

!-----------------------------------------------------------------------
!         update bottom most northeast index of bounding volume
!-----------------------------------------------------------------------

          if (pxyz(1,n) .ge. xu(i)) then
            pijk(1,n) = i + 1
          elseif (pxyz(1,n) .lt. xu(i-1)) then
            pijk(1,n) = i - 1
          endif

          if (pxyz(2,n) .ge. yu(jrow)) then
            pijk(2,n) = jrow + 1
          elseif (pxyz(2,n) .lt. yu(jrow-1)) then
            pijk(2,n) = jrow - 1
          endif

          if (pxyz(3,n) .ge. zt(k)) then
            pijk(3,n) = k + 1
          elseif (pxyz(3,n) .lt. zt(k-1)) then
            pijk(3,n) = k - 1
          endif

!-----------------------------------------------------------------------
!         do not allow any component of the trajectory to enter
!         land. If it does, reset it to its previous value
!         thereby simulating free slip conditions. hey... not perfect,
!         but beats loosing particles in land. Also if the grid has
!         isolated "T,S" cells (ones where all eight surrounding
!         velocities are on land), replace references to "tmask(i,k,j)"
!         by umask(i,k,j) + umask(i-1,k,j) + umask(i,k,j-1) +
!         umask(i-1,k,j-1) to prevent stagnation of particles if this
!         is a problem.
!-----------------------------------------------------------------------

# if defined cyclic
          if (pijk(1,n) .gt. imt) then
            if (tmask(3,k,j) .ne. 0) then
              pijk(1,n) = pijk(1,n) - (imt-2)
              i         = pijk(1,n)
              pxyz(1,n) = xu(i-1) + (pxyz(1,n)-xu(imt))
            else
              pijk(1,n) = i
              pxyz(1,n) = xold
            endif
          endif
          if (pijk(1,n) .lt. 2) then
            if (tmask(imt-2,k,j) .ne. 0) then
              pijk(1,n) = pijk(1,n) + (imt-2)
              i         = pijk(1,n)
              pxyz(1,n) = xu(i-1) + (pxyz(1,n)-xu(1))
            else
              pijk(1,n) = i
              pxyz(1,n) = xold
            endif
          endif
# endif

!-----------------------------------------------------------------------
!         constrain particles vertically to lie within ocean
!-----------------------------------------------------------------------

          if (pijk(3,n) .ne. k) then
            if (pijk(3,n) .eq. 1 .or. pijk(3,n) .gt. km) then
              pxyz(3,n) = zold
              pijk(3,n) = k
            elseif (pijk(3,n) .gt. k .and. tmask(i,k+1,j) .eq. c0) then
              pxyz(3,n) = zold
              pijk(3,n) = k
            endif
          endif

!-----------------------------------------------------------------------
!         constrain particles longitudinally to stay within ocean
!-----------------------------------------------------------------------

          if (pijk(1,n) .ne. i) then
            if (pijk(1,n) .gt. i .and. tmask(i+1,k,j) .eq. c0) then
              pxyz(1,n) = xold
              pijk(1,n) = i
            elseif (pijk(1,n) .lt. i .and. tmask(i-1,k,j) .eq. c0) then
              pxyz(1,n) = xold
              pijk(1,n) = i
            endif
          endif

!-----------------------------------------------------------------------
!         constrain particles latitudinally to stay within ocean
!-----------------------------------------------------------------------

          if (pijk(2,n) .ne. jrow) then
            if (pijk(2,n) .gt. jrow .and. tmask(i,k,j+1) .eq. c0) then
              pxyz(2,n) = yold
              pijk(2,n) = jrow
            elseif (pijk(2,n) .lt. jrow .and. tmask(i,k,j-1) .eq. c0)
     &        then
              pxyz(2,n) = yold
              pijk(2,n) = jrow
            endif
          endif
        endif
      enddo

!-----------------------------------------------------------------------
!     write particle positions
!     note: last positions are also written to "iorest" for restarting
!-----------------------------------------------------------------------

      if ((jrow .eq. jmtm1) .and. trajts) then

        reltim = relyr

        write (stdout,9000) nptraj, stamp, pxyz(1,1), pxyz(2,1)
     &,                       pxyz(3,1)*0.01

        fname    = 'particles.dta'
        call getunit (io, fname
     &,               'unformatted sequential append ieee')

        iotext = 'read (iotraj) reltim'
        write (io) stamp, iotext, expnam
        write (io) reltim

# if defined lyapunov
        iotext ='read (iotraj) nptraj, pxyz, pijk em'
        write (io) stamp, iotext, expnam
        write (io) nptraj, pxyz, pijk, em
        write (stdout,*) ' => Particle trajectories (lyapunov) written'
     &,  ' unformatted to file particles.dta on ts = ',itt, ' ', stamp
# else
        iotext ='read (iotraj) nptraj, pxyz, pijk'
        write (io) stamp, iotext, expnam
        write (io) nptraj, pxyz, pijk
        write (stdout,*) ' => Particle trajectories written'
     &,  ' unformatted to file particles.dta on ts = ',itt,' ', stamp
# endif
        call relunit (io)
        write (stdout,'(/a,a,a,i10,a,a)')
     &  ' => Particle trajectories written to file'
     &, fname,' on ts = ',itt,' ', stamp
      endif

      return
9000  format (1x,'=> Following',i6, ' particles at ',a32
     &,       '. (lon,lat,depth) of 1st particle is ('
     &,f8.2,',',f8.2,',',f8.2,'m)')
      end

#else
      subroutine ptraj (j, jrow)
      return
      end
#endif
