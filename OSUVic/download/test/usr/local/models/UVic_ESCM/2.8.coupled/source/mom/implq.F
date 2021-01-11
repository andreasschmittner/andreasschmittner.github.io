      subroutine implq (q2, q1, tdt, q2bct, j)
#if defined tcvmix && defined implicitvmix

!======================================================================
!      solve the vertical diffusion equation of tke implicitly
!      using the method of inverting a tridiagonal matrix
!      as described in richtmyer and morton.
!      dissipation is included implicitly

!      inputs:

!      q1      = right hand side terms
!      tdt     = 2 * timestep
!      q2bct   = top boundary condition
!      j       = row j

!      outputs:

!      q2      = returned solution

#include "param.h"
#include "ctcmix.h"
#include "levind.h"
#include "vmixc.h"
      common /invt/ a, b, c, d, e, f, g
      dimension a(imt,km), b(imt,km), c(imt,km), d(imt,km)
      dimension e(imt,0:kmp1), f(imt,0:kmp1), g(imt)
      dimension q1(imt,km), q2(imt,km)
      dimension q2bct(imt)

!======================================================================

      do 100 k=1,kmm1
        do 90 i=1,imt
          a(i,k) = eeq(i,k)*tdt
          c(i,k) = ffq(i,k)*tdt
          e(i,k) = tdt*dissp(i,k)
90      continue
100    continue

      do 200 k=1,kmm1
        do 190 i=1,imt
          b(i,k) = c1 + a(i,k) + c(i,k) + e(i,k)
          d(i,k) = q1(i,k)
190     continue
200   continue

      do 300 i = 1,imt
        e(i,kmm1) = c0
        f(i,kmm1) = c0
300   continue

!     b.c. at bottom (q2=q2l=0)

      do 400 i=1,imt
        kz = kmt(i,j)
        if (kz .eq. 0) goto 400
        a(i,kz) = c0
        d(i,kz) = c0
400   continue

      do 500 kk=2,km
        k = km + 1 - kk
          do 490 i=1,imt
            g(i) = c1/(b(i,k)-c(i,k)*e(i,k))
            f(i,k-1) = (d(i,k)+c(i,k)*f(i,k))*g(i)
            e(i,k-1) = a(i,k)*g(i)
490       continue
500   continue

      do 600 i=1,imt
        q2(i,1) = e(i,0)*q2bct(i) + f(i,0)
600   continue

      do 700 k=2,km
        do 690 i=1,imt
          q2(i,k) = e(i,k-1)*q2(i,k-1) + f(i,k-1)
690     continue
700   continue

#endif
      return
      end
