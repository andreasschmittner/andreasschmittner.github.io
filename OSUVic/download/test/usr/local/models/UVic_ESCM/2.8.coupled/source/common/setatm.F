      subroutine setatm

!-----------------------------------------------------------------------
!     set up for simple time varying data atmosphere (just for sbc)

!     based on code by: R. C. Pacanowski and M. Eby
!-----------------------------------------------------------------------
#if defined time_varying_sbc_data

      implicit none

!     this is set up only for periodic monthly data.

!     before adding a surface boundary condition
!      1. check numsbc and ntdbc are large enough (in csbc.h & ctdbc.h)
!      2. add index definition to UVic_ESCM (eq. itaux)
!      3. add code below and to atmos.F

      integer i, id, im, iou, j, k, m, n

# include "param.h"
# include "csbc.h"
# if defined uvic_embm
#  include "atm.h"
# endif
# include "calendar.h"
# include "ctdbc.h"
# include "tmngr.h"
# include "switch.h"

      real c10, c100, C2K, p001, p035, realdays

      c10 = 10.
      c100 = 100.
      C2K = 273.15
      p001 = 0.001
      p035 = 0.035

      do n=1,ntdbc
        ntdrec(n) = 12
        period(n) = .true.
        do m=1, ntdrec(n)
          k = m + 1
          aprec(m,n) = daypm(m)
          if (.not. eqyear .and. nint(aprec(m,n)) .eq. 28) then
            aprec(m,n) = aprec(m,n) + 0.2425
            print*, '=>Warning: adding 0.2425 days to feb for leap year'
          endif
!         create time stamp for the end of each month
          if (k .gt. 12) then
            call mkstmp (dstamp(m,n), 2, 1, 1, 0, 0, 0)
          else
            call mkstmp (dstamp(m,n), 1, k, 1, 0, 0, 0)
          endif
        enddo

        call timeinterpi (ntdrec(n), dstamp(1,n), aprec(1,n), tdrec(1,n)
     &,                  isbcstart(n), period(n))

        call addtime (initial, imodeltime, itemptime)
        call subtime (itemptime, isbcstart(n), itemptime)
        daysbc(n) = realdays(itemptime)
        iprevm(n) = 1
        inextm(n) = 2
        method    = 3
        call timeinterp (daysbc(n), n, tdrec(1,n), aprec(1,n), ntdrec(n)
     &,      period(n), method, inextd(n), iprevd(n), wprev(n)
     &,      rdtdbc(n), inextm(n), iprevm(n))
      enddo

      euler2 = .false.

!     load previous and next data for all boundary conditions
# if defined uvic_sbc_in_memory
!     this is set up only for monthly data.
      do k=1,12
# else
      do k=1,2
# endif
        n = 1
        do m=1,numsbc

          if (n .le. ntdbc) then

            rdtdbc(n) = .true.
# if defined uvic_sbc_in_memory
            id = k
            im = k
# else
            id = iprevd(n)
            im = iprevm(n)
            if (k .eq. 2) then
              id = inextd(n)
              im = inextm(n)
           endif
# endif

            if ( m .eq. itaux ) then
              if (k .eq. 1) print*, 'x component of windstress'
              call get_sbc (n, 'taux_mth.nc', 'taux', id, im, c10, c0)

            elseif ( m .eq. itauy ) then
              if (k .eq. 1) print*, 'y component of windstress'
              call get_sbc (n, 'tauy_mth.nc', 'tauy', id, im, c10, c0)

            elseif ( m .eq. iws ) then
              if (k .eq. 1) print*, 'surface wind speed'
              call get_sbc (n, 'ws_mth.nc', 'ws', id, im, c100, c0)

            elseif ( m .eq. iaca ) then
              if (k .eq. 1) print*, 'atmospheric coalbedo'
              call get_sbc (n, 'a_calb_mth.nc', 'a_calb', id, im, c1
     &,         c0)

# if defined uvic_read_my_stf
              if (k .eq. 1) print*, 'heat flux'
            elseif ( m .eq. ihflx ) then
              call get_sbc (n, 'hflx_mth.nc', 'hflx', id, im, c1, c0)

            elseif ( m .eq. isflx ) then
              if (k .eq. 1) print*, 'salt flux'
              call get_sbc (n, 'sflx_mth.nc', 'sflx', id, im, c1, c0)

# elif defined uvic_replacst
            elseif ( m .eq. isst ) then
              if (k .eq. 1) print*, 'sea surface temperature'
              call get_sbc (n, 'sst_mth.nc', 'sst', id, im, c1, -C2K)

            elseif ( m .eq. isss ) then
              if (k .eq. 1) print*, 'sea surface salinity'
              call get_sbc (n, 'sss_mth.nc', 'sss', id, im, p001, -p035)

# elif defined restorst
            elseif ( m .eq. ihflx ) then
              if (k .eq. 1) print*, 'data for restoring sst'
              call get_sbc (n, 'sst_mth.nc', 'sst', id, im, c1, -C2K)

            elseif ( m .eq. isflx ) then
              if (k .eq. 1) print*, 'data for restoring sss'
              call get_sbc (n, 'sss_mth.nc', 'sss', id, im, p001, -p035)

# endif
# if defined uvic_embm_awind &&  defined uvic_embm
            elseif ( m .eq. iwa ) then
              if (k .eq. 1) print*, 'surface wind angle'
              call get_sbc (n, 'wa_mth.nc', 'wa', id, im, c1, c0)

# endif
# if defined uvic_embm_adv_q &&  defined uvic_embm
            elseif ( m .eq. iwxq ) then
              if (k .eq. 1) print*, 'x component of advecting wind'
              call get_sbc (n, 'wx_mth.nc', 'wx_q', id, im, c100, c0)

            elseif ( m .eq. iwyq ) then
              if (k .eq. 1) print*, 'y component of advecting wind'
              call get_sbc (n, 'wy_mth.nc', 'wy_q', id, im, c100, c0)

# endif
# if defined uvic_embm_adv_t &&  defined uvic_embm
            elseif ( m .eq. iwxt ) then
              if (k .eq. 1) print*, 'x component of advecting wind'
              call get_sbc (n, 'wx_mth.nc', 'wx_t', id, im, c100, c0)

            elseif ( m .eq. iwyt ) then
              if (k .eq. 1) print*, 'y component of advecting wind'
              call get_sbc (n, 'wy_mth.nc', 'wy_t', id, im, c100, c0)

# endif
# if defined shortwave
            elseif ( m .eq. ipsw ) then
              if (k .eq. 1) print*, 'penetrating shortwave radiation'
              call get_sbc (n, 'psw_mth.nc', 'psw', id, im, c1, c0)

# endif
# if defined uvic_mtlm
            elseif ( m .eq. idtr ) then
              if (k .eq. 1) print*, 'diurnal temperature range'
              call get_sbc (n, 'dtr_mth.nc', 'dtr', id, im, c1, c0)

# endif
            endif

          endif

        enddo
      enddo

      return
      end

      subroutine get_sbc (n, file, name, id, im, scalar, offset)

      implicit none

# include "param.h"
# include "csbc.h"
# include "ctdbc.h"

      character (*) :: file, name
      character(120) :: fname, new_file_name
      integer id, im, iou, n, ib(10), ic(10)
      logical exists
      real offset, scalar

      fname = new_file_name (file)
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'get_sbc in setatm.f'
      else
        ib(:) = 1
        ib(3) = id
        ic(:) = imt
        ic(2) = jmt
        ic(3) = 1
        call openfile (fname, iou)
        call getvara (name, iou, imt*jmt, ib, ic, obc(1,1,n,im)
     &,   scalar, offset)
        call closefile (iou)
        where (obc(:,:,:,:) .gt. 1.e30) obc(:,:,:,:) = 0.
      endif
      n = n + 1
#endif

      return
      end
