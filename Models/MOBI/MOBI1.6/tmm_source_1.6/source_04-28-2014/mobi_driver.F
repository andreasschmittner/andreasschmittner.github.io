#include "MOBI_TMM_OPTIONS.h"

      program driver

      implicit none

#include "size.h"
#include "mw.h"

      real z_loc(km), drF_loc(km), dtts_loc
      real tlat_loc, day_loc, relyr_loc
      real temp_loc(km), salt_loc(km)
      real TR_loc(km,nsrc), JTR_loc(km,nsrc)
      
      integer dum1, kmt_loc
      integer n

# if defined O_carbon
      real co2atm
      real winds_loc
#endif      
#  if defined O_npzd_nitrogen
      real sg_bathy_loc(km)
#  endif
#  if defined O_npzd_fe_limitation
      real fe_loc(km) ! only the first three elements are used
#  endif
#  if defined O_embm
	  real swr_loc
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
      real aice_loc, hice_loc, hsno_loc
#   endif
#  endif
      
      dtts_loc = 108000.0

      open(unit=11, file='MOBI.dat', status='old')
      read(unit=11,fmt=*)kmt_loc
      read(unit=11,fmt=*)tlat_loc
      read(unit=11,fmt=*)day_loc
      read(unit=11,fmt=*)relyr_loc
      read(unit=11,fmt=*)temp_loc(1:kmt_loc)
      read(unit=11,fmt=*)salt_loc(1:kmt_loc)
	  do n=1,nsrc
	    read(unit=11,fmt=*)dum1,TR_loc(1:kmt_loc,n)
	    if (dum1.ne.(n+2)) then
	      print*,'ERROR!'
	      stop
	    endif
	  enddo
	  read(unit=11,fmt=*)co2atm
      read(unit=11,fmt=*)winds_loc
      read(unit=11,fmt=*)sg_bathy_loc(1:kmt_loc)
      read(unit=11,fmt=*)fe_loc(1:3)
      read(unit=11,fmt=*)swr_loc
      read(unit=11,fmt=*)aice_loc
      read(unit=11,fmt=*)hice_loc
      read(unit=11,fmt=*)hsno_loc
      close(11)

CSPK these should be in cm      
      open(unit=11, file='grid.dat', status='old')
      read(unit=11,fmt=*)z_loc(1:km)
      read(unit=11,fmt=*)drF_loc(1:km)      
      close(11)

      print*,'kmt_loc=',kmt_loc
      print*,'drF_loc=',drF_loc(:)
      print*,'z_loc=',z_loc(:)

      print*,temp_loc(:)
      print*,salt_loc(:)
      print*,TR_loc(:,:)
      print*,co2atm
      
      call MOBI_INI(z_loc, drF_loc, dtts_loc, S_surf_glob, 
     &              TR_surf_glob)

      write(*,*)'done ini'
            
      call MOBI_CALC(kmt_loc, tlat_loc, day_loc, relyr_loc, 
     &     temp_loc, salt_loc, TR_loc,     
# if defined O_carbon
     &     co2atm, winds_loc,
#endif      
#  if defined O_npzd_nitrogen
     &     sg_bathy_loc
#  endif
#  if defined O_npzd_fe_limitation
     &     fe_loc,
#  endif
#  if defined O_embm
     &     swr_loc,
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
     &     aice_loc, hice_loc, hsno_loc,
#   endif
#  endif
     &  JTR_loc)

	  do n=1,nsrc
		print*,'n=',n	  
		print*,mapt(n+2),JTR_loc(1:kmt_loc,n)
	  enddo
      
      end
