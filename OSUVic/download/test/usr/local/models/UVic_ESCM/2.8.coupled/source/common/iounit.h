!====================== include file "iounit.h" ========================

!     i/o units and related variables

!     taum1disk = disk pointer for tau-1 latitude rows
!     taudisk   = disk pointer for tau   latitude rows
!     taup1disk = disk pointer for tau+1 latitude rows
!     kflds     = disk unit used for two dimensional fields
!     latdisk   = disk units for latitude rows (alternately pointed to
!                by taum1disk, taudisk, and taup1disk)
!     wide_open_mw = logical to indicate that the MW is fully opened.
!              if .true. then jmw = jmt and there are no latitude rows
!              on disk. Instead, they are all in the MW.
!              if .false. then jmw < jmt and all latitude rows are on
!              disk so they must be transferred between the MW and disk.

!     iodoc  = unit for documentation
!     iostab = unit for stability testing
!     iotim  = unit for time means
!     iotim1 = scratch disk (SSD) unit for accumulating time means
!     ionew1 = unit for reading sponge layer data
!     ionew2 = mirror unit of sponge layer data on SSD
#if !defined orlanski
# if defined obc_south || defined obc_north
!     ionew3 = unit for reading obc data (T,S)
!     ionew4 = mirror unit of obc data (T,S) on SSD
!     ionew7 = unit for reading obc data (psi)
!     ionew8 = mirror unit of obc data (psi) on SSD
# endif
# if defined obc_west || defined obc_east
!     ionew5 = unit for reading obc data (T,S)
!     ionew6 = mirror unit of obc data (T,S) on SSD
!     ionew9 = unit for reading obc data (psi)
!     ionew10= mirror unit of obc data (psi) on SSD
# endif
#endif

!     for the following, a control # < 0 implies that unformatted data
!     will be written to a unit selected by the i/o manager "iomngr.F"
!     and given a hardwired name (grep getunit *.F to see names)
!     and formatted data (to stdout) will be written. if a # > 0 and
!      # <> stdout, only unformatted data will be written.

!     iotavg = control # for tracer averages
!     iotmb  = control # for writing tracer meridional budget.
!     iotrmb = control # for term balances for tracer and momentum
!     ioglen = control # for writing global energetics integrals
!     iovmsf = control # for writing meridional stream function
!     iogyre = control # for writing gyre transport.
!     ioprxz = control # for writing x-z sections from latitudes
!     ioext  = control # for writing external mode (stream function)
!     iodsp  = control # for writing diagnostic surface pressure
!     iotsi  = control # for writing time step integrals
!     ioxbt  = control # for writing time averaged xbt data
!     iozmbc = control # for writing zonal mean surf boundary conditions

!     iotext    = 120 character text string for describing the details
!                 of the next unformatted data record.
!     expnam    = 60 character text string for the experiment name
!     runstamp  = 120 character text string for the run stamp

!     when writing unformatted data records in MOM, each data record is
!     proceeded by a header record which was written as:
!     write(unit) stamp, iotext, expnam
!     where stamp is a 32 character specification of the model date &
!     time corresponding to the time step when the data was written and
!     iotext is a 120 character description of what is in the
!     data record and how it is to be read. expnam is a 60 character
!     experiment name which shows which experiment wrote the data.
!     this makes it easy to decipher any unformatted output from the
!     model by using a program similar to the following:

!      program decifr

!-----------------------------------------------------------------------
!      decipher an unformatted file from MOM by showing the header
!      records. the file needs to copied to file "fort.21"
!-----------------------------------------------------------------------

!      character*32 stamp
!      character*120 iotext
!      character*60 expnam

!      iounit = 21
!      rewind iounit
!      do n=1,100000

!        read the header record

!        read (iounit, end=110) stamp, iotext, expnam
!        write (*,'(1x,a32,1x,a80)') stamp, iotext

!        skip the data record

!        read (iounit)
!      enddo
!110   continue
!      write (*,*) " => end of file on fort.",iounit
!      stop
!      end

!     note: all unformatted diagnostic MOM data is handled this way.
!     to insure that data is read properly, verify that arrays are
!     dimensioned correctly by comparing the listed variables against
!     those in the *.h files. (grep -i -n "variable" *.h) Also, most
!     data from MOM is written IEEE 32bit so it is read directly by
!     most workstations. However, when trying to read these IEEE files
!     on the CRAY, they must be assigned IEEE before being read.
!     Some diagnostic data is averaged over time before being written.
!     In these cases, the time "stamp" refers to the last time step
!     at the end of the averaging period. An averaging interval is
!     also written as part of the data. Averaging periods = zero
!     indicate instantaneous data.

      character iotext*120, expnam*60, runstamp*120
      common /iounit_c/ iotext, expnam, runstamp

      integer taum1disk, taudisk, taup1disk, latdisk, kflds
      integer iodoc, iostab, iotavg, iotmb, iotrmb, iotim, iotim1
      integer ioglen, iovmsf, iogyre, ioprxz, ioext, iodsp
      integer iotsi, iozmbc, ionew1, ionew2, ioxbt
      integer ionew3, ionew4, ionew7, ionew8
      integer ionew5, ionew6, ionew9, ionew10

      common /iounit_i/ taum1disk, taudisk, taup1disk
#if defined coarse_grained_parallelism
      common /iounit_i/ latdisk(3), kflds
#else
      common /iounit_i/ latdisk(2), kflds
#endif
      common /iounit_i/ iodoc, iostab, iotavg, iotmb, iotrmb
      common /iounit_i/ iotim, iotim1
      common /iounit_i/ ioglen, iovmsf, iogyre, ioprxz, ioext, iodsp
      common /iounit_i/ iotsi, iozmbc, ionew1, ionew2, ioxbt
#if !defined orlanski
# if defined obc_south || defined obc_north
      common /iounit_i/ ionew3, ionew4, ionew7, ionew8
# endif
# if defined obc_west || defined obc_east
      common /iounit_i/ ionew5, ionew6, ionew9, ionew10
# endif
#endif

      logical wide_open_mw
      common /iounit_l/ wide_open_mw
