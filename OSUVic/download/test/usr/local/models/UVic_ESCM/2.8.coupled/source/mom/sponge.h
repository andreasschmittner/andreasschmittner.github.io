!====================== include file "sponge.h" ========================

!     newtonian damping variables for sponge regions adjacent to
!     artificial southern and northern boundaries for use in limited
!     domain basins. data must be prepared using the "sponge" routines
!     included in PREP_DATA.

!     spngs  = coefficient for  damping T & S back to southern boundary
!     spngn  = coefficient for  damping T & S back to northern boundary

!     iprevd = pointer to previous month on disk
!              (the  month whose mid point has just been passed)
!     inextd = pointer to next month on disk
!              (the next month whose mid point hasn`t been reached yet)
!     iprev  = pointer to memory buffer for the previous month data
!              (the  month whose mid point has just been passed)
!     inext  = pointer to memory buffer for the next month data
!              (the next month whose mid point hasn`t been reached yet)
!     spbuf  = buffer for holding previous and next month disk data
!     spbuf  = sponge buffer data (imt,km,4,2) laid out as follows:
!              T(imt,km,1), S(imt,km,2) for the southern boundary
!              T(imt,km,3), S(imt,km,4) for the northern boundary
!     annlev = (t,f) = (replace seasonal data by annual means, use
!                       seasonal data)
!     spgdpm = period in days for each monthly record
!     tspng  = time at midpoints of monthly records (days)
!     indxsp = index of dataset (needed for the interpolator)
!     readsp = (true,false) = (read, do not read) sponge data
!     wprev  = interpolation weight for previous month`s data

!     symbolically:
!     spdata(at time step) = (1-wprev)*spbuf(inext) + wprev*spbuf(iprev)

      common /cnewti/ inext, iprev, inextd, iprevd, indxsp
      common /cnewt/ wprev, spgdpm(12), tspng(12)
      common /cnewt/ spbuf(imt,km,4,2)
      common /cnewt/ spngs(jmt), spngn(jmt)

!     array spbuf1 is only used locally to increase I/O efficiency
!     when reading and writing sponge data from disk.

      dimension spbuf1(imt,km,4)

      character(32) :: sstamp, stprev, stnext
      character(128) :: opt_sponge
      common /cnewtc/ sstamp, stprev, stnext, opt_sponge
      logical annlev, readsp
      common /cnewtl/ annlev, readsp
