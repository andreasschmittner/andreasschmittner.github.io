!====================== include file "obc_data.h" ======================

!     newtonian damping variables for open boundary regions for use in
!     limited domain basins. data must be prepared using the "mkobc"
!     routines included in PREP_DATA.

!     obcs   = coefficient for  damping T & S back to southern boundary
!     obcn   = coefficient for  damping T & S back to northern boundary
!     obcw   = coefficient for  damping T & S back to western  boundary
!     obce   = coefficient for  damping T & S back to eastern  boundary

!     iprevdobc = pointer to previous month on disk
!              (the  month whose mid point has just been passed)
!     inextdobc = pointer to next month on disk
!              (the next month whose mid point hasn`t been reached yet)
!     iprevobc  = pointer to memory buffer for the previous month data
!              (the  month whose mid point has just been passed)
!     inextobc  = pointer to memory buffer for the next month data
!              (the next month whose mid point hasn`t been reached yet)
!     obbuf_south = buffer for holding previous and next month disk data
!                   T(imt,km,1), S(imt,km,2) for the southern boundary
!     obbuf_north = buffer for holding previous and next month disk data
!                   T(imt,km,1), S(imt,km,2) for the northern boundary
!     obbuf_west  = buffer for holding previous and next month disk data
!                   T(jmt,km,1), S(jmt,km,2) for the western boundary
!     obbuf_east  = buffer for holding previous and next month disk data
!                   T(jmt,km,1), S(jmt,km,2) for the eastern boundary
!     annlevobc = (t,f) = (replace seasonal data by annual means, use
!                       seasonal data)
!     obcdpm    = period in days for each monthly record
!     tobc      = time at midpoints of monthly records (days)
!     indxob    = index of dataset (needed for the interpolator)
!     readob    = (true,false) = (read, do not read) sponge data
!     wprevobc  = interpolation weight for previous month`s data

!     symbolically:
!     obdata(at time step) = (1-wprevobc)*obbuf(inextobc)
!                              + wprevobc*obbuf(iprevobc)

!     ..1  = northern and southern open boundaries
!     ..2  = western  and eastern  open boundaries
!     ..1p = northern and southern open boundaries, psi
!     ..2p = western  and eastern  open boundaries, psi

#if defined obc_south || defined obc_north
      common /onewti/ inextobc1,iprevobc1,inextdobc1,iprevdobc1,indxob1
      common /onewti/ inextobc1p,iprevobc1p,inextdobc1p
      common /onewti/ iprevdobc1p,indxob1p
      common /onewt/  wprevobc1,obc1dpm(12),tobc1(12)
      common /onewt/  wprevobc1p,obc1pdpm(12),tobc1p(12)
#endif
#if defined obc_west || defined obc_east
      common /onewti/ inextobc2,iprevobc2,inextdobc2,iprevdobc2,indxob2
      common /onewti/ inextobc2p,iprevobc2p,inextdobc2p
      common /onewti/ iprevdobc2p,indxob2p
      common /onewt/  wprevobc2,obc2dpm(12),tobc2(12)
      common /onewt/  wprevobc2p,obc2pdpm(12),tobc2p(12)
#endif
#if defined obc_south
      common /onewt/ obcs, obbuf_south(imt,km,2,2)
#endif
#if defined obc_north
      common /onewt/ obcn, obbuf_north(imt,km,2,2)
#endif
#if defined obc_west
      common /onewt/ obcw, obbuf_west(jmt,km,2,2)
#endif
#if defined obc_east
      common /onewt/ obce, obbuf_east(jmt,km,2,2)
#endif
      character(32) :: obcstamp, obctprev, obctnext
      common /onewtc/ obcstamp, obctprev, obctnext
      character(120) :: opt_obc1,opt_obc2,opt_obcpsi1,opt_obcpsi2
      common /onewtc/ opt_obc1,opt_obc2,opt_obcpsi1,opt_obcpsi2
      logical annlevobc, readob1,readob2,readob1p,readob2p
      common /onewtl/ annlevobc,readob1,readob2,readob1p,readob2p

!     psiwall_south = psi buffer data for southern boundary
!     psiwall_north = psi buffer data (imt) for northern boundary
!     psiwall_west  = psi buffer data (jmt) for western  boundary
!     psiwall_east  = psi buffer data (jmt) for eastern  boundary

!     jpsimax       = south of these index psi of all land masses
!                     are set to psimax
!     psimax        = streamfunction value to pass the basin

#if defined obc_south
      common /pnew/ psiwall_south(imt,2)
#endif
#if defined obc_north
      common /pnew/ psiwall_north(imt,2)
#endif
#if defined obc_west
      common /pnew/ psiwall_west(jmt,2)
#endif
#if defined obc_east
      common /pnew/ psiwall_east(jmt,2)
#endif
#if defined obc_west && defined obc_east
      common /pnew/ jpsimax, psimax
#endif
