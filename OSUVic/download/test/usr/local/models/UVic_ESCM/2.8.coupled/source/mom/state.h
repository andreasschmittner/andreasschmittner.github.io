!====================== include file "state.h" =========================

!     variables for equation of state

!     to = reference temperture for level
!     to = reference salinity for level
!     ro0= reference density for level
!     c  = polynomial coefficients for equation of state
!     tmink = min temperature at level k used for polynomial coeffs
!     tmaxk = max temperature at level k used for polynomial coeffs
!     smink = min salinity at level k used for polynomial coeffs
!     smaxk = max salinity at level k used for polynomial coeffs

      common /stater/ ro0(km), to(km), so(km), c(km,9)
      common /stater/ tmink(km), tmaxk(km), smink(km), smaxk(km)
