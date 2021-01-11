!====================== include file "tidal_kv.h" ========================

!     quantities for parameterization of tidal mixing
# if defined uvic_tidal_kv

      real zetar, ogamma, gravrho0r
      common /tdr/ zetar, ogamma, gravrho0r

!     edr = energy dissipation rate due to tides

      real tdr
      common /tdr/ edr(imt,jmt)
#endif
