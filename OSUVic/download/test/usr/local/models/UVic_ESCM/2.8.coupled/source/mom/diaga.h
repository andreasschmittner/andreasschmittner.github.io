!====================== include file "diaga.h" =========================

!     variables used for computing diagnostics:
#if defined uvic_carbon && defined uvic_carbon_14

!     dc14bar  = total volume weighted delta c14
!     dc14     = delta c14

      common /cdiaga_r/ dc14bar, dc14(imt,km,jsmw:jemw)
#endif
#if defined uvic_save_convection

!     totalk   = total number of levels involved in convection
!     vdepth   = ventilation depth (cm)
!     pe       = potential energy lost due to explicit convection (g/s2)

      common /cdiaga_r/ totalk(imt,jsmw:jemw), vdepth(imt,jsmw:jemw)
      common /cdiaga_r/ pe(imt,jsmw:jemw)
#endif
#if defined save_convection

!     convect0  = temperature before explicit convection
!     convect1  = time rate of change of temperature due to
!                 explicit convection. set to epsln over land
!                 for use in identifying land cells

      common /cdiaga_r/ excnv0(imt,km, jsmw:jemw)
      common /cdiaga_r/ excnv1(imt,km,jsmw:jemw)
#endif
#if defined save_mixing_coeff

!     ce = ceoff on east face of cell (1 is for momentum, 2 for tracers)
!     cn = ceoff on north face of cell (1 is for momentum, 2 for tracers)
!     cb = ceoff on bottom face of cell(1 is for momentum,2 for tracers)

      common /cdiaga_r/ ce(imt,km,jsmw:jemw,2)
      common /cdiaga_r/ cn(imt,km,jsmw:jemw,2)
      common /cdiaga_r/ cb(imt,km,jsmw:jemw,2)
#endif
