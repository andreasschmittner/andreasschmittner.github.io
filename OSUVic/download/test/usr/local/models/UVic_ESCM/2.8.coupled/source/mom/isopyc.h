!======================== include file "isopyc.h" ======================
#if defined isopycmix

!     isopycnal diffusion variables:

!     ahisop = isopycnal tracer mixing coefficient (cm**2/sec)
!     drodx  = d(rho)/dx local to east face of T cell
!     drody  = d(rho)/dy local to north face of T cell
!     drodz  = d(rho)/dz local to bottom face of T cell
!     Ai_e   = diffusion coefficient on eastern face of T cell
!     Ai_n   = diffusion coefficient on northern face of T cell
!     Ai_bx  = diffusion coefficient on bottom face of T cell
!     Ai_by  = diffusion coefficient on bottom face of T cell

!     fisop  = structure function for isopycnal diffusion coefficient.
!     slmxr  = reciprocal of maximum allowable slope of isopycnals for
!              small angle approximation

      real K11, K22, K33

      common /cisop/ alphai(imt,km,jmw), betai(imt,km,jmw)

# if defined full_tensor
      common /cisop/ ddxt(imt,km,jmw,2)
# else
      common /cisop/ ddxt(imt,km,jsmw:jemw,2)
# endif
      common /cisop/ ddyt(imt,km,1:jemw,2)
      common /cisop/ ddzt(imt,0:km,jmw,2)

      common /cisop/ Ai_ez(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop/ Ai_nz(imt,km,1:jemw,0:1,0:1)
      common /cisop/ Ai_bx(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop/ Ai_by(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop/ K11(imt,km,jsmw:jemw)
      common /cisop/ K22(imt,km,1:jemw)
      common /cisop/ K33(imt,km,jsmw:jemw)
      common /cisop/ ahisop, fisop(imt,jmt,km), slmxr
# if defined dm_taper
      common /cisop/ del_dm, s_dmr
# endif
# if defined full_tensor
      common /cisop/ delta_iso, s_minus, s_plus
# endif
# if defined gent_mcwilliams
!     adv_vetiso = zonal isopycnal mixing velocity computed at the
!                  center of the eastern face of the "t" cells
!     adv_vntiso = meridional isopycnal mixing velocity computed at
!                  the center of the northern face of the "t" cells
!                  (Note: this includes the cosine as in "adv_vnt")
!     adv_vbtiso = vertical isopycnal mixing velocity computed at the
!                  center of the top face of the "t" cells
!     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
!                  the bottom face of the "t" cells
!     athkdf = isopycnal thickness diffusivity (cm**2/sec)
      common /cisop/ athkdf
      common /cisop/ adv_vetiso(imt,km,jsmw:jemw)
      common /cisop/ adv_vntiso(imt,km,1:jemw)
      common /cisop/ adv_vbtiso(imt,0:km,jsmw:jemw)
      common /cisop/ adv_fbiso(imt,0:km,jsmw:jemw)
# endif

!     statement functions

      drodxe(i,k,j,ip) =    alphai(i+ip,k,j)*ddxt(i,k,j,1) +
     &                      betai(i+ip,k,j)*ddxt(i,k,j,2)
      drodze(i,k,j,ip,kr) = alphai(i+ip,k,j)*ddzt(i+ip,k-1+kr,j,1) +
     &                      betai(i+ip,k,j)*ddzt(i+ip,k-1+kr,j,2)

      drodyn(i,k,j,jq) =    alphai(i,k,j+jq)*ddyt(i,k,j,1) +
     &                      betai(i,k,j+jq)*ddyt(i,k,j,2)
      drodzn(i,k,j,jq,kr) = alphai(i,k,j+jq)*ddzt(i,k-1+kr,j+jq,1) +
     &                      betai(i,k,j+jq)*ddzt(i,k-1+kr,j+jq,2)

      drodxb(i,k,j,ip,kr) = alphai(i,k+kr,j)*ddxt(i-1+ip,k+kr,j,1) +
     &                      betai(i,k+kr,j)*ddxt(i-1+ip,k+kr,j,2)
      drodyb(i,k,j,jq,kr) = alphai(i,k+kr,j)*ddyt(i,k+kr,j-1+jq,1) +
     &                      betai(i,k+kr,j)*ddyt(i,k+kr,j-1+jq,2)
      drodzb(i,k,j,kr) =    alphai(i,k+kr,j)*ddzt(i,k,j,1) +
     &                      betai(i,k+kr,j)*ddzt(i,k,j,2)

# if defined full_tensor
      drodye(i,k,j,ip,jq) = alphai(i+ip,k,j)*ddyt(i+ip,k,j-1+jq,1) +
     &                      betai(i+ip,k,j)*ddyt(i+ip,k,j-1+jq,2)
      drodxn(i,k,j,ip,jq) = alphai(i,k,j+jq)*ddxt(i-1+ip,k,j+jq,1) +
     &                      betai(i,k,j+jq)*ddxt(i-1+ip,k,j+jq,2)
# endif

#endif
