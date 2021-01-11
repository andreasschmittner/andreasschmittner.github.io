!====================== include file "tcslab.h" ========================

!     slab block extensions for turbulent closure

#if !defined multitasking && !defined task
# define task
#endif
      task common /slabs/
     $                   q2(imt,km,nslabs,ntau)
     $,                  vdc(imt,km,nslabs,ntau)
     $,                  vvc(imt,km,nslabs,ntau)
     $,                  vdqc(imt,km,nslabs,ntau)
#  if defined leq
     $,                  q2l(imt,km,nslabs,ntau)
#  endif

!     q2a   = buffer area to hold the updated "n+1" q2 slab data
!     vdca  = buffer area to hold the updated "n+1" kh slab data
!     vvca  = buffer area to hold the updated "n+1" km slab data
!     vdqca = buffer area to hold the updated "n+1" kq slab data
!     q2la  = buffer area to hold the updated "n+1" q2l slab data

      task common /bufout/
     $                   q2a (imt,km)
     $,                  vdca(imt,km)
     $,                  vvca(imt,km)
     $,                  vdqca(imt,km)
#if defined leq
     $,                  q2la(imt,km)
#endif
