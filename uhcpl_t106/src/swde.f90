!+ delta-eddington in a cloudy layer
!+ $Id: swde.f90,v 1.6 1999/07/20 14:21:43 m214003 Exp $

SUBROUTINE swde(kdlon,pgg,pref,prmuz,pto1,pw,pre1,pre2,ptr1,ptr2)

  ! Description:
  !
  ! Computes the reflectivity and transmissivity of a cloudy
  ! layer using the delta-eddington's approximation.
  !
  ! Method:
  !
  ! *swde* is called by *swr*, *sw2s*
  !
  ! Explicit arguments:
  ! pgg    : (ndlon)             ; assymetry factor
  ! pref   : (ndlon)             ; reflectivity of the underlying layer
  ! prmuz  : (ndlon)             ; cosine of solar zenith angle
  ! pto1   : (ndlon)             ; optical thickness
  ! pw     : (ndlon)             ; single scattering albedo
  ! ==== outputs ===
  ! pre1   : (ndlon)             ; layer reflectivity assuming no
  !                              ; reflection from underlying layer
  ! ptr1   : (ndlon)             ; layer transmissivity assuming no
  !                              ; reflection from underlying layer
  ! pre2   : (ndlon)             ; layer reflectivity assuming
  !                              ; reflection from underlying layer
  ! ptr2   : (ndlon)             ; layer transmissivity assuming
  !                              ; reflection from underlying layer
  !
  ! Standard delta-eddington layer calculations.
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, December 1988, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kdlon

  !  Array arguments 
  REAL :: pgg(kdlon), pre1(kdlon), pre2(kdlon), pref(kdlon), prmuz(kdlon), &
&      pto1(kdlon), ptr1(kdlon), ptr2(kdlon), pw(kdlon)

  !  Local scalars: 
  REAL :: za11, za12, za13, za21, za22, za23, zalpha, zam2b, zap2b, zarg, &
&      zarg2, zb21, zb22, zb23, zbeta, zc1a, zc1b, zc2a, zc2b, zdena, zdenb, &
&      zdt, zexkm, zexkp, zexmu0, zff, zgp, zri0a, zri0b, zri0c, zri0d, zri1a, &
&      zri1b, zri1c, zri1d, zrk, zrm2, zrp, ztop, zwcp, zwm, zx1, zx2, zxm2p, &
&      zxp2p
  INTEGER :: jl

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF
!DIR$ VFUNCTION EXPHF
#define EXP(x)  EXPHF(x)
#else
  INTRINSIC EXP
#endif
  INTRINSIC MIN, SQRT


  !  Executable statements 

!-- 1. Delta-eddington calculations

  DO jl = 1, kdlon

!-- 1.1 Set up the delta-modified parameters

    zff = pgg(jl)*pgg(jl)
    zgp = pgg(jl)/(1.+pgg(jl))
    ztop = (1.-pw(jl)*zff)*pto1(jl)
    zwcp = (1-zff)*pw(jl)/(1.-pw(jl)*zff)
    zdt = 2./3.
    zx1 = 1. - zwcp*zgp
    zwm = 1. - zwcp
    zrm2 = prmuz(jl)*prmuz(jl)
    zrk = SQRT(3.*zwm*zx1)
    zx2 = 4.*(1.-zrk*zrk*zrm2)
    zrp = zrk/zx1
    zalpha = 3.*zwcp*zrm2*(1.+zgp*zwm)/zx2
    zbeta = 3.*zwcp*prmuz(jl)*(1.+3.*zgp*zrm2*zwm)/zx2
    zarg = MIN(ztop/prmuz(jl),200.)
    zexmu0 = EXP(-zarg)
    zarg2 = MIN(zrk*ztop,200.)
    zexkp = EXP(zarg2)
    zexkm = 1./zexkp
    zxp2p = 1. + zdt*zrp
    zxm2p = 1. - zdt*zrp
    zap2b = zalpha + zdt*zbeta
    zam2b = zalpha - zdt*zbeta

!-- 1.2 Without reflection from the underlying layer

    za11 = zxp2p
    za12 = zxm2p
    za13 = zap2b
    za22 = zxp2p*zexkp
    za21 = zxm2p*zexkm
    za23 = zam2b*zexmu0
    zdena = za11*za22 - za21*za12
    zc1a = (za22*za13-za12*za23)/zdena
    zc2a = (za11*za23-za21*za13)/zdena
    zri0a = zc1a + zc2a - zalpha
    zri1a = zrp*(zc1a-zc2a) - zbeta
    pre1(jl) = (zri0a-zdt*zri1a)/prmuz(jl)
    zri0b = zc1a*zexkm + zc2a*zexkp - zalpha*zexmu0
    zri1b = zrp*(zc1a*zexkm-zc2a*zexkp) - zbeta*zexmu0
    ptr1(jl) = zexmu0 + (zri0b+zdt*zri1b)/prmuz(jl)

!-- 1.3 With reflection from the underlying layer

    zb21 = za21 - pref(jl)*zxp2p*zexkm
    zb22 = za22 - pref(jl)*zxm2p*zexkp
    zb23 = za23 - pref(jl)*zexmu0*(zap2b-prmuz(jl))
    zdenb = za11*zb22 - zb21*za12
    zc1b = (zb22*za13-za12*zb23)/zdenb
    zc2b = (za11*zb23-zb21*za13)/zdenb
    zri0c = zc1b + zc2b - zalpha
    zri1c = zrp*(zc1b-zc2b) - zbeta
    pre2(jl) = (zri0c-zdt*zri1c)/prmuz(jl)
    zri0d = zc1b*zexkm + zc2b*zexkp - zalpha*zexmu0
    zri1d = zrp*(zc1b*zexkm-zc2b*zexkp) - zbeta*zexmu0
    ptr2(jl) = zexmu0 + (zri0d+zdt*zri1d)/prmuz(jl)
  END DO

  RETURN
END SUBROUTINE swde
