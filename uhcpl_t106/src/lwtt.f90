!+ longwave transmission functions
!+ $Id: lwtt.f90,v 1.11 2000/08/23 10:01:02 m214003 Exp $

SUBROUTINE lwtt(kdlon,kaer,kcfc,pga,pgb,pgc,pgd,puu,ptt,zttnc)

  ! Description:
  !
  ! Computes longwave transmission functions.
  !
  ! Method:
  !
  ! This routine computes the transmission functions for all the
  ! absorbers (h2o, uniformly mixed gases, and o3) in all six spectral
  ! intervals.
  !
  ! *lwtt* is called from *lwvn*, *lwvd*, *lwvb*
  !
  ! For reference see radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Explicit arguments :
  ! ==== inputs ===
  ! knd    :                     ; weighting index
  ! puu    : (kdlon,nua)         ; absorber amounts
  ! ==== outputs ===
  ! ptt    : (kdlon,ntra)        ; transmission functions
  !
  ! 1. Tansmission function by h2o and uniformly mixed gases are
  !    computed using pade approximants and horner's algorithm.
  ! 2. Transmission by o3 is evaluated with malkmus's band model.
  ! 3. Transmission by h2o continuum, cfc's and aerosols follow an
  !    a simple exponential decrease with absorber amount.
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, December 1988, original source
  ! R. Van Dorland, KNMI, June 1992, changed
  ! U. Schlese, DKRZ, May 1993, changed
  ! U. Schlese, DKRZ, June 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_longwave

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kaer, kcfc, kdlon

  !  Array arguments 
  REAL :: pga(kdlon,8,2), pgb(kdlon,8,2), pgc(kdlon,5,2), pgd(kdlon,5,2), &
&      ptt(kdlon,ntra), puu(kdlon,nua), zttnc(kdlon,5)

  !  Local scalars: 
  REAL :: zcoac, zeu, zpu, zsoz, zsq1, zsq2, zto1, zto2, ztoz, zuxy, zuxz, &
&      zvxy, zvxz, zx, zxd, zxi2, zxn, zy, zyi2, zz
  INTEGER :: jl, jk

  !  Intrinsic functions 
#ifdef ECLIB
  REAL :: EXPHF, SQRTHF
!DIR$ VFUNCTION EXPHF, SQRTHF
#define EXP(x)  EXPHF(x)
#define SQRT(x) SQRTHF(x)
#else
  INTRINSIC EXP, SQRT
#endif

!DIR$ NOBOUNDS

  !  Executable statements 

!-- 1.1 Horner's algorithm for h2o and co2 transmission

  DO jk=1,8
    DO jl = 1, kdlon
      zz = SQRT(puu(jl,jk))
      zxd = pgb(jl,jk,1) + zz*(pgb(jl,jk,2)+zz)
      zxn = pga(jl,jk,1) + zz*(pga(jl,jk,2))
      ptt(jl,jk) = zxn/zxd
    END DO
  END DO

!-- 1.2 Horner's algorithm for n2o and ch4 transmission

  DO jk=1,5
    DO jl = 1, kdlon
      zz = SQRT(puu(jl,jk+13))
      zxn = pgc(jl,jk,1) + zz*(pgc(jl,jk,2))
      zxd = pgd(jl,jk,1) + zz*(pgd(jl,jk,2)+zz)
      zttnc(jl,jk) = zxn/zxd
    END DO
  END DO

!-- 2. Continuum, ozone, aerosols and tracegases

  DO jl = 1, kdlon
    ptt(jl,9) = ptt(jl,8)

    ! -  Continuum absorption: e- and p-type

    zpu = puu(jl,10)
    zeu = puu(jl,11)

    ! -  Ozone absorption

    zx = puu(jl,12)
    zy = puu(jl,13)
    zuxy = 4.*zx*zx/(pialf0*zy)
    zvxy = (pialf0*zy)/(zx+zx)
    zsq1 = SQRT(1.+o1h*zuxy) - 1.
    zsq2 = SQRT(1.+o2h*zuxy) - 1.

    zxi2 = puu(jl,28)
    zyi2 = puu(jl,29)
    zuxz = 4.*zxi2*zxi2/(piaod2*zyi2)
    zsoz = SQRT(1.+savod2*zuxz) - 1.
    zvxz = piaod2*zyi2/(2.*zxi2)

    ! Interval 0-350 cm-1 + 1440-1880 cm-1

    ptt(jl,10) = EXP(-puu(jl,23))

    ! Interval 500-800 cm-1

    zcoac = 47.7*(0.017*zpu+zeu) + puu(jl,19) + puu(jl,24)
    ztoz = EXP(-zvxz*zsoz-zcoac)
    ptt(jl,11) = zttnc(jl,1)*ztoz

    ! Interval 800-970 cm-1 + 1110-1250 cm-1

    zcoac = 8.31*(0.0025*zpu+zeu) + puu(jl,20) + puu(jl,25)
    ptt(jl,12) = zttnc(jl,2)*zttnc(jl,4)*EXP(-zcoac)

    ! Interval 970-1110 cm-1

    zcoac = 5.87*(0.0018*zpu+zeu) + puu(jl,21) + puu(jl,26)
    zto1 = EXP(-zvxy*zsq1-zcoac)
    zto2 = EXP(-zvxy*zsq2-zcoac)
    ptt(jl,13) = 0.7554*zto1 + 0.2446*zto2

    ! Interval 350-500 cm-1

    zcoac = 209.*(0.059*zpu+zeu) + puu(jl,27)
    ptt(jl,14) = EXP(-zcoac)

    ! Interval 1250-1440 cm-1 + 1880-2820 cm-1

    zcoac = puu(jl,22) + puu(jl,23)
    ptt(jl,15) = zttnc(jl,3)*zttnc(jl,5)*EXP(-zcoac)
  END DO

  RETURN
END SUBROUTINE lwtt
