MODULE mo_gaussgrid

  USE mo_parameters

  IMPLICIT NONE

  ! ---------------------------------------------------------------
  !
  ! module *mo_gaussgrid* - quantities related to the gaussian grid.
  !
  ! ---------------------------------------------------------------

  REAL :: gw(jphgl)        !   gaussian weights.
  REAL :: gmu(jphgl)       !   mu=sin(gaussian latitudes).
  REAL :: coriol(jpgl)     !   *coriolis parameter:2*omega*mu.
  REAL :: twomu(jpgl)      !   2*mu.
  REAL :: cst(jphgl)       !   *square of cos(latitude):(1-mu**2).
  REAL :: sqcst(jpgl)      !   sqrt(1-mu**2).
  REAL :: rcsth(jphgl)     !   *half the reciprocal of *cst*:1/(2*cst).
  REAL :: racst(jphgl)     !   1./(a*(1.-mu**2).
  REAL :: sinlon(2*jpnlon) !   sin(longitude).
  REAL :: coslon(2*jpnlon) !   cos(longitude).
  REAL :: cxps             !   *-(a**4/(7.04e+37*0.00333*g))
  REAL :: cxuv             !   *-((a**3*1.43)/(7.04e+37*0.00333*g*omega))
  REAL :: cx3ps            !   * ((a**4*0.70)/(7.04e+37*g))
  REAL :: cx3uv            !   * (a**3/(7.04e+37*qmega*g))
  REAL :: rsqcst(jphgl)    !   *1./sqrt(1-mu**2)
  REAL :: budw(jpgl)       !   weights for global budgets *budw=gw/nlon* .
  REAL :: aslm(jpgl)       !   number of land points on each latitude line,
                           !     set by physc, on local domain
CONTAINS


  SUBROUTINE inigau

    ! Description:
    !
    ! Preset constants in *mo_gaussgrid*.
    !
    ! Method:
    !
    ! *inigau* is called from *setdyn*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, Jan 1999, subroutine inigau -> module mo_gaussgrid
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,   ONLY: ngl, nhgl, nlon
    USE mo_constants, ONLY: a, api, g, omega

    IMPLICIT NONE

    !  Local scalars: 
    REAL :: zc, zcma, zcst, zl, zsqcst
    INTEGER :: jgl, jlon

    !  Local arrays: 
    REAL :: zgmu(ngl), zgw(ngl)

    !  Intrinsic functions 
    INTRINSIC COS, SIN, SQRT


    !  Executable statements 

    !-- 1. Compute Gaussian latitudes and weights

    CALL gauaw(zgmu,zgw,ngl)

    DO jgl = 1, nhgl
      gw(jgl) = zgw(jgl)*.5
      gmu(jgl) = zgmu(jgl)

    !-- 2. Derive some other constants

      coriol(2*jgl-1) = 2.*omega*zgmu(jgl)
      coriol(2*jgl) = -2.*omega*zgmu(jgl)
      twomu(2*jgl-1) = 2.*zgmu(jgl)
      twomu(2*jgl) = -2.*zgmu(jgl)
      budw(2*jgl-1) = gw(jgl)/nlon
      budw(2*jgl) = gw(jgl)/nlon

      zcst = 1. - zgmu(jgl)**2
      cst(jgl) = zcst
      zsqcst = SQRT(zcst)
      sqcst(2*jgl-1) = zsqcst
      sqcst(2*jgl) = zsqcst
      rsqcst(jgl) = 1./zsqcst
      rcsth(jgl) = .5/zcst
      racst(jgl) = 1./(a*zcst)
    END DO

    DO jlon = 1, nlon               ! double size for rotated domains 
      zl = 2.*api*(jlon-1.)/nlon    ! on decomposed grid
      sinlon(jlon) = SIN(zl); sinlon(jlon+nlon) = sinlon(jlon)
      coslon(jlon) = COS(zl); coslon(jlon+nlon) = coslon(jlon)
    END DO
    zc = 7.04E+37
    zcma = 0.00333
    cxps = -(a**4/(zc*zcma*g))
    cxuv = -((a**3*1.43)/(zc*zcma*g*omega))
    cx3ps = ((a**4*0.7)/(zc*g))
    cx3uv = (a**3/(zc*g*omega))

    RETURN
  END SUBROUTINE inigau
!------------------------------------------------------------------------------

  SUBROUTINE gauaw (pa, pw, nlat)

    ! Description:
    !
    ! Compute abscissas and weights for gaussian integration.
    !
    ! Method:
    !

    USE mo_constants, ONLY: api

    IMPLICIT NONE

    !  Scalar arguments 
    INTEGER :: nlat

    !  Array arguments 
    REAL    :: pa(nlat), pw(nlat)
    ! *pa*  - array, length at least *k,* to receive abscis abscissas.
    ! *pw*  - array, length at least *k,* to receive weights.


    !  Local scalars: 
    REAL, PARAMETER :: eps = EPSILON(0.0)
    INTEGER, PARAMETER :: itemax = 20

    INTEGER :: iter, ins2, isym, jn, jgl
    REAL    :: za, zw, z, zan
    REAL    :: zk, zkm1, zkm2, zx, zxn, zldn, zmod

    !  Intrinsic functions 
    INTRINSIC ABS, COS, MOD, TAN

    !  Executable statements 

    ins2 = nlat/2+MOD(nlat,2)

    ! Find first approximation of the roots of the
    ! Legendre polynomial of degree nlat
    
    DO jgl = 1, ins2
       z = REAL(4*jgl-1)*api/REAL(4*nlat+2)
       pa(jgl) = COS(z+1./(TAN(z)*REAL(8*nlat**2)))
    END DO

    ! Computes roots and weights
    ! Perform the Newton loop
    ! Find 0 of Legendre polynomial with Newton loop

    DO jgl = 1, ins2

       za = pa(jgl)
    
       DO iter = 1, itemax+1
          zk = 0.0

          ! Newton iteration step
    
          zkm2 = 1.0
          zkm1 = za
          zx = za
          DO jn = 2, nlat
             zk = (REAL(2*jn-1)*zx*zkm1-REAL(jn-1)*zkm2)/REAL(jn)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zldn = (REAL(nlat)*(zkm1-zx*zk))/(1.-zx*zx)
          zmod = -zk/zldn
          zxn = zx+zmod
          zan = zxn
    
          ! computes weight
    
          zkm2 = 1.0
          zkm1 = zxn
          zx = zxn
          DO jn = 2,nlat
             zk = (REAL(2*jn-1)*zx*zkm1-REAL(jn-1)*zkm2)/REAL(jn)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zw = (1.0-zx*zx)/(REAL(nlat*nlat)*zkm1*zkm1)
          za = zan
          IF (ABS(zmod) <= eps) EXIT
       END DO

       pa(jgl) = zan
       pw(jgl) = zw * 2.
    
    ENDDO

!DIR$ IVDEP
!OCL NOVREC

    DO jgl = 1, nlat/2
       isym = nlat-jgl+1
       pa(isym) = -pa(jgl)
       pw(isym) =  pw(jgl)
    ENDDO

  END SUBROUTINE gauaw

END MODULE mo_gaussgrid
