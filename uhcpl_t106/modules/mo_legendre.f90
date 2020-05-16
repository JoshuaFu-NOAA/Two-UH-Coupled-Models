MODULE mo_legendre
  !======================================================
  ! This module gathers all routines and coefficents
  ! required for performing Legendre transformations.
  !======================================================

  IMPLICIT NONE

  !================
  ! Public entities
  !================
  PRIVATE
  !------------------------------------------------------------------
  ! Quantities depending on the current truncation/spatial resolution
  !   set by subroutine inileg
  !------------------------------------------------------------------
  PUBLIC :: pnm     ! Legendre coefficients (scalar transform)
  PUBLIC :: anm     ! Legendre coefficients (?)
  !------------------------------------------------
  ! Quantities depending on the truncation/spatial,
  !   specific to the latitude currently processed  
  ! ..set by legmod
  !------------------------------------------------
  PUBLIC :: pnmd    ! Modified coefficients for direct legendre transform
  PUBLIC :: anmd    ! ...
  PUBLIC :: rnmd    !
  !----------------
  ! ..set by leginv
  !----------------
  PUBLIC :: pnmi    ! Modified coefficients for inverse legendre transform
  PUBLIC :: anmi    ! ...
  PUBLIC :: pnmiuv  !
  PUBLIC :: anmiuv  !
  !------------------
  ! module procedures
  !------------------
  PUBLIC :: inileg  ! Set up polynomials needed for the Legendre transforms.
  PUBLIC :: legmod  ! Calculate modified Legendre polynomials for direct tran.
  PUBLIC :: leginv  ! Calculate modified Legendre polynomials for inverse tra.
  !============================================================================
  !==========================
  ! Declarations of variables
  !==========================
  !------------
  ! coefficents
  !------------
  REAL    ,ALLOCATABLE :: pnm(:,:) ! Legendre polinominals
  REAL    ,ALLOCATABLE :: anm(:,:)
  !--------------------------------------------------
  ! Modified legendre coefficients for one latitudude
  !   set by subroutine legmod (module mo_legendre)
  !--------------------------------------------------
  REAL, TARGET,  ALLOCATABLE :: pnmd(:)   ! direct legendre transform
  REAL, TARGET,  ALLOCATABLE :: anmd(:)
  REAL, TARGET,  ALLOCATABLE :: rnmd(:)
  !------------------------------------------------
  !   set by subroutine leginv (module mo_legendre)
  !------------------------------------------------
  REAL, TARGET,  ALLOCATABLE :: pnmi(:)   ! inverse legendre transform
  REAL, TARGET,  ALLOCATABLE :: anmi(:)
  REAL, TARGET,  ALLOCATABLE :: pnmiuv(:)
  REAL, TARGET,  ALLOCATABLE :: anmiuv(:)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE phcs(pnm, panm, kmp1, kkp1, pmu)

    ! Description:
    !
    ! *phcs* computes the values of the *Legendre polynomials and of their
    ! meridional derivatives at given latitudes for a rhomboidal truncation.
    !
    !
    ! 
    !        ^
    !   *N*  !
    !        !                            .
    !        !            +             .
    !        !          + +           .
    !        !        +   +         .
    !        !      +     +       .
    !        !    +       +     .
    !        !  +         +   .
    ! *KMAX* !+           + .
    !        !            .
    !        !          . .
    !        !        .   .
    !        !      .     .
    !        !    .       .
    !        !  .         .
    !        !________________________________________>
    !                  *MMAX*                      *M*
    !
    !
    ! Method:
    !
    ! *call* *phcs(pnm,panm,kmp1,kkp1,pmu)*
    !
    ! *pnm*    :*legendre polynomials values.
    ! *panm*   :(mu**2-1)*dpnm/dmu.
    ! *kmp1*   :mmax+1.
    ! *kkp1*   :kmax+1.
    ! *pmu*    :value at which *pnm* and *panm* are computed.
    !
    ! The *Legendre polynomials are defined as follows:
    !     * p(n,m)(mu)=sqrt((2n+1)*fact(n-m)/fact(n+m))    *
    !     *            /(fact(n)*2**(n+1))                 *
    !     *            *(1-mu**2)**(m/2)                   *
    !     *            *d**(n+m)(mu**2-1)**n/(dmu)**(n+m)  *
    !
    ! with *fact(k)=k*(k-1)**1
    !
    ! They are computed with the following numerically stable
    ! recurrence relation (Belousov,1962):
    !
    !     * p(n,m)(mu)=c(n,m)*p(n-2,m-2)                   *
    !     *           -d(n,m)*mu*p(n-1,m-2)                *
    !     *           +e(n,m)*p(n-1,m)                     *
    !
    ! with
    !     *c(n,m)=sqrt((2n+1)*(m+n-1)*(m+n-3)              *
    !     *           /(2n-3)/(m+n  )/(m+n-2))             *
    !
    !     *d(n,m)=sqrt((2n+1)*(m+n-1)*(n-m+1)              *
    !     *           /(2n-1)/(m+n  )/(m+n-2))             *
    !
    !     *e(n,m)=sqrt((2n+1)*(n-m)                        *
    !     *           /(2n-1)/(n+m))                       *
    !
    ! The derivatives *(panm)* are then computed as follows:
    !
    !     *pa(n,m)=n*f(n+1,m)*p(n+1,m)                     *
    !     *       -(n+1)*f(n,m)*p(n-1,m)                   *
    !
    ! with:
    !
    !     *f(n,m)=sqrt((n**2-m**2)/(4*n**2-1))             *
    !
    ! Results.:
    ! The *Legendre polynomials and their derivatives are stored
    ! column-wise. The following normalisation is used:
    ! Integral over [-1,+1] of *(pnm**2)* =.5
    !
    ! References:
    ! Belousov,S.L.,1962:Tables of normalised associated Legendre
    ! polynomials.(Mathematical tables series,
    ! Vol 18, PergamonPpress, New York, USA) 379pp
    !

    !  Scalar arguments 
    REAL :: pmu
    INTEGER :: kkp1, kmp1

    !  Array arguments 
    REAL :: panm(:), pnm(:)

    !  Local scalars: 
    REAL :: z2mm1, z2q2, zan, zateta, zcnm, zcos2, zcosfak, zcospar, zcostet, &
         &      zdnm, zenm, zmm1, zn, zn2, zn2m1, znm1, zp, zp2, zq, zq2m1, zsinfak, &
         &      zsinpar, zsintet, zsqp, zteta, zw, zwm2, zwm2q2, zwnm1, zwq
    INTEGER :: ik, inmax, inmaxm, ito, iton, jk, jm, jn

    !  Local arrays: 
    REAL :: ztemp(3,kmp1+kkp1)

    !  Intrinsic functions 
    INTRINSIC ACOS, COS, SIN, SQRT


    !  Executable statements 

    ztemp(:,:) = 0.

    ! 1. Initiate recurrence by computing:
    !      *p(0,0)* *p(1,1)* *pa(0,0)* *pa(1,1)*.

    inmax = kkp1 + kmp1

    zcos2 = SQRT(1.-pmu**2)
    zteta = ACOS(pmu)
    zan = 1.

    ztemp(1,1) = .5
    DO jn = 2, inmax
      zsinpar = 0.
      zcospar = 0.
      zp   = jn - 1.
      zp2  = zp*zp
      zsqp = 1./SQRT(zp2+zp)
      zan  = zan*SQRT(1.-1./(4*zp2))
      zcosfak = 1.
      zsinfak = zp*zsqp

      ik = jn
      DO jk = 1, ik, 2
        zq = jk - 1.
        zateta = (zp-zq)*zteta
        zcostet = COS(zateta)
        zsintet = SIN(zateta)
        IF (jn==jk) zcostet = zcostet*.5
        IF (jk/=1) THEN
           zcosfak = (zq-1.)/zq*(2*zp-zq+2.)/(2*zp-zq+1.)*zcosfak
           zsinfak = zcosfak*(zp-zq)*zsqp
        END IF
        zcospar = zcospar + zcostet*zcosfak
        zsinpar = zsinpar + zsintet*zsinfak
      END DO

      ztemp(1,jn) = zan*zcospar
      ztemp(2,jn-1) = zan*zsinpar
    END DO

    pnm(1) = .5
    pnm(1+kkp1) = ztemp(2,1)
    panm(1) = 0.
    panm(1+kkp1) = pmu*ztemp(2,1)

    ! 2. Complete recurrence

    ! 2.1 First 2 columns.

    DO jn = 2, kkp1
      pnm(jn) = ztemp(1,jn)
      pnm(jn+kkp1) = ztemp(2,jn)
      zn2   = 2.*jn
      zn2m1 = zn2 - 1.
      zn    = jn
      znm1  = jn - 1.
      panm(jn) = znm1*(pmu*ztemp(1,jn)-SQRT(zn2m1/(zn2-3))*ztemp(1,jn-1))
      panm(jn+kkp1) = zn*pmu*ztemp(2,jn) &
                    - SQRT(((zn2+1)*(zn**2-1.))/zn2m1)*ztemp(2,jn-1)
    END DO

    ! 2.2 Other columns.

    DO jm = 3, kmp1
      zmm1  = jm - 1.
      z2mm1 = zmm1*2
      ztemp(3,1) = SQRT(1.+1./z2mm1)*zcos2*ztemp(2,1)
      ito = (jm-1)*kkp1
      pnm(ito+1)  = ztemp(3,1)
      panm(ito+1) = zmm1*pmu*ztemp(3,1)
      inmaxm = inmax - jm

      DO jn = 2, inmaxm
        iton   = ito + jn
        znm1   = jn - 1.
        zq     = z2mm1 + znm1 - 1.
        zwm2   = zq + znm1
        zw     = zwm2 + 2
        zwq    = zw*zq
        zq2m1  = zq*zq - 1.
        zwm2q2 = zwm2*zq2m1
        zwnm1  = zw*znm1
        z2q2   = zq2m1*2
        zcnm   = SQRT((zwq*(zq-2.))/(zwm2q2-z2q2))
        zdnm   = SQRT((zwq*(znm1+1.))/zwm2q2)
        zenm   = SQRT(zwnm1/((zq+1.)*zwm2))
        ztemp(3,jn) = zcnm*ztemp(1,jn) - pmu*(zdnm*ztemp(1,jn+1)-zenm* &
                      ztemp(3,jn-1))
        pnm(iton)  = ztemp(3,jn)
        panm(iton) = (zmm1+znm1)*pmu*ztemp(3,jn) - SQRT(zwnm1*(zq+1.)/zwm2)* &
                      ztemp(3,jn-1)
      END DO

      DO jn = 1, inmax
        ztemp(1,jn) = ztemp(2,jn)
        ztemp(2,jn) = ztemp(3,jn)
      END DO

    END DO

    RETURN
  END SUBROUTINE phcs

!------------------------------------------------------------------------------
  SUBROUTINE inileg

    ! Description:
    !
    ! Set up polynomials needed for the Legendre transforms.
    !
    ! Method:
    !
    ! *pointer* dimensions made dynamic, and i/o unit numbers changed
    ! to comply with the i/o scheme.
    !
    ! This subroutine computes, normalises and writes suitably
    ! modified *Legendre polynomials needed for the *Legendre
    ! transforms.
    !
    ! *inileg* is called from *control*
    !
    ! Externals:
    !
    ! *phcs*      called to compute the *Legendre polynomials
    !             and their meridional derivatives.
    ! *reord*     reorders the *Legendre polynomials
    !             and their meridional derivatives.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, March 1982, original source
    ! J. K. Gibson, ECMWF, April 82, changed
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, September 1999, changes for parallel version
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,        ONLY: nhgl, nmp1, nnp1
    USE mo_gaussgrid,      ONLY: cst, gmu
    USE mo_constants,      ONLY: a
    USE mo_decomposition,  ONLY: dc => local_decomposition

    !  Local scalars: 
    REAL :: za2, zmu, zrcst
    INTEGER :: jgl

    !  Local arrays: 
    REAL :: zanm(dc%lnsp), zanmt(nmp1*nnp1), zpnm(dc%lnsp), zpnmt(nmp1*nnp1)


    !  Executable statements 

    ALLOCATE (pnm(dc%lnsp,nhgl)); pnm=0.0
    ALLOCATE (anm(dc%lnsp,nhgl)); anm=0.0

    ALLOCATE (pnmd   (dc%lnsp))
    ALLOCATE (anmd   (dc%lnsp))
    ALLOCATE (rnmd   (dc%lnsp))
    ALLOCATE (pnmi   (dc%lnsp))
    ALLOCATE (anmi   (dc%lnsp))
    ALLOCATE (pnmiuv (dc%lnsp))
    ALLOCATE (anmiuv (dc%lnsp))

    ! 1. Set up some constants and initiate scan

    za2 = a**2

    DO jgl = 1, nhgl
      zmu   = gmu(jgl)
      zrcst = 1. / cst(jgl)

      ! 2. Compute, reorder and normalise Legendre
      !    polynomials and their meridional derivatives.

      ! 2.1 Compute *Legendre polynomials

      CALL phcs(zpnmt, zanmt, nmp1, nnp1, zmu)

      ! *phcs* computes polynomials for a rhomboidal truncation.

      ! 2.2 Reorder *Legendre polynomials

      CALL reord(zpnm, zpnmt, zanm, zanmt)

      ! 2.3 Normalise *Legendre polynomials

      zpnm(:) = zpnm(:) * 2.0

      ! *anm* is divided by -(1.-mu**2) in order to get *d(pnm)/dmu.*

      zanm(:) = -zanm(:) * 2.0 * zrcst

      ! 3. Write Legendre polynomials to module fields

      pnm(:,jgl) = zpnm(:)
      anm(:,jgl) = zanm(:)

    END DO

  END SUBROUTINE inileg

!------------------------------------------------------------------------------
  SUBROUTINE reord(ppnm, ptp, phnm, pth)

    ! Description:
    !
    ! Reorders the *Legendre polynomials and their derivatives generated
    ! on a rhomboidal truncation by *phcs* to get them on the truncation
    ! used by the model.
    !
    ! Method:
    !
    ! *reord* is called from *inicom*. The polynomials are
    ! received *(ptp,pth)* and returned *(ppnm,phnm)* as subroutine
    ! arguments. The other parameters used are obtained from modules.
    !
    ! Results:
    ! *pnm* and *hnm* are returned in the natural order (i.e
    ! column wise ) on the truncation used by the model.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, March 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_control,       ONLY: nnp1
    USE mo_decomposition, ONLY: dc => local_decomposition

    !  Array arguments 
    REAL :: phnm(:), ppnm(:), pth(:), ptp(:)

    !  Local scalars: 
    INTEGER :: imp, inp, jm, jn
    INTEGER :: j


    !  Executable statements 

    ! 1. Reorders polynomials for a pentagonal truncation

    DO j = 1, dc%nlm    ! loop over zonal wavenumbers on this pe
      jm   = dc%lm(j)   ! zonal wave number
      inp  = dc%nlnp(j) ! lengt of m column
      imp  = dc%nlmp(j) ! offset of m column on this pe
      DO jn = 1, inp
        ppnm(imp+jn) = ptp((jm)*nnp1+jn)
        phnm(imp+jn) = pth((jm)*nnp1+jn)
      END DO
    END DO

  END SUBROUTINE reord
!------------------------------------------------------------------------------
  SUBROUTINE legmod(kpass, pnmt, anmt, rnmt)

    ! Description:
    !
    ! Calculate modified Legendre polynomials for direct tran.
    ! (grid-point to spherical harmonics)
    !
    ! Method:
    !
    ! *legmod* computes the modified Legendre polynomials 'pnmd', 
    ! 'anmd', 'rnmd' (one latitudinal index if kpass is present) for direc
    ! Legendre transforms using the normalised Legendre polynomials
    ! 'pnm', 'anm' stored in module 'mo_legendre'
    !
    ! *call* legmod(kpass)
    !
    ! *kpass*     The value of the main loop control variable
    !             (half latitudinal index)
    !             when *legmod* is called
    !
    ! Authors:
    !
    ! D. W. Dent, ECMWF, May 1984, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A.Rhodin        MPI, Sep 1999, changes for parallel version
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_constants,     ONLY: a
    USE mo_gaussgrid,     ONLY: gw
    USE mo_decomposition, ONLY: dc => local_decomposition

    !  Arguments 
    INTEGER , INTENT(in)  , OPTIONAL :: kpass      ! latitude index
    REAL    , INTENT(out) , OPTIONAL :: pnmt(:,:)  ! legendre coeffitients
    REAL    , INTENT(out) , OPTIONAL :: anmt(:,:)  !   present if kpass
    REAL    , INTENT(out) , OPTIONAL :: rnmt(:,:)  !   is not present

    !  Local scalars: 
    REAL    :: fnnp, z2w, zsa2
    INTEGER :: imp, innp, inp, jn, j

    !  Intrinsic functions 
    INTRINSIC REAL


    !  Executable statements 

    !-- 1. Compute modified Legendre polynomials

    !      *pnmd*=2*w*pnmd
    !      *anmd*=2*w*anmd
    !      *rnmd*=2*w*(-n(n+1)/a**2)*pnmd

    IF (PRESENT (kpass)) THEN

      z2w = 2.*gw(kpass)

      pnmd = pnm(:,kpass) * z2w
      anmd = anm(:,kpass) * z2w

      zsa2 = 1./a**2

      DO j = 1, dc%nlm

        inp  = dc%nlnp(j)
        imp  = dc%nlmp(j)
        innp = dc%lm  (j)
        DO jn = 1, inp
          fnnp = REAL(innp+jn)
          rnmd(imp+jn) = -zsa2*pnmd(imp+jn)*(fnnp-1.)*fnnp
        END DO
      END DO

    ELSE 

      DO j=1, dc%nlat/2    
        z2w = 2. * gw(j)
        pnmt(j,:) = pnm(:,j) * z2w
        anmt(j,:) = anm(:,j) * z2w
      END DO

      zsa2 = 1. / a**2

      DO j = 1, dc%nlm
        inp  = dc%nlnp(j)
        imp  = dc%nlmp(j)
        innp = dc%lm  (j)
        DO jn = 1, inp
          fnnp = REAL(innp+jn)
          rnmt(:,imp+jn) = -zsa2*pnmt(:,imp+jn)*(fnnp-1.)*(fnnp)
        END DO
      END DO

    ENDIF

  END SUBROUTINE legmod
!------------------------------------------------------------------------------
  SUBROUTINE leginv(kpass, pnmit, anmit, pnmiuvt, anmiuvt)

    ! Description:
    !
    ! Calculate modified Legendre polynomials for inverse tra.
    !
    ! Method:
    !
    ! *leginv* computes the modified Legendre polynomials for inver
    ! Legendre transforms using the normalised Legendre polynomials
    ! 'pnm', 'anm' stored in module 'mo_legendre'.
    !
    ! *call* *leginv(kpass)*
    !
    ! *kpass*      The value of the main loop control variable
    !              (half latitudinal index)
    !
    ! Authors:
    !
    ! D. W. Dent, ECMWF, May 1984, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A.Rhodin        MPI, Feb 2000, changes for parallel version
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_constants,     ONLY: a
    USE mo_gaussgrid,     ONLY: cst
    USE mo_decomposition, ONLY: dc => local_decomposition

    !  Arguments 
    INTEGER , INTENT(in)  , OPTIONAL :: kpass
    REAL    , INTENT(out) , OPTIONAL :: pnmit  (:,:)  !
    REAL    , INTENT(out) , OPTIONAL :: anmit  (:,:)  !
    REAL    , INTENT(out) , OPTIONAL :: pnmiuvt(:,:)  ! legendre coeffitients
    REAL    , INTENT(out) , OPTIONAL :: anmiuvt(:,:)  ! ..if kpass is not present

    !  Local scalars: 
    REAL    :: fnnp, zcst, zsa
    INTEGER :: imp, innp, inp, is, jm, jn, nhgl

    !  Intrinsic functions 
    INTRINSIC REAL


    !  Executable statements 

    !-- 1. Compute modified Legendre polynomials
    !      suitable for inverse Legendre transforms of u and v

    !      *pnmiuv*=pnmi*a*m/(n(n+1)).
    !      *anmiuv*=anmi*a*(1.-mu**2)/(n(n+1)).

    IF(PRESENT(kpass)) THEN

      pnmi = pnm(:,kpass)
      anmi = anm(:,kpass)

      zcst = cst(kpass)

      DO jm = 1, dc%nlm

        inp  = dc%nlnp(jm)
        imp  = dc%nlmp(jm)
        innp = dc%lm  (jm)

        IF (innp==0) THEN
          is = 2
        ELSE
          is = 1
        END IF

        DO jn = is, inp
          fnnp = REAL(innp+jn)
          pnmiuv(imp+jn) = pnmi(imp+jn)*a*innp / ((fnnp-1.)*fnnp)
          anmiuv(imp+jn) = anmi(imp+jn)*a*zcst / ((fnnp-1.)*fnnp)
        END DO

      END DO

      IF (dc%nlnm0 > 0) THEN
        pnmiuv(dc%nlmp(1)+1) = 0.
        anmiuv(dc%nlmp(1)+1) = 0.
      ENDIF

      !-- 2. Compute inverse legendre polynomials *anmi*=*anmi/a*

      zsa = 1./a

      anmi(:) = anmi(:)*zsa

    ELSE

      nhgl = dc%nlat/2

      DO jm = 1, dc%nlm

        inp  = dc%nlnp(jm)
        imp  = dc%nlmp(jm)
        innp = dc%lm  (jm)

        IF (innp==0) THEN
          is = 2
        ELSE
          is = 1
        END IF

        DO jn = is, inp
          fnnp = REAL(innp+jn)
          pnmiuvt(:,imp+jn) = pnm(imp+jn,:)*a*innp       / ((fnnp-1.)*fnnp)
          anmiuvt(:,imp+jn) = anm(imp+jn,:)*a*cst(:nhgl) / ((fnnp-1.)*fnnp)
        END DO
      END DO

      IF (dc%nlnm0 > 0) THEN
        pnmiuvt(:,dc%nlmp(1)+1) = 0.
        anmiuvt(:,dc%nlmp(1)+1) = 0.
      ENDIF

      !-- 2. Compute inverse legendre polynomials *anmi*=*anmi/a*

      zsa = 1. / a
      anmit(:,:) = TRANSPOSE(anm(:,:)) * zsa
      pnmit(:,:) = TRANSPOSE(pnm(:,:))

    ENDIF

  END SUBROUTINE leginv

END MODULE mo_legendre
