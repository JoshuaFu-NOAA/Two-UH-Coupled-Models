!+ inverse Legendre transforms
!+ $Id: lti.f90,v 1.23 2000/02/23 08:45:46 m214003 Exp $

!OCL NOALIAS

#ifdef CRAY
#define ddot  sdot
#define dgemm sgemm
#endif


SUBROUTINE lti

  ! Description:
  !
  ! Inverse Legendre transforms
  !
  ! Method:
  !
  ! This subroutine performs inverse *legendre transforms
  !
  ! *lti* is called from *scan2*
  !
  ! Results:
  ! *lti* computes the *fourier components
  ! for the current latitude line.
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, December 1984, original source
  ! U. Schlese, DKRZ, December 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_memory_f,      ONLY: fad, fadu0, fatp, fatpm, fau, fau0, fav, favo, &
                              fsd, fsdu0, fstp, fstpm, fsu, fsu0, fsv, fsvo
  USE mo_legendre,      ONLY: leginv
  USE mo_decomposition, ONLY: lc => local_decomposition

  IMPLICIT NONE

  !  Local loop bounds

  INTEGER          :: nllev, nllevp1, nlmp1, nlnm0
  INTEGER ,POINTER :: lm(:), nlmp(:), nlnp(:)

  !  Global bounds
  INTEGER          :: nhgl

  !  Local scalars: 

  INTEGER :: j, ik, ims, inn, inp, ins, irow, iwrk, jl, jlev, jm
  LOGICAL :: lotypa

  !  Local arrays: 

  REAL :: pnmit   (lc% nlat/2 ,lc% lnsp) ,&
          anmit   (lc% nlat/2 ,lc% lnsp) ,&
          pnmiuvt (lc% nlat/2 ,lc% lnsp) ,&
          anmiuvt (lc% nlat/2 ,lc% lnsp)

  REAL, TARGET :: faud (lc% nllev, 2, lc% nlm, lc% nlat/2), &
                  faur (lc% nllev, 2, lc% nlm, lc% nlat/2), &
                  favd (lc% nllev, 2, lc% nlm, lc% nlat/2), &
                  favr (lc% nllev, 2, lc% nlm, lc% nlat/2), &
                  fsud (lc% nllev, 2, lc% nlm, lc% nlat/2), &
                  fsur (lc% nllev, 2, lc% nlm, lc% nlat/2), &
                  fsvd (lc% nllev, 2, lc% nlm, lc% nlat/2), &
                  fsvr (lc% nllev, 2, lc% nlm, lc% nlat/2)

  REAL, POINTER :: fd (:,:,:,:), fdu0(:,:),     ftp(:,:,:,:), ftpm(:,:,:,:), fu0(:,:), &
                   fud(:,:,:,:), fur (:,:,:,:), fvd(:,:,:,:), fvo (:,:,:,:), fvr(:,:,:,:)

  !  External functions

  REAL, EXTERNAL :: ddot

  !  External subroutines 

  EXTERNAL dgemm

  !  Intrinsic functions 

  INTRINSIC MOD


  !  Executable statements 

!-- Set local loop bounds

  nllev   =  lc% nllev    ! number of levels
  nllevp1 =  lc% nllevp1  ! number of levels + 1
  nlmp1   =  lc% nlm      ! number of m wave numbers
  nlnm0   =  lc% nlnm0    ! number of coefficients for m=0
  lm      => lc% lm       ! actual local wave numbers m handled
  nlmp    => lc% nlmp     ! displacement of the first point of columns
  nlnp    => lc% nlnp     ! number of points on each column
  nhgl    =  lc% nlat/2   ! global halv number of gaussian latitudes

!-- Set initial values for transforms (zero coefficients (1,m=0 n=0)))

  IF (nlnm0 > 0) lvo(:,1,1) = 0.0
  IF (nlnm0 > 0) ld (:,1,1) = 0.0

! -- Control scan

!-- Calculate modified Legendre polynomials

    CALL leginv(pnmit=pnmit ,anmit=anmit ,pnmiuvt=pnmiuvt ,anmiuvt=anmiuvt)

!-- Inverse *Legendre transforms

  DO iwrk = 0, 1

    lotypa = MOD(iwrk,2) == 0
    IF (lotypa) THEN
      fd   => fsd  (:,:,:,:)
      ftp  => fstp (:,:,:,:)
      ftpm => fatpm(:,:,:,:)
      fvo  => fsvo (:,:,:,:)
      fu0  => fsu0 (:,:)
      fdu0 => fadu0(:,:)
      fur  => faur
      fvr  => fsvr
      fud  => fsud
      fvd  => favd
      inn  =  (nlnm0+1)/2
      ik = 1
    ELSE
      fd   => fad  (:,:,:,:)
      ftp  => fatp (:,:,:,:)
      ftpm => fstpm(:,:,:,:)
      fvo  => favo (:,:,:,:)
      fu0  => fau0 (:,:)
      fdu0 => fsdu0(:,:)
      fur  => fsur
      fvr  => favr
      fud  => faud
      fvd  => fsvd
      inn  =  nlnm0/2
      ik = 2
    END IF

    DO j = 1, nlmp1
      jm = lm(j)+1               ! m + 1
      ims  = nlmp(j) + 1         ! offset to local m columns  (spectral coef.)
      inp  = nlnp(j)             ! column length
      IF (lotypa) THEN
        ins  = (inp+1)/2
      ELSE
        ins  = inp/2
        ims  = ims  + 1
      END IF

      IF (ins>0) THEN
        IF(nllev   > 0) THEN
          CALL dgemm ('N','T',2*nllev,nhgl,ins,1.0,ld (1,1,ims),4*nllev,pnmit  (1,ims),2*nhgl,0.0,fd (1,1,j,1),2*nllev  *nlmp1)
          CALL dgemm ('N','T',2*nllev,nhgl,ins,1.0,ld (1,1,ims),4*nllev,pnmiuvt(1,ims),2*nhgl,0.0,fud(1,1,j,1),2*nllev  *nlmp1)
          CALL dgemm ('N','T',2*nllev,nhgl,ins,1.0,ld (1,1,ims),4*nllev,anmiuvt(1,ims),2*nhgl,0.0,fvd(1,1,j,1),2*nllev  *nlmp1)
          CALL dgemm ('N','T',2*nllev,nhgl,ins,1.0,lvo(1,1,ims),4*nllev,pnmit  (1,ims),2*nhgl,0.0,fvo(1,1,j,1),2*nllev  *nlmp1)
          CALL dgemm ('N','T',2*nllev,nhgl,ins,1.0,lvo(1,1,ims),4*nllev,anmiuvt(1,ims),2*nhgl,0.0,fur(1,1,j,1),2*nllev  *nlmp1)
          CALL dgemm ('N','T',2*nllev,nhgl,ins,1.0,lvo(1,1,ims),4*nllev,pnmiuvt(1,ims),2*nhgl,0.0,fvr(1,1,j,1),2*nllev  *nlmp1)
        ENDIF
        IF(nllevp1 > 0) THEN
          CALL dgemm ('N','T',2*nllevp1,nhgl,ins,1.0,ltp(1,1,ims),4*nllevp1,pnmit(1,ims),2*nhgl,0.0,ftp (1,1,j,1),2*nllevp1*nlmp1)
          CALL dgemm ('N','T',2*nllevp1,nhgl,ins,1.0,ltp(1,1,ims),4*nllevp1,anmit(1,ims),2*nhgl,0.0,ftpm(1,1,j,1),2*nllevp1*nlmp1)
        ENDIF
      ELSE
        fd  (:,:,j,:) = 0.
        fud (:,:,j,:) = 0.
        fvd (:,:,j,:) = 0.
        fvo (:,:,j,:) = 0.
        fur (:,:,j,:) = 0.
        fvr (:,:,j,:) = 0.
        ftp (:,:,j,:) = 0.
        ftpm(:,:,j,:) = 0.
      ENDIF
    END DO

!-- 2. Transform mean wind

    IF (inn > 0 ) THEN
      DO irow = 1, nhgl
        DO jlev = 1, nllev
          fu0 (jlev,irow) = ddot(inn,lu0(jlev,ik),2*nllev,pnmit(irow,ik),2*nhgl)
          fdu0(jlev,irow) = ddot(inn,lu0(jlev,ik),2*nllev,anmit(irow,ik),2*nhgl)
        END DO
      END DO
    END IF

  ! Repeat for southern row

  END DO ! iwrk = 0, 1

!-- 3. Combine rotational and divergent parts of u and v

  DO irow = 1, nhgl
    DO j = 1, nlmp1
      jm = lm(j)+1               ! m + 1
!DIR$ IVDEP
      DO jl = 1, nllev
        fsu(jl,1,j,irow) =  fsur(jl,1,j,irow) + fsud(jl,2,j,irow)
        fsu(jl,2,j,irow) =  fsur(jl,2,j,irow) - fsud(jl,1,j,irow)
        fau(jl,1,j,irow) =  faur(jl,1,j,irow) + faud(jl,2,j,irow)
        fau(jl,2,j,irow) =  faur(jl,2,j,irow) - faud(jl,1,j,irow)
        fsv(jl,1,j,irow) =  fsvr(jl,2,j,irow) - fsvd(jl,1,j,irow)
        fsv(jl,2,j,irow) = -fsvr(jl,1,j,irow) - fsvd(jl,2,j,irow)
        fav(jl,1,j,irow) =  favr(jl,2,j,irow) - favd(jl,1,j,irow)
        fav(jl,2,j,irow) = -favr(jl,1,j,irow) - favd(jl,2,j,irow)
      END DO
    END DO
  END DO

END SUBROUTINE lti
