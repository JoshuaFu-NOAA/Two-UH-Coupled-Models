!+ computes fourier components from their symmetric and antisymmetric parts
!+ $Id: sym2.f90,v 1.15 1999/09/29 12:56:27 m214003 Exp $

SUBROUTINE sym2

  ! Description:
  !
  ! This subroutine computes *fourier components from
  ! their symmetric and antisymmetric parts and
  ! retrieves u and v from their rotational and divergent parts.
  !
  ! Method:
  !
  ! *sym2* is called from *scan1sl*.
  !
  ! Reference:
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, March 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_memory_f,      ONLY: fad, fadu0, fatp, fatpm, fau, fau0, fav, favo, &
                              fsd, fsdu0, fstp, fstpm, fsu, fsu0, fsv, fsvo
  USE mo_buffer_fft,    ONLY: lalps, ld, lt, lu, lv, lvo, ldalpsm, ldtm, &
                              lu0, ldu0, lul
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: jlev, jm, jlan, jlas, nlat
  INTEGER :: n2mp1, nlev, nlevp1, nlp2, nmp1


  !  Executable statements 

  ! Array bounds on local PE

  nmp1   = dc% nlm
  n2mp1  = dc% nlm  * 2
  nlev   = dc% nllev
  nlevp1 = dc% nllevp1
  nlp2   = dc% nlon + 2
  nlat   = dc% nlat

  ! row indices

  DO jlan=1,nlat/2         ! northern latitude or half row index
    jlas= nlat-jlan+1      ! southern latitude row index

!-- 1. Compute *Fourier components

!-- 1.1 Northern hemisphere

    DO jlev = 1, nlev

!DIR$ IVDEP
!OCL NOVREC, NOALIAS
      DO jm = 1, nmp1
        lvo (2*jm-1,jlev,jlan) = fsvo(jlev,1,jm,jlan) + favo(jlev,1,jm,jlan)
        lvo (2*jm  ,jlev,jlan) = fsvo(jlev,2,jm,jlan) + favo(jlev,2,jm,jlan)
        ld  (2*jm-1,jlev,jlan) = fsd(jlev,1,jm,jlan)  + fad(jlev,1,jm,jlan)
        ld  (2*jm  ,jlev,jlan) = fsd(jlev,2,jm,jlan)  + fad(jlev,2,jm,jlan)
        lt  (2*jm-1,jlev,jlan) = fstp(jlev,1,jm,jlan) + fatp(jlev,1,jm,jlan)
        lt  (2*jm  ,jlev,jlan) = fstp(jlev,2,jm,jlan) + fatp(jlev,2,jm,jlan)
        lu  (2*jm-1,jlev,jlan) = fsu(jlev,1,jm,jlan)  + fau(jlev,1,jm,jlan)
        lu  (2*jm  ,jlev,jlan) = fsu(jlev,2,jm,jlan)  + fau(jlev,2,jm,jlan)
        lv  (2*jm-1,jlev,jlan) = fsv(jlev,1,jm,jlan)  + fav(jlev,1,jm,jlan)
        lv  (2*jm,  jlev,jlan) = fsv(jlev,2,jm,jlan)  + fav(jlev,2,jm,jlan)
        ldtm(2*jm-1,jlev,jlan) = fstpm(jlev,1,jm,jlan)+ fatpm(jlev,1,jm,jlan)
        ldtm(2*jm  ,jlev,jlan) = fstpm(jlev,2,jm,jlan)+ fatpm(jlev,2,jm,jlan)
      END DO

    END DO

    IF (nlevp1>nlev) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jm = 1, nmp1
       lalps  (2*jm-1,jlan) = fstp (nlevp1,1,jm,jlan) + fatp (nlevp1,1,jm,jlan)
       lalps  (2*jm  ,jlan) = fstp (nlevp1,2,jm,jlan) + fatp (nlevp1,2,jm,jlan)
       ldalpsm(2*jm-1,jlan) = fstpm(nlevp1,1,jm,jlan) + fatpm(nlevp1,1,jm,jlan)
       ldalpsm(2*jm  ,jlan) = fstpm(nlevp1,2,jm,jlan) + fatpm(nlevp1,2,jm,jlan)
      END DO
    ENDIF

!-- 1.2 Southern hemisphere

    DO jlev = 1, nlev

!DIR$ IVDEP
!OCL NOVREC, NOALIAS
      DO jm = 1, nmp1
        lvo (2*jm-1,jlev,jlas) = fsvo(jlev,1,jm,jlan) - favo(jlev,1,jm,jlan)
        lvo (2*jm  ,jlev,jlas) = fsvo(jlev,2,jm,jlan) - favo(jlev,2,jm,jlan)
        ld  (2*jm-1,jlev,jlas) = fsd(jlev,1,jm,jlan)  - fad(jlev,1,jm,jlan)
        ld  (2*jm  ,jlev,jlas) = fsd(jlev,2,jm,jlan)  - fad(jlev,2,jm,jlan)
        lt  (2*jm-1,jlev,jlas) = fstp(jlev,1,jm,jlan) - fatp(jlev,1,jm,jlan)
        lt  (2*jm  ,jlev,jlas) = fstp(jlev,2,jm,jlan) - fatp(jlev,2,jm,jlan)
        lu  (2*jm-1,jlev,jlas) = fsu(jlev,1,jm,jlan)  - fau(jlev,1,jm,jlan)
        lu  (2*jm  ,jlev,jlas) = fsu(jlev,2,jm,jlan)  - fau(jlev,2,jm,jlan)
        lv  (2*jm-1,jlev,jlas) = fsv(jlev,1,jm,jlan)  - fav(jlev,1,jm,jlan)
        lv  (2*jm,  jlev,jlas) = fsv(jlev,2,jm,jlan)  - fav(jlev,2,jm,jlan)
        ldtm(2*jm-1,jlev,jlas) = fstpm(jlev,1,jm,jlan)- fatpm(jlev,1,jm,jlan)
        ldtm(2*jm  ,jlev,jlas) = fstpm(jlev,2,jm,jlan)- fatpm(jlev,2,jm,jlan)
      END DO

    END DO

    IF (nlevp1>nlev) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jm = 1, nmp1
       lalps  (2*jm-1,jlas) = fstp (nlevp1,1,jm,jlan) - fatp (nlevp1,1,jm,jlan)
       lalps  (2*jm  ,jlas) = fstp (nlevp1,2,jm,jlan) - fatp (nlevp1,2,jm,jlan)
       ldalpsm(2*jm-1,jlas) = fstpm(nlevp1,1,jm,jlan) - fatpm(nlevp1,1,jm,jlan)
       ldalpsm(2*jm  ,jlas) = fstpm(nlevp1,2,jm,jlan) - fatpm(nlevp1,2,jm,jlan)
      END DO
    ENDIF

    ! arrays with spectral coefficients m=0 only

    IF (dc% nlnm0 > 0) then
      lu0 (:,jlas) = fsu0 (:,jlan) - fau0 (:,jlan)
      ldu0(:,jlas) = fsdu0(:,jlan) - fadu0(:,jlan)
      lu0 (:,jlan) = fsu0 (:,jlan) + fau0 (:,jlan)
      ldu0(:,jlan) = fsdu0(:,jlan) + fadu0(:,jlan)
      lul (:,jlan) = lu (1,:,jlan)
      lul (:,jlas) = lu (1,:,jlas)
    ENDIF
  END DO

END SUBROUTINE sym2
