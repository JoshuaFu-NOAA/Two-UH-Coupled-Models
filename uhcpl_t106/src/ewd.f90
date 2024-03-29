!+ computes east west derivatives.
!+ $Id: ewd.f90,v 1.13 1999/08/02 10:39:46 m214089 Exp $

SUBROUTINE ewd

  ! Description:
  !
  ! Computes east west derivatives.
  !
  ! Method:
  !
  ! This subroutine computes east-west derivatives
  ! divided by(a*(1.-mu**2)).
  !
  ! *ewd* is called from *scan1sl*.
  !
  ! The east- west derivation is reduced in *fourier space to
  ! a multiplication by i*m.
  !
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_buffer_fft,    ONLY: falps, ft, fdalpsl, fdtl
  USE mo_gaussgrid,     ONLY: racst
  USE mo_truncation,    ONLY: am
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  REAL    :: zmracst, zracst
  INTEGER :: irow1, jlev, jm, jlat
  INTEGER :: n2mp1, nlev, nlevp1, nlp2, nmp1, nlat
  INTEGER ,pointer :: glat (:)

  !  Executable statements 

!-- array bounds on local PE (Fourier space)

  nlev   =  dc% nflev
  nlevp1 =  dc% nflevp1
  nlp2   =  dc% nlon + 2
  nmp1   =  dc% nm + 1
  n2mp1  = (dc% nm + 1) * 2
  nlat   =  dc% nflat
  glat   => dc% glat

!-- 1. Compute east-west derivatives and blank end of the fields.

  DO jlat = 1, nlat
    irow1 = glat (min(jlat, nlat-jlat+1))

    zracst = racst(irow1)

    DO jlev = 1, nlev

!DIR$ IVDEP
!OCL NOVREC
      DO jm = 1, nmp1
        zmracst = am(jm)*zracst
        fdtl(2*jm-1,jlev,jlat) = -zmracst*ft(2*jm  ,jlev,jlat)
        fdtl(2*jm  ,jlev,jlat) =  zmracst*ft(2*jm-1,jlev,jlat)
      END DO

      fdtl(n2mp1+1:,jlev,jlat) = 0.

    END DO

    IF (nlevp1>nlev) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO jm = 1, nmp1
        zmracst = am(jm)*zracst
        fdalpsl(2*jm-1,jlat) = -zmracst*falps(2*jm  ,jlat)
        fdalpsl(2*jm  ,jlat) =  zmracst*falps(2*jm-1,jlat)
      END DO
      fdalpsl(n2mp1+1:,jlat) = 0.
    END IF

  END DO

END SUBROUTINE ewd
