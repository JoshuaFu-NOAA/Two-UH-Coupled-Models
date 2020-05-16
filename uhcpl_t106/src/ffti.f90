!+ inverse Fourier transforms.
!+ $Id: ffti.f90,v 1.24 1999/11/03 07:47:34 m214030 Exp $

SUBROUTINE ffti

  ! Description:
  !
  ! Inverse Fourier transforms.
  !
  ! Method:
  !
  ! This subroutine calls the fast *Fourier transforms
  ! subroutine (*fft991cy*) in order to perform the inverse *Fourier
  ! transforms.
  !
  ! *ffti* is called from *scan1sl*.
  ! The parameters needed for the *Fourier transforms are obtained
  ! through module *mo_fft* (they were initialized in *inifft*).
  !
  ! *fft991cy*    fast *Fourier transforms.
  !
  ! 1-appendix *b1:organisation of the spectral model.
  ! m.j       5/10/81.
  ! 2-description of fft991cy
  ! c.t       ????????
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_buffer_fft,    ONLY: fftz, nvar
  USE mo_fft,           ONLY: nfax, trig, fft991cy
  USE mo_decomposition, ONLY: dc => local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: inc, isign
  INTEGER :: nlon, nlp2, nlev, nlat
  LOGICAL :: col_1d
!  INTEGER ::  nvar

  !  Local arrays: 
!  REAL :: zffti((7*nlev+3)*nlp2), zwork((7*nlev+3)*nlp2)
! f95
!  REAL :: zwork(size(fftz))
  REAL :: zwork((dc%nlon+2) * dc%nflevp1 * dc%nflat * nvar)

!-- 2. Inverse *Fourier transforms

!-- 2.1 Set constants

  inc    = 1
  isign  = 1
  nlon   = dc% nlon
  nlp2   = nlon + 2
!  nvar  = size(fftz,4)
  nlev   = size(fftz,2)
  nlat   = dc% nflat
  col_1d = dc% col_1d

!-- 2.2 fft(*vo*,*d*,*t*,*alps*,*u*,*v*,*dtl*,*dtm*,*dalpsl*,*dalpsm*)

  IF (.not.col_1d) THEN
    CALL fft991cy(fftz,zwork,trig,nfax,inc,nlp2,nlon,nvar*nlev*nlat,isign)
  ENDIF

END SUBROUTINE ffti
