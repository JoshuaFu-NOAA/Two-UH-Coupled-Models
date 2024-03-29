MODULE mo_truncation

  USE mo_parameters, ONLY: jphgl, jpmp1, jpnlev
  USE mo_doctor,     ONLY: nout

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmp, nnp, am, ntrn, ntrm, ntrk
  PUBLIC :: scpar               ! initialization routine

  ! ---------------------------------------------------------------
  !
  ! module *mo_truncation* - quantities related to the spectral truncation.
  !
  ! ---------------------------------------------------------------

  INTEGER, ALLOCATABLE :: ntrm(:)  !  max zonal wave number.
  INTEGER, ALLOCATABLE :: ntrk(:)  !  max meridional wave number.
  INTEGER :: ntrn(jpnlev)  !  max meridional wave number for m=0.
  INTEGER :: nmp(jpmp1)    !  displacement of the first point of columns
                           !  for computations in spectral space.
  INTEGER :: nnp(jpmp1)    !  number of points on each column.
  REAL    :: am(jpmp1)     !  float(m).i.e jm-1 in the jm loops.

CONTAINS

  SUBROUTINE scpar (nm, nn, nk)

    USE mo_mpi, ONLY: p_pe, p_io

    INTEGER ,intent(in) :: nm ! max zonal wave number
    INTEGER ,intent(in) :: nn ! max meridional wave number for m=0
    INTEGER ,intent(in) :: nk ! max meridional wave number

    ! Description:
    !
    ! Computes parameters used for computations in spectral space.
    !
    ! Method:
    !
    ! This subroutine computes some parameters related to the
    ! truncation and used for the computations in spectral space and
    ! for the *Legendre transforms.
    !
    ! *scpar* is called from *initialise*
    !
    ! Results:
    ! The results are stored in arrays in module *mo_truncation*
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, March 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

!   USE mo_control,    only: nkp1, nmp1, nn
!   USE mo_truncation, only: am, nmp, nnp

    !  Local scalars: 
    INTEGER :: jm

    !  Intrinsic functions 
    INTRINSIC MIN

    !  Executable statements 

!-- 0. These parameters may be overwritten later in 'setdyn'

    ntrm (:) = nm
    ntrn (:) = nn
    ntrk (:) = nk

!-- 1. Preliminary computations

!-- 2. Compute parameters

    DO jm = 1, nm+1
      nnp(jm) = MIN(nk+1-jm,nn) + 1
    END DO

    nmp(1) = 0
    DO jm = 2, nm+1
      nmp(jm) = nmp(jm-1) + nnp(jm-1)
    END DO

!-- 3. Fill arrays *am* and *annp1*

    DO jm = 1, nm+1
      am(jm) = jm - 1.
    END DO

    IF (p_Pe == p_io) THEN
       WRITE (nout, '(a)') &
            ' Number of points on each column (NNP): '  
       WRITE (nout, '(11i7)') nnp(1:nm+1)
       WRITE (nout, '(a)') &
            ' Displacement of the first point of columns (NMP): '  
       WRITE (nout, '(11i7)') nmp(1:nm+1)
    END IF

  END SUBROUTINE scpar

END MODULE mo_truncation
