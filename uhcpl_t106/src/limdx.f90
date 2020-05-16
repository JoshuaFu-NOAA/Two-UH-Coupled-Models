!+ limit the derivative estimates for data
!+ $Id: limdx.f90,v 1.6 1999/08/27 17:00:27 m214030 Exp $

SUBROUTINE limdx(pidim,ibeg,len,dx,f,fxl,fxr)

  ! Description:
  !
  ! Limit the derivative estimates for data.
  !
  ! Method:
  !
  ! Limit the derivative estimates for data on an equally spaced grid
  ! so they satisfy the SCM0 condition, that is, the spline will be
  ! monotonic, but only C0 continuous on the domain
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception, ONLY: finish
  USE mo_doctor,    ONLY: nout
  USE mo_grid,      ONLY: pbpts

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: dx
  INTEGER, INTENT (IN) :: ibeg, len, pidim

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: f(pidim)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: fxl(pidim), fxr(pidim)

  ! pidim   Length of f, fxl, and fxr.
  ! ibeg    First interval of grid for which derivatives are computed.
  ! len     Number of grid intervals for which derivatives are computed.
  !         (There are pidim - 1 intervals between the pidim gridpoints
  !         represented in f, fxl, and fxr.)
  ! dx      Value of grid spacing.
  ! f       Values on equally spaced grid from which derivatives fxl and
  !         fxr were computed.
  ! fxl     fxl(i) is the limited derivative at the left  edge of interval
  ! fxr     fxr(i) is the limited derivative at the right edge of interval

  !  Local scalars: 
  REAL :: rdx                ! 1./dx
  INTEGER :: i               ! index
  INTEGER :: iend            ! index to end work on vector

  !  Local arrays: 
  REAL :: deli(pbpts)        ! simple linear derivative

  !  External Subroutines 
  EXTERNAL scm0


  !  Executable Statements 

  IF (pidim>pbpts) THEN
    WRITE (nout,'(A)') 'LIMDX: Local work array DELI not dimensioned large enough'
    WRITE (nout,'(A,I5)') '       Increase local parameter pbpts to ',pidim
    CALL finish('limdx','Run terminated.')
  END IF

  iend = ibeg + len - 1
  rdx = 1./dx

  DO i = ibeg, iend
    deli(i) = (f(i+1)-f(i))*rdx
  END DO

  ! Limiter

  CALL scm0(len,deli(ibeg),fxl(ibeg),fxr(ibeg))

  RETURN
END SUBROUTINE limdx
