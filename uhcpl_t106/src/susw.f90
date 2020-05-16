!+ initialize module shortwave
!+ $Id: susw.f90,v 1.5 1999/04/26 07:47:38 m214003 Exp $

SUBROUTINE susw

  ! Description:
  !
  ! Initialize shortwave, the module that contains coefficients
  ! needed to run the shortwave radiation subroutines
  !
  ! Method:
  !
  ! Reference:
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

  USE mo_constants
  USE mo_shortwave

  IMPLICIT NONE

  !  Local parameters: 
  REAL, PARAMETER :: zpdh2o = 0.8, zpdumg = 0.75, zprh2o = 30000., zprumg = 30000.

  !  Local scalars: 
  REAL :: zh2o, zumg


  !  Executable statements 

  rpdh1 = zpdh2o + 1.
  rpdu1 = zpdumg + 1.
  zh2o = 1./(10.*g*rpdh1)
  zumg = 1./(10.*g*rpdu1)
  rpnu = zumg/(zprumg**zpdumg)
  rpnh = zh2o/(zprh2o**zpdh2o)

  rswcp = 0.002*rswce

  RETURN
END SUBROUTINE susw
