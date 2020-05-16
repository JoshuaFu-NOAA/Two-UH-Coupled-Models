!+ copy data to the longitude extensions of the extended array
!+ $Id: extx.f90,v 1.11 2000/01/13 11:50:12 m214003 Exp $

!OCL NOVREC

SUBROUTINE extx(pkcnst,pkdim,fb)

  ! Description:
  !
  ! Copy data to the longitude extensions of the extended array
  !
  ! Method:
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

  USE mo_grid,          ONLY: plond, platd, nxpt, plon, i2pi
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pkcnst ! Dimension construct for 3-D arrays
  INTEGER, INTENT (IN) :: pkdim  ! Vertical dimension

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: fb(plond,pkdim*pkcnst,platd) ! Constituents array

  !  Local scalars: 
  INTEGER :: i, j, k, m


  !  Executable statements 

  IF (dc%nprocb > 1) RETURN

  DO m=1,pkcnst

     ! Fill west edge points.
     
     IF (nxpt >= 1) THEN
        DO j = 1, platd
           DO k = 1, pkdim
              DO i = 1, nxpt                 
                 fb(i,k+(m-1)*pkdim,j) = fb(i+plon,k+(m-1)*pkdim,j)
              END DO
           END DO
        END DO
     END IF
     
     ! Fill east edge points
     
     DO j = 1, platd
        DO k = 1, pkdim
           DO i = i2pi, plond
              fb(i,k+(m-1)*pkdim,j) = fb(i-plon,k+(m-1)*pkdim,j)
           END DO
        END DO
     END DO
     
  END DO

  RETURN
END SUBROUTINE extx
