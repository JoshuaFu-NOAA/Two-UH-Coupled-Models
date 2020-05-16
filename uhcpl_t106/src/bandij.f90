!+ Calculate longitude and latitude indices
!+ $Id: bandij.f90,v 1.26 2000/03/20 14:17:17 m214003 Exp $

SUBROUTINE bandij(dlam,phib,lamp,phip,iband,jband,nxpt_a)

  ! Description:
  !
  ! Calculate longitude and latitude indices
  !
  ! Method:
  !
  ! Calculate longitude and latitude indices that identify the
  ! intervals in the extended grid that contain the departure points.
  ! Upon entry, all dep. points should be within jintmx intervals of the
  ! Northern- and Southern-most model latitudes. Note: the algorithm
  ! relies on certain relationships of the intervals in the Gaussian grid.
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! parallel version:  T. Diehl, July 1999
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_grid,          ONLY: pgls, platd, plev, istart, plon, nxpt, dphibr
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_doctor,        ONLY: nerr
  USE mo_exception,     ONLY: finish

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL,    INTENT (IN) :: dlam

  !  Array arguments with intent(In):
  REAL,    INTENT (IN) :: lamp(pgls), phib(platd), phip(pgls)
  INTEGER, INTENT (IN) :: nxpt_a(plev)  

  !  Array arguments with intent(Out):
  INTEGER, INTENT (OUT) :: iband(pgls), jband(pgls)

  ! dlam    Length of increment in equally spaced longitude grid (radians)
  ! phib    Latitude values for the extended grid.
  ! lamp    Longitude coordinates of the points.  It is assumed that
  !         0.0 .le. lamp(i) .lt. 2*pi .
  ! phip    Latitude coordinates of the points.
  ! iband   Longitude index of the points.  This index points into
  !         the extended arrays, e.g.,
  !         lam(iband(i)) .le. lamp(i) .lt. lam(iband(i)+1) .
  ! jband   Latitude index of the points.  This index points into
  !         the extended arrays, e.g.,
  !         phib(jband(i)) .le. phip(i) .lt. phib(jband(i)+1) .

  !  Local scalars: 
  INTEGER :: i           ! index
  INTEGER :: k, ieloc(1)
  LOGICAL :: lerror(pgls)

  !  Intrinsic functions 
  INTRINSIC INT, MERGE, ANY, MAXLOC


  !  Executable statements 

  ! Longitude indices.

  DO i = 1, pgls

     !   iband(i) = istart + CEILING(lamp(i)/dlam) - dc%glons(1)
     iband(i) = istart + INT(lamp(i)/dlam) - dc%glons(1) &
              + MERGE(0,1,lamp(i)<0.0)

     ! For dynamic sltpads need to check the level maximum

     k = (i-1)/plon + 1  ! determine current level

     !   lerror(i) = (iband(i) <  istart-nxpt_a(k)) .OR. &
     !               (iband(i) >= plon+1+nxpt+nxpt_a(k))
     !
     ! the original check was not sufficient
     !   cf. subroutine cubxdr:
     !   arrays are indexed with i-1 to i+2
     !
     !   if only one PE is present in e-W direction the test is skipped.
     !   nxpt must be at least 2 in this case
     !
     lerror(i) = .false.
     IF (dc% nprocb > 1) THEN
        lerror(i) = (iband(i) <  istart   -nxpt_a(k) +1) .OR. &
                    (iband(i) >  plon+nxpt+nxpt_a(k) -2)
     ENDIF

  END DO

  IF (ANY(lerror)) THEN
     ieloc = MAXLOC((/(i,i=1,pgls)/),lerror)
     WRITE (nerr,*) 'bandij: set B = ', dc%set_b, &
                    ' set A = ', dc%set_a, &
                    ' bandi = ', iband(ieloc(1))
     WRITE (nerr,*) '       Increase nxpt (in mo_grid) and/or isave (in sltb1)'
     CALL finish('bandij','Point out of bounds in E-W.')
  END IF

  ! Latitude indices.

  jband(:) = INT ((phip(:) - phib(1))*dphibr + 1.)

#ifndef NAG
  WHERE(phip(:) >= phib(jband(:)+1)) jband(:) = jband(:) + 1
#else
  DO i = 1, pgls
     IF (phip(i) >= phib(jband(i)+1)) jband(i) = jband(i) + 1
  END DO
#endif

  !
  !  lerror(:) = jband(:) < 1 .OR. jband(:) > platd
  !
  ! the original check was not sufficient
  !   cf. subroutine herxin:
  !   arrays are indexed with jband-1 to jband+2
  !
  lerror(:) = jband(:) < 2 .OR. jband(:) > platd -2

  IF (ANY(lerror)) THEN
     ieloc = MAXLOC((/(i,i=1,pgls)/),lerror)
     WRITE (nerr,*) 'bandij: set B = ', dc%set_b, &
                    ' set A = ', dc%set_a, &
                    ' bandij = ',jband(ieloc(1))
     CALL finish('bandij','Point out of bounds in N-S.')
  END IF

  RETURN
END SUBROUTINE bandij

