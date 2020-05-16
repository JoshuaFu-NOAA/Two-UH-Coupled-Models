!+ interpolate 2-d field to departure point
!+ $Id: hrvint.f90,v 1.6 1999/09/27 13:25:53 k202057 Exp $

SUBROUTINE hrvint(pf,fb,fxl,fxr,x,wdy,xdp,ydp,idp,jdp,fint,fdp)

  ! Description:
  !
  ! Interpolate 2-d field to departure point
  !
  ! Method:
  !
  ! Interpolate 2-d field to departure point using equivalent of tensor
  ! product Lagrange cubic interpolation. For economy, code actually
  ! uses Hermite cubic with Lagrange cubic derivative estimate in x
  ! and Lagrange cubic in y.
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

  USE mo_grid, ONLY: plond, plev, platd, plon
  USE mo_slt,  ONLY: ppdy

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pf

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: fb(plond,plev,pf,platd), fxl(plond,plev,pf,platd), &
&         fxr(plond,plev,pf,platd), wdy(4,2,platd), x(plond), &
&         xdp(plon,plev),ydp(plon,plev)
  INTEGER :: idp(plon,plev), jdp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: fint(plon,plev,ppdy,pf), fdp(plon,plev,pf)

  ! pf      Number of fields being interpolated.
  ! fb      Extended array of data to be interpolated.
  ! fxl     x-derivatives at the left  edge of each interval containing
  !         the departure point.
  ! fxr     x-derivatives at the right edge of each interval containing
  !         the departure point.
  ! x       Equally spaced x grid values in extended arrays.
  ! wdy     Weights for Lagrange cubic interpolation on the
  !         unequally spaced y-grid.
  ! xdp     xdp(i,k) is the x-coordinate of the departure point that
  !         corresponds to global grid point (i,k) in the latitude slice
  !         being forecasted.
  ! ydp     ydp(i,k) is the y-coordinate of the departure point that
  !         corresponds to global grid point (i,k) in the latitude slice
  !         being forecasted.
  ! idp     idp(i,k) is the index of the x-interval that contains the
  !         departure point corresponding to global grid point (i,k) in
  !         the latitude slice being forecasted.
  !         Note that
  !         x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
  ! jdp     jdp(i,k) is the index of the y-interval that contains the
  !         departure point corresponding to global grid point (i,k) in
  !         the latitude slice being forecasted.
  !         Suppose yb contains the y-coordinates of the extended array
  !         and ydp(i,k) is the y-coordinate of the departure point
  !         corresponding to grid point (i,k).  Then,
  !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
  ! fint    WORK ARRAY, results not used on return
  ! fdp     Value of field at the departure points.

  !  External subroutines 
  EXTERNAL herxin, lagyin


  !  Executable statements 

  ! Hermite cubic interpolation to the x-coordinate of each
  ! departure point at each y-coordinate required to compute the
  ! y-interpolants.

  CALL herxin(pf,fb,fxl,fxr,x,xdp,idp,jdp,fint)

  ! Lagrange cubic interpolation in y

  CALL lagyin(pf,fint,wdy,ydp,jdp,fdp)

  RETURN
END SUBROUTINE hrvint

