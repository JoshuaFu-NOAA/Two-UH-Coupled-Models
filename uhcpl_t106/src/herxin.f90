!+ interpolate to its x value at each latitude
!+ $Id: herxin.f90,v 1.10 2000/03/20 14:17:18 m214003 Exp $

SUBROUTINE herxin(pf,fb,fxl,fxr,x,xdp,idp,jdp,fint)

  ! Description:
  !
  ! Interpolate to its x value at each latitude
  !
  ! Method:
  !  
  ! For each departure point in the latitude slice being forecasted,
  ! interpolate (using equally spaced Hermite cubic formulas) to its
  ! x value at each latitude required for later interpolation in the y
  ! direction.
  !
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! Modified:          U. Schlese,  DKRZ - Hamburg, May 1994
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_grid,          ONLY: plond, plev, platd, plon, nxpt, istart
  USE mo_slt,           ONLY: ppdy
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pf

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: fb(plond,plev,pf,platd), &
                       fxl(plond,plev,pf,platd), fxr(plond,plev,pf,platd), &
                       x(plond), xdp(plon,plev)

  INTEGER, INTENT (IN) :: idp(plon,plev), jdp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: fint(plon,plev,ppdy,pf)

  ! pf      Number of fields being interpolated.
  ! fb      extended array of data to be interpolated.
  ! fxl     x derivatives at the left edge of each interval containing
  !         the departure point
  ! fxr     x derivatives at the right edge of each interval containing
  !         the departure point
  ! x       Equally spaced x grid values in extended arrays.
  ! xdp     xdp(i,k) is the x-coordinate (extended grid) of the
  !         departure point that corresponds to global grid point (i,k)
  !         in the latitude slice being forecasted.
  ! idp     idp(i,k) is the index of the x-interval (extended grid) that
  !         contains the departure point corresponding to global grid
  !         point (i,k) in the latitude slice being forecasted.
  !         Note that x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
  ! jdp     jdp(i,k) is the index of the y-interval (extended grid) that
  !         contains the departure point corresponding to global grid
  !         point (i,k) in the latitude slice being forecasted.
  !         Suppose yb contains the y-coordinates of the extended array
  !         and ydp(i,k) is the y-coordinate of the departure point
  !         corresponding to grid point (i,k).  Then,
  !         yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
  ! fint    (fint(i,k,j,n),j=1,ppdy) contains the x interpolants at each
  !         latitude needed for the y derivative estimates at the
  !         endpoints of the interval that contains the departure point
  !         for grid point (i,k).  The last index of fint allows for
  !         interpolation of multiple fields.

  !  Local scalars: 
  REAL :: dx        ! x-increment
  REAL :: rdx       ! 1./dx
  ! interpolation coefficient
  REAL :: xl, xr
  REAL :: pi, xx
  INTEGER :: i, k, m

  !  Local arrays: 
  ! interpolation coefficient
  REAL    :: dhl(plon,plev), dhr(plon,plev), hl(plon,plev), hr(plon,plev)
  INTEGER :: id2(plon,plev), ioff


  !  Executable statements 
! dx = x(nxpt+2) - x(nxpt+1)
  ! parallel version
  pi = 4.0*ATAN(1.0)
  dx = 2.0*pi/REAL(dc%nlon)
  ioff = dc%glons(1) - istart

  rdx = 1./dx
  DO k = 1, plev
    DO i = 1, plon
      xx = REAL(idp(i,k)+ioff)*dx         ! same as in grdxy
      xl = ( xx           -xdp(i,k))*rdx
!     xl = ( x(idp(i,k)+1)-xdp(i,k))*rdx  ! old version, relies on extended grid
      xr = 1. - xl
      hl(i,k) = (3.0-2.0*xl)*xl**2
      hr(i,k) = (3.0-2.0*xr)*xr**2
      dhl(i,k) = -dx*(xl-1.)*xl**2
      dhr(i,k) = dx*(xr-1.)*xr**2
    END DO
  END DO

  ! access inner domain only, if 1 PE in e-W direction

  id2 = idp
  IF (dc% nprocb == 1) THEN
    WHERE (id2 >  nxpt + plon) id2 = id2 - plon
    WHERE (id2 <= nxpt       ) id2 = id2 + plon
  ENDIF

  ! x interpolation at each latitude needed for y interpolation.
  ! Once for each field.

!vdir noloopchg
  DO m = 1, pf
    DO k = 1, plev
      DO i = 1, plon
        fint(i,k,1,m) = fb(id2(i,k),k,m,jdp(i,k)-1)*hl(i,k) +   &
                        fb(id2(i,k)+1,k,m,jdp(i,k)-1)*hr(i,k) + &
                        fxl(id2(i,k),k,m,jdp(i,k)-1)*dhl(i,k) + &
                        fxr(id2(i,k),k,m,jdp(i,k)-1)*dhr(i,k)
        fint(i,k,2,m) = fb(id2(i,k),k,m,jdp(i,k))*hl(i,k) +     &
                        fb(id2(i,k)+1,k,m,jdp(i,k))*hr(i,k) +   &
                        fxl(id2(i,k),k,m,jdp(i,k))* &
                        dhl(i,k) + fxr(id2(i,k),k,m,jdp(i,k))*dhr(i,k)
        fint(i,k,3,m) = fb(id2(i,k),k,m,jdp(i,k)+1)*hl(i,k) +   &
                        fb(id2(i,k)+1,k,m,jdp(i,k)+1)*hr(i,k) + &
                        fxl(id2(i,k),k,m,jdp(i,k)+1)*dhl(i,k) + &
                        fxr(id2(i,k),k,m,jdp(i,k)+1)*dhr(i,k)
        fint(i,k,4,m) = fb(id2(i,k),k,m,jdp(i,k)+2)*hl(i,k) +   &
                        fb(id2(i,k)+1,k,m,jdp(i,k)+2)*hr(i,k) + &
                        fxl(id2(i,k),k,m,jdp(i,k)+2)*dhl(i,k) + &
                        fxr(id2(i,k),k,m,jdp(i,k)+2)*dhr(i,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE herxin
