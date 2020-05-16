!+ interpolate field on vertical slice to vertical departure point
!+ $Id: herzin.f90,v 1.8 1999/08/27 17:00:23 m214030 Exp $

SUBROUTINE herzin(pkdim,pf,f,fst,fsb,sig,dsig,sigdp,kdp,fdp)

  ! Description:
  !
  ! Interpolate field on vertical slice to vertical departure point
  !
  ! Method:
  !
  ! Interpolate field on vertical slice to vertical departure point using
  ! Hermite cubic interpolation.
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

  USE mo_grid, ONLY: plon, plev

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pf, pkdim

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: dsig(pkdim), f(plon,pkdim,pf), fsb(plon,pkdim,pf),  &
&                      fst(plon,pkdim,pf), sig(pkdim), sigdp(plon,plev)
  INTEGER, INTENT (IN) :: kdp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: fdp(plon,plev,pf)

  ! pkdim   Vertical dimension of vertical slice arrays.
  ! pf      Number of fields being interpolated.
  ! f       Vertical slice of data to be interpolated.
  ! fst     z-derivatives at the top edge of each interval contained in f
  ! fsb     z-derivatives at the bot edge of each interval contained in f
  ! sig     Sigma values corresponding to the vertical grid
  ! dsig    Increment in sigma value for each interval in vertical grid.
  ! sigdp   Sigma value at the trajectory midpoint or endpoint for each
  !         gridpoint in a vertical slice from the global grid.
  ! kdp     Vertical index for each gridpoint.  This index points into a
  !         vertical slice array whose vertical grid is given by sig.
  ! E.g.,   sig(kdp(i,j)) .le. sigdp(i,j) .lt. sig(kdp(i,j)+1) .
  ! fdp     Value of field at the trajectory midpoints or endpoints.

  !  Local scalars: 
  ! vert interval containing the dep. pt.
  REAL :: dzk, zb, zt
  INTEGER :: i, k, m

  !  Local arrays: 
  ! interpolation coefficients
  REAL :: dhb(plon,plev), dht(plon,plev), hb(plon,plev), ht(plon,plev)

!DIR$ NOBOUNDS

  !  Executable statements 

  DO k = 1, plev
    DO i = 1, plon
      dzk      = dsig(kdp(i,k))
      zt       = (sig(kdp(i,k)+1)-sigdp(i,k))/dzk
      zb       = 1. - zt
      ht(i,k)  = (3.0-2.0*zt)*zt**2
      hb(i,k)  = (3.0-2.0*zb)*zb**2
      dht(i,k) = -dzk*(zt-1.)*zt**2
      dhb(i,k) =  dzk*(zb-1.)*zb**2
    END DO
  END DO

  ! Loop over fields.

  DO m = 1, pf
    DO k = 1, plev
      DO i = 1, plon
        fdp(i,k,m) = f(i,kdp(i,k)  ,m)*ht(i,k) + fst(i,kdp(i,k),m)*dht(i,k) + &
&                    f(i,kdp(i,k)+1,m)*hb(i,k) + fsb(i,kdp(i,k),m)*dhb(i,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE herzin
