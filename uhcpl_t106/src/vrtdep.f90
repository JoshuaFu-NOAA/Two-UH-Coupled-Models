!+ compute vertical departure point and departure point index.
!+ $Id: vrtdep.f90,v 1.5 1999/08/27 17:00:38 m214030 Exp $

SUBROUTINE vrtdep(pmap,dt,iterdp,wb,wst,wsb,sig,sigh,dsigh,kdpmpf,kdpmph, &
&      sigmp,sigdp,kdp)

  ! Description:
  !
  ! Compute vertical departure point and departure point index.
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

  USE mo_grid, ONLY: plevp1, plev, plon

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: dt
  INTEGER, INTENT (IN) :: iterdp, pmap

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: dsigh(plevp1), sig(plev), sigh(plevp1), &
&                      wb(plon,plevp1), wsb(plon,plevp1), wst(plon,plevp1)
  INTEGER, INTENT (IN) :: kdpmpf(pmap), kdpmph(pmap)

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: sigmp(plon,plev)

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: sigdp(plon,plev)
  INTEGER, INTENT (OUT) :: kdp(plon,plev)

  ! pmap    Dimension of kdpmap arrays
  ! dt      Time interval that parameterizes the parcel trajectory.
  ! iterdp  Number of iterations used for departure point calculation.
  ! wb      Vertical velocity component (sigma dot).
  ! wst     z-derivs at the top edge of each interval contained in wb
  ! wsb     z-derivs at the bot edge of each interval contained in wb
  ! sig     Sigma values at the full-index levels.
  ! sigh    Half-index sigma levels including sigh(1) = sigma(1/2) = 0.0
  !         sigh(plev+1) = sigma(plev+1/2) = 1.0 .  Note that in general
  !         sigh(k) .lt. sig(k)  where sig(k) is the sigma value at the
  !         k_th full-index level.
  ! dsigh   Increment in half-index sigma levels.
  ! kdpmpf  Array of indices of the model full levels which are mapped
  !         into an artificial evenly spaced vertical grid.  Used to aid
  !         in search for vertical position of departure point
  ! kdpmph  Array of indices of the model half levels which are mapped
  !         into an artificial evenly spaced vertical grid.  Used to aid
  !         in search for vertical position of departure point
  ! sigmp   Sigma value at the trajectory midpoint for each gridpoint
  !         in a vertical slice from the global grid.  On entry sigmp is
  !         an initial guess.
  ! sigdp   Sigma value at the trajectory endpoint for each gridpoint
  !         in a vertical slice from the global grid.
  ! kdp     Vertical index for each gridpoint.  This index points into a
  !         vertical slice array whose vertical grid is given by sig.
  !         E.g.,   sig(kdp(i,k)) .le. sigdp(i,k) .lt. sig(kdp(i,k)+1).

  !  Local scalars: 
  INTEGER :: i, iter, k

  !  Local arrays: 
  REAL :: wmp(plon,plev)

  !  External subroutines 
  EXTERNAL herzin, kdpfnd, vdplim


  !  Executable statements 

  ! Loop over departure point iterates.

  DO iter = 1, iterdp

    ! Compute midpoint indices in half-index sigma-level arrays (use kdp
    ! as temporary storage).

    CALL kdpfnd(plevp1,pmap,sigh,sigmp,kdpmph,kdp)

    ! Interpolate sigma dot field to trajectory midpoints using Hermite
    ! cubic interpolant.

    CALL herzin(plevp1,1,wb,wst,wsb,sigh,dsigh,sigmp,kdp,wmp)

    ! Update estimate of trajectory midpoint.

    DO k = 1, plev
      DO i = 1, plon
        sigmp(i,k) = sig(k) - .5*dt*wmp(i,k)
      END DO
    END DO

    ! Restrict vertical midpoints to be between the top and bottom half-
    ! index sigma levels.

    CALL vdplim(plevp1,sigh,sigmp)
  END DO

  ! Compute trajectory endpoints.

  DO k = 1, plev
    DO i = 1, plon
      sigdp(i,k) = sig(k) - dt*wmp(i,k)
    END DO
  END DO

  ! Restrict vertical departure points to be between the top and bottom
  ! full-index sigma levels.

  CALL vdplim(plev,sig,sigdp)

  ! Vertical indices for trajectory endpoints that point into full-index
  ! sigma level arrays.

  CALL kdpfnd(plev,pmap,sig,sigdp,kdpmpf,kdp)

  RETURN
END SUBROUTINE vrtdep
