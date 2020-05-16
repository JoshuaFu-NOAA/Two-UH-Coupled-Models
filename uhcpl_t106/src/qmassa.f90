!+ calculate contribution of current latitude to mass of constituents
!  being advected by slt.
!+ $Id: qmassa.f90,v 1.8 1999/09/08 16:38:41 m214030 Exp $

SUBROUTINE qmassa(cwava,w,q3,pdel,hw1lat,jlat)

  ! Description:
  !
  ! Calculate contribution of current latitude to mass of constituents
  ! being advected by slt.
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          P. Rasch, D. Williamson, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! parallel version:  T. Diehl, DKRZ, July 1999
  !                    A. Rhodin, MPI, Sept 1999
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_grid,          ONLY: plond, plev, pcnst, plon, i1
  USE mo_global_op,     ONLY: sum_zonal_sl

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  REAL, INTENT (IN) :: cwava                ! normalization factor l/(g*plon)
  REAL, INTENT (IN) :: w                    ! gaussian weight this latitude
  INTEGER, INTENT (IN) :: jlat ! local row index N->S

  !  Array arguments with intent(In):
  REAL, INTENT (IN) :: pdel(plond,plev)     ! pressure diff between interfaces
  REAL, INTENT (IN) :: q3(plond,plev,pcnst) ! constituents

  !  Array arguments with intent(Out):
  REAL, INTENT (OUT) :: hw1lat(pcnst)       ! accumulator

  !  Local scalars: 
  INTEGER :: m

  !  Executable statements 

  ! longitude, level, constituent indices

  DO m = 1, pcnst
    hw1lat(m) = sum_zonal_sl(&
                  q3(i1:plon+i1-1,:,m)*pdel(i1:plon+i1-1,:),jlat=jlat)
!   more efficient but not identical to serial version:
!   hw1lat(m) = sum_zonal (                                             &
!     sum ( q3(i1:plon+i1-1,:,m)*pdel(i1:plon+i1-1,:), dim=2), jlat=jlat)
  END DO

     ! The 0.5 factor arises because gaussian weights sum to 2

  DO m = 1, pcnst
     hw1lat(m) = cwava*w*hw1lat(m)*0.5
  END DO

END SUBROUTINE qmassa
