!+ computes maximum winds
!+ $Id: maxwind.f90,v 1.12 1999/10/14 07:29:56 m214003 Exp $

SUBROUTINE maxwind(ulz, vmaxz)

  ! Description:
  !
  ! Computes maximum winds for horizontal diffusion
  ! and diagnostics.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters
  USE mo_control
  USE mo_gaussgrid
  USE mo_semi_impl,     ONLY: nulev, uvmax, vcheck, vmax 
  USE mo_start_dataset
  USE mo_constants
  USE mo_doctor,        ONLY: nout
  USE mo_decomposition, ONLY: dc => local_decomposition
  USE mo_global_op,     ONLY: maxval_latit

  IMPLICIT NONE

  !  Array arguments 
  REAL :: ulz  (nlev, dc%nglat)
  REAL :: vmaxz(nlev, dc%nglat)

  !  Local array bounds
  INTEGER :: nglon, nglpx, nglat

  ! Local arrays
  INTEGER :: nmaxlev(1)   

  !  Local scalars: 
!  INTEGER :: jlat, klat

  !  Executable statements 

  !  Local array bounds
  nglon = dc%nglon ! number of longitudes
  nglpx = dc%nglpx ! number of longitudes allocated
  nglat = dc%nglat ! number of latitudes

!-- 1. Maximum level winds for diffusion

  vmax(:) = maxval_latit (vmaxz(:,:))

!-- 2. Maximum wind for diagnostics

  nmaxlev = MAXLOC(vmax)                    ! vertical index of max
  nulev   = nmaxlev(1)                      ! for postatd
  uvmax   = vmax(nulev)                     ! max value (absolute)
!  jlat = ilmax(nulev)                      ! latitude index of max
!  klat = min(2*jlat-1,2*(ngl+1-jlat))      ! ping pong index
!  ulat = ASIN(twomu(klat)*0.5)*180./api    ! for postatd
!  ulm = ulz(nulev,jlat)/sqcst(klat)        ! for postatd

  ! Check for high windspeeds

  IF (uvmax > vcheck) THEN
    WRITE (nout,'(a,f6.0,a)') ' WARNING! high wind speed: ', uvmax, ' m/s'
!    WRITE (nout,*) ' Level: ', nulev, ' Latitude: ', jlat, ' NSTEP= ', nstep
    WRITE (nout,*) ' Level: ', nulev, ' NSTEP= ', nstep
  END IF

END SUBROUTINE maxwind
