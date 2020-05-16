!+ controls the initial start-up of the experiment
!+ $Id: start.f90,v 1.11 1999/09/16 08:07:37 k202057 Exp $

SUBROUTINE start

  ! Description:
  !
  ! Controls the setting up of work files and the initial
  ! start-up of the experiment.
  !
  ! Externals:
  ! *scan2*     - scan 2 through latitudes to perform
  !               inverse *Legendre transforms
  ! *helmo*     - computes matrix for inversion of helmholtz equat
  ! *inhysi*    - computes constants for vertical
  !               part of semi-implicit time scheme
  ! *inisu0*    - computes initial spectral components
  !               of the zonal mean wind
  ! *ioinitial* - reads initial fields and sets up work files
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, November 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,       ONLY: nresum
  USE mo_start_dataset, ONLY: nstart, nstep
  USE mo_doctor,        ONLY: nout
  USE mo_mpi,           ONLY: p_pe, p_io  

  IMPLICIT NONE

  !  External subroutines 
  EXTERNAL inhysi, ioinitial, scan2


  !  Executable statements 

!-- 1. Read initial fields

!-- 1.1 Read and position files and set up work files

  CALL ioinitial

  IF (p_pe == p_io) THEN
     WRITE (nout,'(a,/)') ' File setup done and spectral arrays allocated.'
  END IF

!-- 2. Update constants for semi implicit scheme

!-- 2.2 Initialise constants for the vertical part of the
!       semi-implicit scheme

  nstep = 0
  CALL inhysi

!-- 2.4 Set forecast step to -1 to distinguish as start up
!       and convert spectral coefficients to fourier coefficients

  nstep = -1
  CALL scan2

!-- 3. Reset switches

  nstep = 0
  nresum = 0
  nstart = 0

END SUBROUTINE start
