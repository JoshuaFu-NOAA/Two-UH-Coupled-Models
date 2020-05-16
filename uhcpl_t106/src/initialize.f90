SUBROUTINE initialize

  ! Description:
  !
  ! Set up constants in various modules.
  !
  ! Method:
  !
  ! This subroutine initializes all the variables and arrays
  ! in modules.
  !
  ! *initialize* is called from *control*
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G AND M.J, ECMWF, December 1982, changed
  ! U. Schlese, DKRZ, in 1994, and 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_io,           ONLY: IO_init
  USE mo_tracer,       ONLY: initrac
  USE mo_parallel_slt, ONLY: setup_overlap

  IMPLICIT NONE

  !  External subroutines 
  EXTERNAL inidoc, inictl, setdyn, setphys, setrad, suslt


  !  Executable statements 

  !-- 1. Set control variables

  !-- 1.1 Initialize netCDF IO

  CALL IO_init

  !-- 1.2 Preset constants in *mo_doctor*

  CALL inidoc

  !-- 1.3 Preset values for tracer transport

  CALL initrac

  !-- 1.4 Set general control variables 

  CALL inictl

  !-- 2. Compute decomposition 

  CALL init_decomposition

  !-- 3. Preset, modify and derive values needed in the
  !      dynamics and the initialisation and call helmo the first time.

  CALL setdyn

  !-- 4. Preset, modify and derive values needed in the physics.

  CALL setphys

  !-- 5. Preset, modify and derive values needed in the radiation.

  CALL setrad

  !-- 6. Preset values needed in the slt-scheme

  CALL suslt

  !-- 7. Determine overlap information for the slt

  CALL setup_overlap

END SUBROUTINE initialize
