!+ preset constants in physical modules.
!+ $Id: setphys.f90,v 1.13 1999/08/17 12:58:03 m214003 Exp $

SUBROUTINE setphys

  ! Description:
  !
  ! Preset constants in physical modules.
  !
  ! Method:
  !
  ! This subroutine preset, modifies and derive
  ! constants in moduls used in the parameterisation
  ! of diabatic processes except the constants internal
  ! to the radiation "black box".
  !
  ! *setphys* is called from *inimod*.
  !
  ! Externals:
  ! *iniphy*   to set constants in *mo_physc2*
  ! *inidia*   to set constants in *mo_diagnostics*
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters
  USE mo_mpi
  USE mo_doctor
  USE mo_control
  USE mo_start_dataset
  USE mo_constants
  USE mo_physc1
  USE mo_param_switches
  USE mo_cumulus_flux
  USE mo_convect_tables, ONLY: set_lookup_tables
  USE mo_exception
  USE mo_diagnostics,    ONLY: inidia

  IMPLICIT NONE

  !  Local scalars:
  REAL    :: zl 
  INTEGER :: ilrad, jlon

  !  Intrinsic functions 

  INTRINSIC COS, SIN

  INCLUDE 'physctl.inc'

  !  Executable statements 

!-- 1. Preset constants

!-- 1.1 Preset constants in *comph1*

  ilrad = 128
  DO jlon = 1, ilrad
    zl = 2.*api*(jlon-1.)/ilrad
    sinrad(jlon) = SIN(zl)
    cosrad(jlon) = COS(zl)
  END DO

!-- 1.2 Preset constants in *compsw*

  lphys = .TRUE.
  lvdiff = .TRUE.
  lconv = .TRUE.
  lcond = .TRUE.
  lsurf = .TRUE.
  lmfpen = .TRUE.
  lmfscv = .TRUE.
  lmfmid = .TRUE.
  lmfdd = .TRUE.
  lmfdudv = .TRUE.
  IF (nn==21) THEN
    lgwdrag = .FALSE.
  ELSE
    lgwdrag = .TRUE.
  END IF
  lice = .TRUE.

  ndiapfr = -240

!-- 1.3 Preset constants in *mo_diagnostics*

  IF ( .NOT. lres) CALL inidia

!-- 1.4 Initialise lookup tables for CUADJTQ

  CALL set_lookup_tables

!-- 2. Read namelist
  
  IF (p_parallel) THEN
     IF (p_parallel_io) THEN
        READ (nin,physctl)
     ENDIF
     CALL p_bcast (lphys, p_io)
     CALL p_bcast (lvdiff, p_io)
     CALL p_bcast (lcond, p_io)
     CALL p_bcast (lsurf, p_io)
     CALL p_bcast (lconv, p_io)
     CALL p_bcast (lmfpen, p_io)
     CALL p_bcast (lgwdrag, p_io)
     CALL p_bcast (lice, p_io)
     CALL p_bcast (ndiapfr, p_io)
  ELSE     
     READ (nin,physctl)
  ENDIF

!-- 3. Modify constants

!-- 3.1 Modify constants in *compsw*

  IF ( .NOT. lphys) THEN
    lvdiff = .FALSE.
    lgwdrag = .FALSE.
    lconv = .FALSE.
    lcond = .FALSE.
    lsurf = .FALSE.
    lice = .FALSE.
  END IF

  IF ( .NOT. lconv) THEN
    lmfpen = .FALSE.
    lmfscv = .FALSE.
    lmfmid = .FALSE.
    lmfdd = .FALSE.
    lmfdudv = .FALSE.
  END IF

  IF ( .NOT. lmfpen) THEN
    lconv = .FALSE.
    lmfscv = .FALSE.
    lmfmid = .FALSE.
    lmfdd = .FALSE.
    lmfdudv = .FALSE.
  END IF

  IF (ndiapfr<0) ndiapfr = -ndiapfr*3600./dtime + 0.1

!-- 3.2 Set up constants in *comph2*

  CALL iniphy

!-- 4. Write namelist

#ifdef DEBUG
  WRITE (nerr,physctl)
#endif

  RETURN
END SUBROUTINE setphys
