SUBROUTINE control

  ! Description:
  !
  ! Control routine for new model.
  !
  ! Method:
  !
  ! This subroutine controls the running of the model.
  ! 
  ! *control* is called from the main program (*master*).
  !
  ! Externals:
  !   *iniloc*     called to initialize interface routine tables.
  !   *initialize* called to initialize modules.
  !   *inileg*     called to compute *Legendre polynomials and
  !                the parameters needed for the *Legendre transforms.
  !   *start*      supervises a newstart
  !   *restart*    supervise a restart.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G and M.J, ECMWF, December 1982, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! I. Kirchner, MPI, November 1998, modify nudging
  ! I. Kirchner, MPI, February 1999, nmi
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! 
  ! for more details see file AUTHORS
  ! 

  USE mo_sst,           ONLY: readsst, readice
  USE mo_start_dataset, ONLY: initsd, lres
  USE mo_control,       ONLY: lnmi, lalai, lavgrat, lmidatm, lnudge, ltdiag, &
                              lamip2, lcouple, ldsst, lstratiform
  USE mo_midatm,        ONLY: readozone
  USE mo_machine,       ONLY: machine_setup
  USE mo_nudging,       ONLY: NudgingInit, NDG_INI_MEM
  USE mo_diag_tendency, ONLY: DIAG_Init, IDIAG_INIT
  USE mo_nmi,           ONLY: NMI_Init
  USE mo_legendre,      ONLY: inileg
  USE mo_alai,          ONLY: readalai
  USE mo_avgrat,        ONLY: readavgrat
  USE mo_tracer,        ONLY: reademi
  USE mo_couple,        ONLY: oce_init
  USE mo_dsst,          ONLY: diurnal_init
  USE mo_stratiform,    ONLY: meso_init

  IMPLICIT NONE

  !  External subroutines 
  EXTERNAL drive, initialize, inipost, init_decomposition, init_memory, &
           labrun, readfld, restart, start


  !  Executable statements 

  !-- 0. Print machine specific values

  CALL machine_setup

  !-- 1. Set i/o units and buffer indices

  CALL initsd

  !-- 2. Initialise modules and parallel decomposition

  CALL initialize

  !-- 3. Initialise memory

  CALL init_memory

  !-- 4. Compute *Legendre polynomials and parameters needed
  !      for the *Legendre transforms.

  CALL inileg

  IF (lres) THEN

     !-- 5. Supervise restart

     CALL restart
  ELSE

     !-- 6. Perform initialisation and set up work files

     CALL start
  END IF

  !-- 6.0 initialise aslm (sum of land points per latitude)

  CALL iniaslm

  ! Setup for NMI

  IF (lnmi) CALL NMI_Init

  !-- 6.1 Initialize postprocessing

  CALL inipost

  ! Read optional sst-file

  CALL readsst

  ! Read optional sea ice for AMIP2 type experiment

  IF (lamip2) CALL readice

  ! Read optional lai-file

  IF (lalai) CALL readalai

  ! Read optional vegetation-file

  IF (lavgrat) CALL readavgrat

  ! Read optional ozone-file

  IF (lmidatm .OR. lamip2) THEN
    CALL readozone
  ENDIF

  ! Read optional fields

  CALL reademi

  CALL readfld

  CALL labrun

  ! Setup nudging fields

  IF (lnudge) CALL NudgingInit(NDG_INI_MEM)

  ! Setup of tendency diagnostics

  IF (ltdiag) CALL DIAG_Init(IDIAG_INIT)

  ! Setup ocean

  IF (lcouple) CALL oce_init
  IF (ldsst) CALL diurnal_init
  IF (lstratiform) CALL meso_init

  !-- 7. Start execution

  CALL drive

  RETURN

CONTAINS

  SUBROUTINE iniaslm

    USE mo_memory_g3b,    ONLY: g3b, memory_info, get_info, get_entry
    USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
    USE mo_decomposition, ONLY: dcg => global_decomposition, &
                                 lc => local_decomposition
    USE mo_transpose,     ONLY: gather_gp
    USE mo_gaussgrid,     ONLY: aslm

    IMPLICIT NONE

    TYPE (memory_info) :: info

    REAL, POINTER :: z3d(:,:,:), ptr3d(:,:,:)

    INTEGER :: jlat, irow


    CALL get_entry (g3b, 'SLM', ptr3d)
    CALL get_info (g3b, 'SLM', info)
    ALLOCATE (z3d(info%gdim_1, info%gdim_2, info%gdim_3))
    CALL gather_gp (z3d, ptr3d, dcg)

    IF (p_pe == p_io) THEN
      DO jlat = 1, lc%nlat
        irow = MIN(2*jlat-1,2*(lc%nlat+1-jlat)) ! global ping pong index
        aslm(irow) = SUM(z3d(1:lc%nlon,jlat,1))
      END DO
    END IF

    DEALLOCATE (z3d)

    CALL p_bcast (aslm, p_io)

  END SUBROUTINE iniaslm

END SUBROUTINE control
