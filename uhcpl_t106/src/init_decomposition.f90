SUBROUTINE init_decomposition
!
!+ $Id: init_decomposition.f90,v 1.22 2000/03/24 14:37:45 m214030 Exp $
!

  USE mo_doctor,        ONLY: nout
  USE mo_control,       ONLY: nproca, nprocb, ngl, nlon, nlev, nm, nn, nk, &
                              lcolumn
  USE mo_exception,     ONLY: finish
  USE mo_mpi,           ONLY: p_nprocs, p_pe, p_io, p_parallel
  USE mo_decomposition, ONLY: local_decomposition, &
                              global_decomposition, &
                              print_decomposition, &
                              decompose, debug_seriell
  USE mo_transpose,     ONLY: indx
  USE mo_column,        ONLY: setcolumn, lat_1d, lon_1d !, inicolumn

  IMPLICIT NONE

  INTEGER :: p,i
  
  LOGICAL :: lrot, lfull_m
  INTEGER :: debug

  lrot          = .TRUE.  ! true: no rotation of longitudes
  lfull_m       = .FALSE. ! true: full m-columns per PE
  debug_seriell = .FALSE. ! true: same results as ser.version if nprocb == 1

  ! debug = 0,1 : PE0 takes full domain (always no rotation)
  !         -1  : no special treatment of PE 0
  !          0  : gather from PE 0
  !          1  : gather from PEs > 0

  IF (p_nprocs == nproca*nprocb) THEN
     debug = -1
  ELSE IF ( p_nprocs == nproca*nprocb+1) THEN
     debug = 0
  ELSE
     CALL finish ('init_decomposition', &
&       'Total number of PEs on runtime doesn''t fit nproca*nprocb(+1)')
  END IF

  IF (p_parallel) THEN 
     WRITE (nout,'(/,a,i4,a,i3,a,i3)') &
          ' Total number of PEs: ', p_nprocs, &
          ' set A: ', nproca, ' set B: ', nprocb
  ENDIF

  ALLOCATE (global_decomposition(1:p_nprocs))

  ! read column model namelist 'columnctl'

!  IF (lcolumn) CALL inicolumn (lcolumn)

  ! derive decomposition

  IF (lcolumn) THEN
    CALL decompose (global_decomposition, nproca, nprocb, ngl, nlon, nlev, &
         nm, nn, nk, norot=lrot, debug=debug, lfull_m=lfull_m, &
         lats_1d=lat_1d(1), lons_1d=lon_1d(1))
  ELSE
    CALL decompose (global_decomposition, nproca, nprocb, ngl, nlon, nlev, &
         nm, nn, nk, norot=lrot, debug=debug, lfull_m=lfull_m)
  ENDIF

  ! keep index values, not id values

  DO p = 1, p_nprocs
     global_decomposition(p)%mapmesh(:,:) = &
&          indx(global_decomposition(p)%mapmesh(:,:),global_decomposition)
  END DO

  ! copy global decomposition table entry to local decomposition

  DO i = 1, p_nprocs
     IF (global_decomposition(i)% pe == p_pe) THEN
        local_decomposition = global_decomposition(i)
     END IF
  END DO

  ! decomposition printout

  IF (p_pe == p_io .AND. p_parallel) THEN
     DO i = 1, p_nprocs
        WRITE (nout,'(78("-"))')
        CALL print_decomposition(global_decomposition(i))
     END DO
     WRITE (nout,'(78("-"),/)')
  END IF

  ! initialize column model

  CALL setcolumn

END SUBROUTINE init_decomposition
