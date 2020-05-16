MODULE mo_call_trans
  !
  !+ $Id: mo_call_trans.f90,v 1.21 1999/11/10 14:18:33 m214030 Exp $
  !
  ! This module holds the routines to invoke the transpositions
  ! with the fields to be transformed as actual parameters
  !
  ! Aditionally the following tests are made
  !
  !  In debug mode (PE 0 handles the whole domain, the other PE's
  !  handle the decomposed domain) the fields on PE0 are compared to those
  !  gathered from the other PE's .
  !
  ! Authors:
  !
  ! A. Rhodin, MPI, August 1999, original source
  !
  USE mo_decomposition, ONLY: dc  => global_decomposition, &
                              dcl =>  local_decomposition, &
                              debug_parallel, any_col_1d
  USE mo_mpi,           ONLY: p_io
  USE mo_doctor,        ONLY: nerr
  USE mo_test_trans,    ONLY: test_spectral, test_legendre, test_gridpoint, &
                              test_symasym,  test_zonmean,  test_row

  USE mo_scan_buffer

  IMPLICIT NONE

  PRIVATE
  !
  ! inverse transpositions
  !
  PUBLIC :: spectral_to_legendre
  PUBLIC :: legendre_to_fourier
  PUBLIC :: fourier_to_gridpoint
  !
  ! direct transpositions
  !
  PUBLIC :: gridpoint_to_fourier
  PUBLIC :: fourier_to_legendre
  PUBLIC :: legendre_to_spectral
  !
  ! test routines
  !
  PUBLIC :: test_memory_f
  PUBLIC :: test_memory_gp
  PUBLIC :: test_scan_buffer
  PUBLIC :: test_row_buffer

CONTAINS
  !==============================================================================
  SUBROUTINE legendre_to_spectral
    USE mo_memory_ls, ONLY: ld, ltp, lvo, lu0 ! Legendre space
    USE mo_memory_sp, ONLY: sd, stp, svo, su0 ! spectral space
    USE mo_transpose, ONLY: tr_ls_sp          ! transposition routine
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (ld, 'ld')
       CALL test_legendre (ltp,'ltp')
       CALL test_legendre (lvo,'lvo')
       CALL test_legendre (lu0,'lu0')
       IF(dcl%pe == p_io) &
            WRITE (nerr,*) 'Test before legendre_to_spectral suceeded.'
    ENDIF
    !
    ! transposition: Legendre -> spectral
    !
    CALL tr_ls_sp (dc, 1, ld, sd, ltp, stp, lvo, svo, lu0, su0)
    !
    ! compare spectral fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_spectral (sd, 'sd')
       CALL test_spectral (stp,'stp')
       CALL test_spectral (svo,'svo')
       CALL test_spectral (su0,'su0')
       IF(dcl%pe == p_io) &
            WRITE (nerr,*) 'Test after legendre_to_spectral suceeded.'
    ENDIF
  END SUBROUTINE legendre_to_spectral
  !------------------------------------------------------------------------------
  SUBROUTINE spectral_to_legendre
    USE mo_memory_ls, ONLY: ld, ltp, lvo, lu0 ! Legendre space
    USE mo_memory_sp, ONLY: sd, stp, svo, su0 ! spectral space
    USE mo_transpose, ONLY: tr_ls_sp          ! transposition routine
    !
    ! compare spectral fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_spectral (sd, 'sd')
       CALL test_spectral (stp,'stp')
       CALL test_spectral (svo,'svo')
       CALL test_spectral (su0,'su0')
       IF(dcl%pe == p_io) &
            WRITE (nerr,*) 'Test before spectral_to_legendre suceeded.'
    ENDIF
    !
    ! transposition: spectral -> Legendre
    !
    CALL tr_ls_sp (dc, -1, ld, sd, ltp, stp, lvo, svo, lu0, su0)
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (ld, 'ld')
       CALL test_legendre (ltp,'ltp')
       CALL test_legendre (lvo,'lvo')
       CALL test_legendre (lu0,'lu0')
       IF(dcl%pe == p_io) &
            WRITE (nerr,*) 'Test after spectral_to_legendre suceeded.'
    ENDIF
  END SUBROUTINE spectral_to_legendre
  !============================================================================
  SUBROUTINE legendre_to_fourier
    USE mo_buffer_fft, ONLY: fftz, fftl,& ! Fourier and Legendre space buffer
         fbm0, lbm0   ! buffer for zonal means (m=0)
    USE mo_transpose,  ONLY: tr_fs_ls     ! transposition routine
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (fftl, 'fftl')
       IF(dcl%pe == p_io) &
            WRITE (nerr,*) 'Test before legendre_to_fourier suceeded.'
    ENDIF
    !
    ! transposition Fourier <- Legendre
    !
    CALL tr_fs_ls (dc, -1, fftz, fftl, fbm0, lbm0)
  END SUBROUTINE legendre_to_fourier
  !----------------------------------------------------------------------------
  SUBROUTINE fourier_to_legendre
    USE mo_buffer_fft, ONLY: fftz, fftl,& ! Fourier and Legendre space buffer
         fbm0, lbm0   ! buffer for zonal means (m=0)
    USE mo_transpose,  ONLY: tr_fs_ls     ! transposition routine
    !
    ! transposition Fourier <- Legendre
    !
    CALL tr_fs_ls (dc, 1, fftz(:,:,:,:6), fftl(:,:,:,:6), &
         fbm0(:,:,1:1),  lbm0(:,:,1:1))
    !
    ! compare legendre fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       CALL test_legendre (fftl(:,:,:,:), 'fftl')
       IF(dcl%pe == p_io) &
            WRITE (nerr,*) 'Test after fourier_to_legendre suceeded.'
    ENDIF
  END SUBROUTINE fourier_to_legendre
  !============================================================================
  SUBROUTINE fourier_to_gridpoint
    USE mo_buffer_fft,  ONLY: fftz, fbm0
    USE mo_scan_buffer, ONLY: d_scb, t_scb, u_scb, v_scb, vo_scb, dtm_scb, &
         dtl_scb, alps_scb, dalpsl_scb, dalpsm_scb,   &
         u0_scb, du0_scb, ul_scb
    USE mo_transpose,   ONLY: tr_gp_fs
    INTEGER :: nlev, nlon
    nlev=dcl% nlev
    nlon=dcl% nlon
    !
    ! transposition grid point <- Fourier
    !
    CALL tr_gp_fs (dc,-1,d_scb,t_scb,u_scb,v_scb,vo_scb,dtm_scb,dtl_scb, &
         sf1=dalpsl_scb, sf2=dalpsm_scb, sf3=alps_scb,   &
         zm1=ul_scb, zm2=u0_scb, zm3=du0_scb,            &
         fs=fftz, fs0=fbm0)
    !
    ! compare gridpoint fields in debug mode
    !
    IF (debug_parallel >= 0) THEN
       IF(any_col_1d) THEN
         IF(dcl%pe == p_io) THEN
           WRITE (nerr,*) 'Test after fourier_to_gridpoint skipped'
           WRITE (nerr,*) 'not tested: d_scb t_scb u_scb v_scb vo_scb dtm_scb'
           WRITE (nerr,*) '            dtl_scb dalpsl_scb dalpsm_scb alps_scb'
         ENDIF
       ELSE
         CALL test_gridpoint (d_scb,     'd_scb')
         CALL test_gridpoint (t_scb,     't_scb')
         CALL test_gridpoint (u_scb,     'u_scb')
         CALL test_gridpoint (v_scb,     'v_scb')
         CALL test_gridpoint (vo_scb,    'vo_scb')
         CALL test_gridpoint (dtm_scb,   'dtm_scb')
         CALL test_gridpoint (dtl_scb,   'dtl_scb')
         CALL test_gridpoint (dalpsl_scb,'dalpsl_scb')
         CALL test_gridpoint (dalpsm_scb,'dalpsm_scb')
         CALL test_gridpoint (alps_scb,  'alps_scb')
         IF(dcl%pe == p_io) &
            WRITE (nerr,*) 'Test after fourier_to_gridpoint suceeded.'
       ENDIF
    ENDIF
  END SUBROUTINE fourier_to_gridpoint
  !------------------------------------------------------------------------------
  SUBROUTINE gridpoint_to_fourier
    USE mo_buffer_fft,  ONLY: fftz, fbm0
    USE mo_scan_buffer, ONLY: rh_scb, dm_scb, vom_scb, vol_scb, &
         u0_scb, du0_scb, ul_scb
    USE mo_memory_g1a,  ONLY: alpsm1, dm1, tm1, vom1
    USE mo_transpose,   ONLY: tr_gp_fs
    INTEGER :: nlev, nlon
    nlev=dcl% nlev
    nlon=dcl% nlon
    !
    ! transposition grid point -> Fourier
    !
    CALL tr_gp_fs (dc, 1,dm1,dm_scb,tm1,rh_scb,vol_scb,vom_scb,vom1,&
         sf3=alpsm1,&
         zm1=ul_scb, zm2=u0_scb, zm3=du0_scb,            &
         fs=fftz, fs0=fbm0)
  END SUBROUTINE gridpoint_to_fourier
  !==============================================================================
  SUBROUTINE test_scan_buffer (text)
    USE mo_decomposition, ONLY: debug_parallel
    !  USE mo_scan_buffer,   ONLY: dtm_scb, dtl_scb, dalpsl_scb, dalpsm_scb,     &
    !                              vo_scb, d_scb, t_scb, alps_scb, u_scb, v_scb, &
    !                              vol_scb, vom_scb, rh_scb, qe_scb, xe_scb,     &
    !                              te_scb, alpse_scb, u0_scb, du0_scb, ul_scb,   &
    !                              alnpr_scb, alpha_scb, vervel_scb
    CHARACTER (len=*) ,INTENT(in) :: text
    IF (debug_parallel>=0) THEN
       CALL test_gridpoint (dtm_scb,   'dtm_scb')  
       CALL test_gridpoint (dtl_scb,   'dtl_scb')
       CALL test_gridpoint (dalpsl_scb,'dalpsl_scb')
       CALL test_gridpoint (dalpsm_scb,'dalpsm_scb')
       CALL test_gridpoint (vo_scb,    'vo_scb')
       CALL test_gridpoint (d_scb,     'd_scb')
       CALL test_gridpoint (t_scb,     't_scb')
       CALL test_gridpoint (alps_scb,  'alps_scb')
       CALL test_gridpoint (u_scb,     'u_scb')
       CALL test_gridpoint (v_scb,     'v_scb')
       CALL test_gridpoint (vol_scb,   'vol_scb')
       CALL test_gridpoint (vom_scb,   'vom_scb')
       CALL test_gridpoint (rh_scb,    'rh_scb')
       CALL test_gridpoint (qe_scb,    'qe_scb')
       CALL test_gridpoint (xe_scb,    'xe_scb')
       CALL test_gridpoint (te_scb,    'te_scb')
       CALL test_gridpoint (alpse_scb, 'alpse_scb')
       CALL test_zonmean   (u0_scb,    'u0_scb', abort=.FALSE.)! not
       CALL test_zonmean   (du0_scb,   'du0_scb',abort=.FALSE.)!(lon,[lev],lat)
       CALL test_zonmean   (ul_scb,    'ul_scb', abort=.FALSE.)! but (lev,lat)
       CALL test_gridpoint (alnpr_scb, 'alnpr_scb')
       CALL test_gridpoint (alpha_scb, 'alpha_scb')
       CALL test_gridpoint (vervel_scb,'vervel_scb')
       IF(dcl%pe == p_io) WRITE (nerr,*) 'Test on scan_buffer suceeded ',text
    ENDIF

    !  REAL, ALLOCATABLE :: xte_scb    (:,:,:,:)

  END SUBROUTINE test_scan_buffer
  !------------------------------------------------------------------------------
  SUBROUTINE test_memory_f (text)
    USE mo_memory_f,      ONLY: f
    USE mo_linked_list,   ONLY: list_element
    USE mo_decomposition, ONLY: debug_parallel
    CHARACTER (len=*) ,INTENT(in) :: text
    TYPE (list_element) ,POINTER :: e
    IF (debug_parallel>=0) THEN
       e => f% first_list_element
       DO
          IF(.NOT.ASSOCIATED(e)) EXIT
          CALL test_symasym (e% field% ptr, e% field% info% name)
          IF(dcl%pe==p_io) &
               WRITE(nerr,*) 'Test on memory_f suceeded: ', e% field% info% name
          e => e% next_list_element
       END DO
       IF(dcl%pe == p_io) WRITE (nerr,*) 'Test on memory_f suceeded ',text
    ENDIF
  END SUBROUTINE test_memory_f
  !------------------------------------------------------------------------------
  SUBROUTINE test_memory_gp (gp, text)
    USE mo_linked_list,   ONLY: list, list_element
    USE mo_decomposition, ONLY: debug_parallel
    TYPE(list)        ,INTENT(in) :: gp
    CHARACTER (len=*) ,INTENT(in) :: text
    TYPE (list_element) ,POINTER :: e
    IF (debug_parallel>=0) THEN
       e => gp% first_list_element
       DO
          IF(.NOT.ASSOCIATED(e))      EXIT
          CALL test_gridpoint (e% field% ptr(:,:,:,:), e% field% info% name)
          IF(dcl%pe==p_io) &
               WRITE(nerr,*) 'Test on memory_gp suceeded: ', e% field% info% name
          e => e% next_list_element
       END DO
       IF(dcl%pe == p_io) WRITE (nerr,*) 'Test on memory_gp suceeded ',text
    ENDIF
  END SUBROUTINE test_memory_gp
  !------------------------------------------------------------------------------
  SUBROUTINE test_row_buffer (j, text)
    USE mo_sc1,        ONLY: vo, d, t, alps, u, v, vol, vom, rh, qe, xe, &
         te, alpse, alnpr, alpha, vervel
    INTEGER          ,INTENT(in) :: j    ! row index
    CHARACTER(len=*) ,INTENT(in) :: text ! comment
    CALL test_row (vo    ,j,'VO '     //text)
    CALL test_row (d     ,j,'D '      //text)
    CALL test_row (t     ,j,'T '      //text)
    CALL test_row (alps  ,j,'ALPS'    //text)
    CALL test_row (u     ,j,'U '      //text)
    CALL test_row (v     ,j,'V '      //text)
    CALL test_row (vol   ,j,'VOL '    //text)
    CALL test_row (vom   ,j,'VOM '    //text)
    CALL test_row (rh    ,j,'RH '     //text)
    CALL test_row (qe    ,j,'QE '     //text)
    CALL test_row (xe    ,j,'XE '     //text)
    CALL test_row (te    ,j,'TE '     //text)
    CALL test_row (alpse ,j,'ALPSE '  //text)
    CALL test_row (alnpr ,j,'ALNPR '  //text)
    CALL test_row (alpha ,j,'ALPHA '  //text)
    CALL test_row (vervel,j,'VERVEL ' //text)
    IF(dcl%pe == p_io) WRITE (nerr,*) 'Test on row_buffer (sc1) suceeded ',text
  END SUBROUTINE test_row_buffer
  !==============================================================================
END MODULE mo_call_trans
