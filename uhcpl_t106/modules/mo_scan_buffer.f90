MODULE mo_scan_buffer

  IMPLICIT NONE

  !                                        ! set in  > used in

  REAL, ALLOCATABLE :: dtm_scb    (:,:,:)  ! ffti    ! dyn
  REAL, ALLOCATABLE :: dtl_scb    (:,:,:)  ! ffti    ! dyn
  REAL, ALLOCATABLE :: dalpsl_scb (:,:)    ! ffti    ! dyn
  REAL, ALLOCATABLE :: dalpsm_scb (:,:)    ! ffti    ! dyn
  REAL, ALLOCATABLE :: dm_scb     (:,:,:)  ! scan1sl ! fftd

  REAL, ALLOCATABLE :: vo_scb     (:,:,:)  ! ffti    ! dyn,scan1sl,statd,tf2n
  REAL, ALLOCATABLE :: d_scb      (:,:,:)  ! ffti
  REAL, ALLOCATABLE :: t_scb      (:,:,:)  ! ffti
  REAL, ALLOCATABLE :: alps_scb   (:,:)    ! ffti
  REAL, ALLOCATABLE :: u_scb      (:,:,:)  ! ffti
  REAL, ALLOCATABLE :: v_scb      (:,:,:)  ! ffti
  REAL, ALLOCATABLE :: vol_scb    (:,:,:)            ! fftd
  REAL, ALLOCATABLE :: vom_scb    (:,:,:)            ! fftd
  REAL, ALLOCATABLE :: rh_scb     (:,:,:)            ! fftd
  REAL, ALLOCATABLE :: qe_scb     (:,:,:)
  REAL, ALLOCATABLE :: xe_scb     (:,:,:)
  REAL, ALLOCATABLE :: xte_scb    (:,:,:,:)
  REAL, ALLOCATABLE :: te_scb     (:,:,:)
  REAL, ALLOCATABLE :: alpse_scb  (:,:)
  REAL, ALLOCATABLE :: u0_scb     (:,:)              ! fftd
  REAL, ALLOCATABLE :: du0_scb    (:,:)              ! fftd
  REAL, ALLOCATABLE :: ul_scb     (:,:)              ! fftd
  REAL, ALLOCATABLE :: alnpr_scb  (:,:,:)
  REAL, ALLOCATABLE :: alpha_scb  (:,:,:)
  REAL, ALLOCATABLE :: vervel_scb (:,:,:)

CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE m_bufscan

    USE mo_tracer,        ONLY: ntrac
    USE mo_decomposition, ONLY: dc=>local_decomposition

    LOGICAL, SAVE :: lnot_used = .TRUE.

    INTEGER :: ngl, nlev, nlon, nlpx


    IF (lnot_used) THEN

       ngl  = dc% nglat
       nlev = dc% nlev
       nlon = dc% nglon
       nlpx = dc% nglpx
       ! zero for test_scan_buffer
       ALLOCATE (dtm_scb    (nlpx,nlev,ngl))       ;dtm_scb    = 0.
       ALLOCATE (dtl_scb    (nlpx,nlev,ngl))       ;dtl_scb    = 0.
       ALLOCATE (dalpsl_scb (nlpx,ngl))            ;dalpsl_scb = 0.
       ALLOCATE (dalpsm_scb (nlpx,ngl))            ;dalpsm_scb = 0.
       ALLOCATE (dm_scb     (nlpx,nlev,ngl))       ;dm_scb     = 0.

       ALLOCATE (vo_scb     (nlpx,nlev,ngl))       ;vo_scb     = 0.
       ALLOCATE (d_scb      (nlpx,nlev,ngl))       ;d_scb      = 0.
       ALLOCATE (t_scb      (nlpx,nlev,ngl))       ;t_scb      = 0.
       ALLOCATE (alps_scb   (nlpx,ngl))            ;alps_scb   = 0.
       ALLOCATE (u_scb      (nlpx,nlev,ngl))       ;u_scb      = 0.
       ALLOCATE (v_scb      (nlpx,nlev,ngl))       ;v_scb      = 0.
       ALLOCATE (vol_scb    (nlpx,nlev,ngl))       ;vol_scb    = 0.
       ALLOCATE (vom_scb    (nlpx,nlev,ngl))       ;vom_scb    = 0.
       ALLOCATE (rh_scb     (nlpx,nlev,ngl))       ;rh_scb     = 0.
       ALLOCATE (qe_scb     (nlpx,nlev,ngl))       ;qe_scb     = 0.
       ALLOCATE (xe_scb     (nlpx,nlev,ngl))       ;xe_scb     = 0.
       ALLOCATE (xte_scb    (nlon,nlev,ntrac,ngl)) ;xte_scb    = 0.
       ALLOCATE (te_scb     (nlpx,nlev,ngl))       ;te_scb     = 0.
       ALLOCATE (alpse_scb  (nlpx,ngl))            ;alpse_scb  = 0.
       ALLOCATE (u0_scb     (nlev,ngl))            ;u0_scb     = 0.
       ALLOCATE (du0_scb    (nlev,ngl))            ;du0_scb    = 0.
       ALLOCATE (ul_scb     (nlev,ngl))            ;ul_scb     = 0.
       ALLOCATE (alnpr_scb  (nlpx,nlev,ngl))       ;alnpr_scb  = 0.
       ALLOCATE (alpha_scb  (nlpx,nlev,ngl))       ;alpha_scb  = 0.
       ALLOCATE (vervel_scb (nlpx,nlev,ngl))       ;vervel_scb = 0.

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE m_bufscan
  !------------------------------------------------------------------------------
  SUBROUTINE m_bufsc1

    USE mo_tracer,        ONLY: ntrac
    USE mo_decomposition, ONLY: dc=>local_decomposition
    USE mo_sc1,           ONLY: alnpr, alpha, alps, alpse, d, du0, qe, rh, t, &
                                te, u, u0, ul, v, vervel, vo, vol, vom, xe, xte

    LOGICAL, SAVE :: lnot_used = .TRUE.

    INTEGER       :: nlev, nlon, nlpx


    IF (lnot_used) THEN

       nlev = dc% nlev
       nlon = dc% nglon
       nlpx = dc% nglpx

       ALLOCATE (vo(nlpx,nlev))
       ALLOCATE (d(nlpx,nlev))	
       ALLOCATE (t(nlpx,nlev))
       ALLOCATE (alps(nlpx))
       ALLOCATE (u(nlpx,nlev))
       ALLOCATE (v(nlpx,nlev))
       ALLOCATE (vol(nlpx,nlev))
       ALLOCATE (vom(nlpx,nlev))
       ALLOCATE (rh(nlpx,nlev))
       ALLOCATE (qe(nlpx,nlev))
       ALLOCATE (xe(nlpx,nlev))
       ALLOCATE (xte(nlon,nlev,ntrac))
       ALLOCATE (te(nlpx,nlev))
       ALLOCATE (alpse(nlpx))
       ALLOCATE (u0(nlev))
       ALLOCATE (du0(nlev))
       ALLOCATE (ul(nlev))
       ALLOCATE (alnpr(nlpx,nlev))
       ALLOCATE (alpha(nlpx,nlev))
       ALLOCATE (vervel(nlpx,nlev))

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE m_bufsc1
  !------------------------------------------------------------------------------
  SUBROUTINE m_buftrow (krow)

    USE mo_tracer, ONLY: ntrac
    USE mo_sc1,    ONLY: alnpr, alpha, alps, alpse, d, du0, qe, rh, t, te, &
                         u, u0, ul, v, vervel, vo, vol, vom, xe, xte

    INTEGER, INTENT(IN) :: krow


    vo(:,:)     = vo_scb(:,:,krow)
    d(:,:)      = d_scb(:,:,krow)
    t(:,:)      = t_scb(:,:,krow)
    alps(:)     = alps_scb(:,krow)
    u(:,:)      = u_scb(:,:,krow)
    v(:,:)      = v_scb(:,:,krow)
    vol(:,:)    = vol_scb(:,:,krow)
    vom(:,:)    = vom_scb(:,:,krow)
    rh(:,:)     = rh_scb(:,:,krow)
    qe(:,:)     = qe_scb(:,:,krow)
    xe(:,:)     = xe_scb(:,:,krow)
    IF (ntrac > 0) THEN
       xte(:,:,:)  = xte_scb(:,:,:,krow)
    ENDIF
    te(:,:)     = te_scb(:,:,krow)
    alpse(:)    = alpse_scb(:,krow)
    u0(:)       = u0_scb(:,krow)
    du0(:)      = du0_scb(:,krow)
    ul(:)       = ul_scb(:,krow)
    alnpr(:,:)  = alnpr_scb(:,:,krow)
    alpha(:,:)  = alpha_scb(:,:,krow)
    vervel(:,:) = vervel_scb(:,:,krow)

  END SUBROUTINE m_buftrow

END MODULE mo_scan_buffer
