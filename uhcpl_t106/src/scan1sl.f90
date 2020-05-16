SUBROUTINE scan1sl(pmap,kdpmpf,gauw,cwava,kdpmph,lam,phi,dphi,sinlam,coslam,   &
                   lbasdy,lbasdz,lbassd,lbasiy,detam,detai,dlam,etamid,etaint, &
                   kftype,ub,vb,fb)

  ! Description:
  !
  ! First three scans over the latitude lines in controlling the computations
  ! in Fourier space and the adiabatic part of grid point
  ! calculations including semi-lagrangian transport
  !
  ! Method:
  !
  ! This subroutine reads symmetric and antisymmetric parts
  ! of *Fourier coefficients produced by the scans controlling the
  ! inverse *Legendre transforms. It then scans over the latitude
  ! lines three times to perform the 2nd part of computations in
  ! *Fourier space and computes adiabatic tendencies in grid-point
  ! space including semi-lagrangian transport.
  !
  ! *scan1sl* is called from *stepon*
  !
  ! Externals.:
  ! *ffti*      inverse *fourier transforms.
  ! *tf1*       time filter (part 1)
  ! *tf2*       time filter (part 2)
  ! *dyn*       adiabatic tendencies (except slt-variables)
  ! *sltb1*     semi lagrangian transport
  ! *outgpi*    output grid point variables
  ! *gpc*       grid point computations.
  ! *statd*     statistics for dynamic
  ! *prestat*   prepare statistics and budgets.
  ! *postatd*   complete statistics for dynamics.
  ! *maxwind*   compute maximum wind
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, February 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! A. Rhodin, MPI, July 1999, parallel version    (Transpositions)
  ! I. Kirchner, MPI, May 2000, revision of tendency diagnostics
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_io_tables,     ONLY: ng3xp
  USE mo_constants,     ONLY: a, g, dayl
  USE mo_control,       ONLY: l4ptime, lptime, ltdiag, lwtime, ngl, nrow, &
                              nstop, numfl1, numfl2, twodt, nlon, nresum, &
                              ntbase, ncbase, dtime, lamip2, lcolumn, lens
  USE mo_grid,          ONLY: j1,pcnst,plat,platd,plon,plond,plato2, &
                              plev, plevp1
  USE mo_parallel_slt,  ONLY: nxpt_a
  USE mo_diag_tendency, ONLY: diag_fftd, diag_spectrans
  USE mo_diagnostics,   ONLY: ldiap
  USE mo_doctor,        ONLY: nout
  USE mo_field,         ONLY: field1, field2
  USE mo_gaussgrid,     ONLY: budw
  USE mo_memory_gl,     ONLY: lammp, phimp, q, sigmp, x, xt
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xm1, xtm1
  USE mo_memory_g1b,    ONLY: alpsf, df, qf, tf, vof, xf, xtf
  USE mo_memory_g2a,    ONLY: um1, vm1
  USE mo_memory_g2b,    ONLY: uf, vf
  USE mo_memory_g3b,    ONLY: aps, copy_g3b_to_g3a
  USE mo_grib,          ONLY: outgpi, outgpli
  USE mo_hyb,           ONLY: apsurf
  USE mo_scan_buffer,   ONLY: m_bufscan, m_bufsc1, m_buftrow,               &
                              vo_scb, d_scb, t_scb, alps_scb, u_scb, v_scb, &
                              qe_scb, xte_scb, xe_scb,                      &
                              rh_scb, dm_scb, vom_scb, vol_scb,             &
                              u0_scb, du0_scb, ul_scb
  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_sst,           ONLY: sst
  USE mo_start_dataset, ONLY: nstart, nstep
  USE mo_stat_global,   ONLY: ldiad
  USE mo_tmp_buffer,    ONLY: dm
  USE mo_tracer,        ONLY: ntrac, prestatr, trastat, nstr, source
  USE mo_legendre,      ONLY: rnmd, legmod
  USE mo_year,          ONLY: cd2dat

  USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
  USE mo_decomposition, ONLY: dc=>local_decomposition, &
                              gl_dc => global_decomposition
  USE mo_transpose,     ONLY: gather_gp
  USE mo_CALL_trans,    ONLY: legendre_to_spectral, legendre_to_fourier, &
                              fourier_to_gridpoint, fourier_to_legendre, &
                              gridpoint_to_fourier,                      &
                              test_scan_buffer, test_memory_f
  USE mo_sc1,           ONLY: rh, ul, vol, vom, u0, du0
  USE mo_physc1,        ONLY: cdisse, czen1
  USE mo_global_op,     ONLY: sum_latit_sl
  USE mo_column,        ONLY: get_col_ffti, get_col_tran, cal_col_expl

  USE mo_test_trans

  IMPLICIT NONE

  !  Local parameters: 
  INTEGER, PARAMETER :: itermn = 1, itermx = 4

  !  Scalar arguments 
  REAL    :: cwava, dlam
  INTEGER :: pmap

  !  Array arguments 
  REAL :: coslam(plon), detai(plevp1), detam(plev), dphi(platd,2), &
          etaint(plevp1), etamid(plev), gauw(plat),                &
          lam(plond), lbasdy(4,2,platd,2), lbasdz(4,2,plev),       &
          lbasiy(4,2,platd,2), lbassd(4,2,plevp1), phi(platd,2),   &
          sinlam(plon),                                            &
          ub(plond,plev,platd,2),                                  &
          vb(plond,plev,platd,2),                                  &
          fb(plond,plev,pcnst,platd,2)

  INTEGER :: kdpmpf(pmap), kdpmph(pmap), kftype(pcnst)

  !  Local scalars: 
  REAL    :: hwps, hwpsm1, pscor, psm1cor, ra, zpref, ztmst, ztodt
  INTEGER :: igprow, ilat, jlat, jglat, iter, j, jcen, m, ih, k, i
  INTEGER :: icday, iy, im, id

  ! Bounds on theis PE
  INTEGER :: lnlat, lnlon, lnlpx
  LOGICAL :: col_1d

  !  Local arrays: 
  REAL :: uvb(plond,plev,2,platd,2), alpha(pcnst), hw1(pcnst),             &
          hw1lat(pcnst,plat), hw1p(pcnst), hw2(pcnst), hw2lat(pcnst,plat), &
          hw3(pcnst), hw3lat(pcnst,plat), pdel(plond,plev,plat),           &
          qfcst(plond,plev,pcnst,platd,2), qxl(plond,plev,pcnst,platd,2),  &
          qxr(plond,plev,pcnst,platd,2), uxl(plond,plev,2,platd,2),        &
          uxr(plond,plev,2,platd,2),                                       &
          zpsm1(nlon), zpsm1z(ngl), zpsz(ngl),                             &
          wfld (plon, plevp1, plat),  &! dyn, sltb1
          ulz  (      plev,   plat),  &! dyn, maxwind
          vmaxz(      plev,   plat)    ! dyn, maxwind

  REAL, POINTER :: zalpsm1(:,:)
  REAL, POINTER :: zaps(:,:)

  ! dlam:    increment in x-direction
  ! sinlam:  sin(lamda)
  ! coslam:  cos(lamda)
  ! ub:      u-component of wind
  ! vb:      v-component of wind
  ! fb:      constituent fields
  ! pdel:    pressure differences
  ! wfld:    vertical wind (eta-dot)
  ! uxl:     x-deriv of u/v on left
  ! uxr:     x-deriv of u/v on right
  ! qxl:     x-deriv of const. on L.
  ! qxr:     x-deriv of const. on R.

#ifdef SLDIAG
  REAL :: hqfcst(plond,plev,pcnst,plat)
#endif

  LOGICAL :: ldo_slt

  !  External subroutines 
  EXTERNAL dyn, ewd, sym2, fftd, ffti, gpc, ltd, maxwind,        &
           postatd, postatp, postrad, posts1, prerad, prestat,   &
           qmassa, qmassd,setzeroi, si2, slt2, sltb1,            &
           sltfill, statd, statpz, sym1, tf1, tf2,               &
           make_ensembles

  !  Intrinsic functions 
  INTRINSIC EXP, MOD, SUM


  !  Executable statements 

  ! bounds on theis PE

  lnlon  = dc%nglon
  lnlpx  = dc%nglpx
  lnlat  = dc%nglat
  col_1d = dc%col_1d

!!!!
  ldo_slt = .NOT. col_1d
  !  ldo_slt = .false.
!!!!

  qxl(:,:,:,:,:) = 0.0
  qxr(:,:,:,:,:) = 0.0

  uxl(:,:,:,:,:) = 0.0

  uxr(:,:,:,:,:) = 0.0

  IF (nstep/=nstart) THEN
    ztmst = twodt
  ELSE
    ztmst = 0.5*twodt
  END IF
  ztodt = ztmst
  ra = 1./a

  !-- 1.1 Set parameters

  IF (nstep<=0) THEN
    iter = itermx
  ELSE
    iter = itermn
  END IF

  !-- 1.2 Preset spectral components

  lvo(:,:,:) = 0.0
  ld(:,:,:)  = 0.0
  ltp(:,:,:) = 0.0
  lu0(:,:)   = 0.0

  CALL prestat

!!!!
  ldiad = .FALSE.  ! switch off statistics
  ldiap = .FALSE.
!!!!

  CALL prerad

  !  Check some orbital parameter

  IF (lamip2 .AND. nstep <= nresum+2) THEN
    icday = ncbase + (ntbase+nstep*dtime)/dayl+1.e-6
    CALL cd2dat(icday,id,im,iy)

    WRITE(nout,*) 'date: ',id,'.',im,'.',iy,'  ncbase=   ',ncbase, &
                  '- century day: ',icday,' nstep= ',nstep
    WRITE(nout,*) 'czen1= ',czen1,'  cdisse= ',cdisse
  ENDIF

  IF (ntrac > 0) THEN
    IF (nstep == nresum) CALL prestatr
  END IF

  ! Prepare scan buffer (see bufsc1)

  IF (nstep == nresum) CALL m_bufscan

  ! Allocate space for variables to be saved for 2nd subscan

  CALL m_bufsc1

  ! calculations in legendre space and fft taken out of main loop:

  !-- 2.3 Second part of computations in *fourier space
  !       (symmetric/antisymmetric recombination,
  !        zonal derivatives)

  !-- 2.3.1 Computation of *Fourier components from their symmetric
  !         and antisymmetric parts and retrievial of u and v.

  CALL test_memory_f    ('before sym2')

  CALL sym2

  CALL test_scan_buffer ('after sym2')

  CALL legendre_to_fourier

  !-- 2.3.2 East-west derivatives

  CALL ewd

  !-- 2.4 Inverse fast *fourier transforms

  CALL ffti

  CALL fourier_to_gridpoint

  CALL get_col_ffti

  CALL test_scan_buffer ('after fourier_to_gridpoint')

  !-- 2. Scan control (part 1)

  !-- 2.5 First part of grid point computations

  !-- 2.5.1 Explicit dynamics for eulerian variables

  ! Compute t-dt values at model start.

  IF (nstep==nstart) THEN
    vom1   (:,:,:) =   vo_scb (:,:,:)
    dm1    (:,:,:) =    d_scb (:,:,:)
    qm1    (:,:,:) =    q     (:,:,:)
    xm1    (:,:,:) =    x     (:,:,:)
    tm1    (:,:,:) =    t_scb (:,:,:)
    alpsm1 (:  ,:) = alps_scb (:  ,:)
    um1    (:,:,:) =    u_scb (:,:,:)
    vm1    (:,:,:) =    v_scb (:,:,:)
    IF (ntrac>0) xtm1(:,:,:,:) = xt(:,:,:,:)
  END IF

  ! Blank tendencies

  qe_scb  (:,:,:)   = 0.0
  xe_scb  (:,:,:)   = 0.0
  xte_scb (:,:,:,:) = 0.0

  ! Eulerian advection, energy conversion term
  ! and computation of vertical winds (wfld)

  CALL dyn (wfld, ulz, vmaxz)

  CALL test_scan_buffer ('after dyn')

  ! Second part of time filter

  CALL tf2

  CALL test_scan_buffer ('after tf2')

  !-- 2.5.2  Zonal mass of dry air

  ! has to be adapted to Set A, Set B system ...

  IF (p_pe == p_io) THEN
    ALLOCATE (zalpsm1(nlon, ngl))
    ALLOCATE (zaps   (nlon, ngl))
  END IF

  CALL gather_gp (zalpsm1, alpsm1, gl_dc)
  CALL gather_gp (zaps,    aps,    gl_dc)

  ! index of zpsm1z, zpsz changed: igprow -> local index north to south:

  IF (p_pe == p_io) THEN
    IF (.NOT.col_1d) THEN
      DO jlat = 1, ngl
        !       jglat  = dc%glat(jlat)                 ! global index north -> south
        jglat  = jlat                        !!! sum calculated on whole domain
        igprow = MIN(2*jglat-1,2*(ngl+1-jglat)) ! global ping pong index
        zpsm1(:) = EXP(zalpsm1(1:nlon,jlat)) 
        zpsm1z(igprow) = budw(igprow)/g*SUM(zpsm1)
        zpsz(igprow) = budw(igprow)/g*SUM(zaps(1:nlon,jlat))
      END DO
    ELSE
      zpsm1z = 0.
      zpsz   = 0.
    END IF
    DEALLOCATE (zalpsm1)
    DEALLOCATE (zaps)
  END IF

  CALL p_bcast (zpsm1z, p_io)
  CALL p_bcast (zpsz, p_io)

  IF (col_1d) THEN
    zpsm1z = 0.
    zpsz   = 0.
  END IF

  CALL test_scan_buffer ('before main loop of first subscan')

  !-- 2.1 Start of main loop of first subscan

  DO jlat  = 1, lnlat                       ! local  index north -> south
    jglat  = dc%glat(jlat)                 ! global index north -> south
    igprow = MIN(2*jglat-1,2*(ngl+1-jglat)) ! global ping pong index

    nrow(1) = igprow ! global ping pong index
    nrow(2) = jlat   ! local  index north -> south

    IF (ldiad) CALL statd

    !-- 2.7 End of main loop for first subscan

  END DO

  CALL test_scan_buffer ('after main loop of first subscan')

  !-- 2.5.3 Fill 3d arrays for semi-lagrangian transport

  IF (ldo_slt) THEN
    CALL sltfill(ub,vb,fb,pdel)
  END IF

  !-- 2.8 Scans for slt computations

  IF (ldo_slt) THEN

    ! Initialize extended arrays:
    ! fill latitude and longitude extensions and compute x-derivatives

    CALL sltini(coslam,sinlam,dlam,ub,vb,fb,uxl,uxr,qxl,qxr)

    ! Initialize moisture mass integrals.

    hw1lat(:,:) = 0.0
    hw2lat(:,:) = 0.0
    hw3lat(:,:) = 0.0

    DO ih = 1,2
      DO j = 1,platd
        DO k = 1,plev
          DO i = 1,plond
            uvb(i,k,1,j,ih) = ub(i,k,j,ih)
            uvb(i,k,2,j,ih) = vb(i,k,j,ih)
          END DO
        END DO
      END DO
    END DO

    DO ilat = 1, plat  ! continuous local latitude index for SLT (S->N)
      jlat = plat+1-ilat  ! latitude index for wfld (wfld structure: N->S)
      IF (ilat <= plato2) THEN
        ih = 1
        jcen = j1 - 1 + ilat  ! latitude index for S
      ELSE
        ih = 2
        jcen = j1 - 1 + ilat -plato2  ! latitude index for N
      END IF

      ! Calculate zonal mass of constituents before advected by slt.

      CALL qmassa(cwava,gauw(ilat),fb(1,1,1,jcen,ih),pdel(1,1,ilat), &
           hw1lat(:,jlat),jlat)

      ! Call slt interface routine.

#ifdef SLDIAG
      CALL sltb1(pmap,ih,jcen,ilat,ztodt,ra,iter,                    &
                 uvb(1,1,1,1,ih),                                    &
                 uxl(1,1,1,1,ih),uxr(1,1,1,1,ih),wfld(1,1,jlat),     &
                 fb(1,1,1,1,ih),qxl(1,1,1,1,ih),qxr(1,1,1,1,ih),     &
                 lam,phi(1,ih),dphi(1,ih),etamid,etaint,detam,detai, &
                 lbasdy(1,1,1,ih),lbasdz,lbassd,lbasiy(1,1,1,ih),    &
                 kdpmpf,kdpmph,                                      &
                 lammp(:,:,ilat),phimp(:,:,ilat),sigmp(:,:,ilat),    &
                 qfcst(1,1,1,jcen,ih),nxpt_a(1,jcen,ih),hqfcst(1,1,1,ilat))
#else
      CALL sltb1(pmap,ih,jcen,ilat,ztodt,ra,iter,                    &
                 uvb(1,1,1,1,ih),                                    &
                 uxl(1,1,1,1,ih),uxr(1,1,1,1,ih),wfld(1,1,jlat),     &
                 fb(1,1,1,1,ih),qxl(1,1,1,1,ih),qxr(1,1,1,1,ih),     &
                 lam,phi(1,ih),dphi(1,ih),etamid,etaint,detam,detai, &
                 lbasdy(1,1,1,ih),lbasdz,lbassd,lbasiy(1,1,1,ih),    &
                 kdpmpf,kdpmph,                                      &
                 lammp(:,:,ilat),phimp(:,:,ilat),sigmp(:,:,ilat),    &
                 qfcst(1,1,1,jcen,ih),nxpt_a(1,jcen,ih))
#endif

      ! Calculate zonal mass of constituents after advected by slt.

      CALL qmassa(cwava,gauw(ilat),qfcst(1,1,1,jcen,ih),pdel(1,1,ilat), &
                  hw2lat(:,jlat),jlat)

      ! Compute contribution to global integral hw3.

      CALL qmassd(cwava,gauw(ilat),fb(1,1,1,jcen,ih),qfcst(1,1,1,jcen,ih), &
                  pdel(1,1,ilat),etamid,kftype,hw3lat(:,jlat),jlat)

    END DO

    ! Global mass of dry air

    hwpsm1 = SUM(zpsm1z(1:ngl))
    hwps   = SUM(zpsz(1:ngl))

    hw1(:) = sum_latit_sl ( hw1lat(:,:))
    hw2(:) = sum_latit_sl ( hw2lat(:,:))
    hw3(:) = sum_latit_sl ( hw3lat(:,:))

    DO m = 1, pcnst
      IF (hw3(m)>0.) THEN
        alpha(m) = (hw1(m)-hw2(m))/hw3(m)
      ELSE
        alpha(m) = 0.
      END IF
      hw1p(m) = hw1(m)/hwps
    END DO

    ! Coefficients for air mass correction

    hwpsm1  = hwpsm1*g
    hwps    = hwps*g
    IF (lamip2) THEN
      zpref = apsurf
    ELSE
      zpref = apsurf + 286.
    END IF
    psm1cor = zpref/hwpsm1
    pscor   = zpref/hwps

    IF (ldiad) THEN
      WRITE (nout,'(/,a)') ' Check of air mass correction:' 
      WRITE (nout,'(a,f25.15,9x,f25.15)') '   hw1     = ', hw1
      WRITE (nout,'(a,f25.15,9x,f25.15)') '   hw2     = ', hw2
      WRITE (nout,'(a,f25.15,9x,f25.15)') '   hw3     = ', hw3
      WRITE (nout,'(a,f25.15,9x,f25.15)') '   alpha   = ', alpha
      WRITE (nout,'(a,f25.15,9x,f25.15)') '   hw1/ps  = ', hw1p
      WRITE (nout,'(a,f25.15,a9,f25.15)') '   psm1    = ', hwpsm1, &
                                          ' ps    = ', hwps              
      WRITE (nout,'(a,f25.15,a9,f25.15)') '   psm1cor = ', psm1cor, &
                                          ' pscor = ', pscor

    END IF

  END IF

  CALL test_scan_buffer ('before start of final subscan')

  !-- 3. Final subscan

  !-- 3.4.1 Mass conservation for semi-lagrangian transport

  CALL test_scan_buffer ('before slt2')

  IF (ldo_slt) THEN

#ifdef SLDIAG
    CALL slt2(fb,qfcst,alpha,psm1cor,pscor,etamid,kftype,hqfcst)
#else
    CALL slt2(fb,qfcst,alpha,psm1cor,pscor,etamid,kftype)
#endif

  END IF

  CALL get_col_tran

  CALL test_scan_buffer ('after slt2')

  !-- 3.2 Start of mainloop for final subscan

  DO jlat  = 1, lnlat                         ! local  index north -> south
    jglat  = dc%glat(jlat)                   ! global index north -> south
    igprow = MIN(2*jglat-1, 2*(ngl+1-jglat)) ! global ping pong index

    nrow(1) = igprow ! global ping pong index
    nrow(2) = jlat   ! local  index north -> south

    !-- 3.2.2 Retrieve latitude slice buffer from first subscan

    CALL m_buftrow(jlat)

    !-- 3.4 2nd part of grid point computations

    !-- 3.4.2 Physics and semi-implicit adjustment in grid space

    CALL tf1

    CALL gpc

    !-- explicit integration in column model

    CALL cal_col_expl (ztodt, jlat)

    !-- 4. Grid point contributions to the semi implicit

    CALL si1

    q (:,:,jlat)   =   qm1(:,:,jlat)
    x (:,:,jlat)   =   xm1(:,:,jlat) 
    xt(:,:,:,jlat) =  xtm1(:,:,:,jlat)

    ! copy row buffer to scan buffer to be used by direct fft

    rh_scb  (:,:,jlat) = rh     (:,:)
    dm_scb  (:,:,jlat) = dm     (:,:)
    vom_scb (:,:,jlat) = vom    (:,:)
    vol_scb (:,:,jlat) = vol    (:,:)

    u0_scb  (:  ,jlat) = u0     (:)
    du0_scb (:  ,jlat) = du0    (:)
    ul_scb  (:  ,jlat) = ul     (:)

    DEALLOCATE (dm)

  END DO

  CALL test_scan_buffer ('after end of final subscan')

  !-- 3.5 Direct fast *fourier transforms

  IF (ltdiag) CALL DIAG_fftd
  CALL gridpoint_to_fourier 

  CALL fftd

  !-- 3.6 Semi implicit adjustment (part 2)

  CALL si2

  ! second part of tendency diagnostics called after SI2

  IF (ltdiag) CALL DIAG_SpecTrans

  CALL fourier_to_legendre

  !-- 3.8 Compute contributions to symmetric and
  !       antisymmetric *fourier components.

  CALL sym1

  !-- 4. Direct *Legendre transforms

  CALL ltd

  !-- Transposition Legendre -> spectral

  CALL legendre_to_spectral  

  !-- 5. Complete the scans

  !-- 5.1 Release space

  rnmd(:) = 0.0 ! Legendre fields

  !-- 5.2 Complete statistics

  CALL maxwind(ulz,vmaxz)
  CALL statpz
  CALL postrad
  IF (ldiad) CALL postatd
  IF (ldiap) CALL postatp
  IF (ntrac>0) THEN
    IF (lwtime .OR. (nstep==nstop)) CALL trastat
  END IF

  ! Postprocessing of grid point fields

  IF (lptime) THEN
    IF (.NOT.lcolumn) CALL outgpi
    IF (.NOT.lcolumn) CALL outgpli
    CALL setzeroi
  END IF

  IF (l4ptime) THEN
    ! LK        CALL OUTGPX
    ! LK        CALL SETZEROX
  END IF

  vom1(:,:,:)   = vof(:,:,:)  
  dm1(:,:,:)    = df(:,:,:)   
  tm1(:,:,:)    = tf(:,:,:)   
  alpsm1(:,:)   = alpsf(:,:)  
  qm1(:,:,:)    = qf(:,:,:)   
  xm1(:,:,:)    = xf(:,:,:)   
  xtm1(:,:,:,:) = xtf(:,:,:,:)

  um1(:,:,:)    = uf(:,:,:)
  vm1(:,:,:)    = vf(:,:,:)

  ! Copy g3b buffer to g3a buffer

  CALL copy_g3b_to_g3a

  !-- 5.3 Release buffers and complete I/O processes

  IF (lwtime .OR. nstep==nstop) THEN
    IF (numfl2 > 0)     DEALLOCATE (field2)
    IF (numfl1 > 0)     DEALLOCATE (field1)
    IF (nstr   > 0)     DEALLOCATE (source)
    IF (ALLOCATED(sst)) DEALLOCATE (sst)
  END IF

  IF (lwtime) THEN
    IF (.NOT.lcolumn) THEN
       CALL posts1
!ik++
       IF (lens) THEN
          ! generates ensembles
          CALL make_ensembles
       END IF
!ik++
    END IF
  END IF

END SUBROUTINE scan1sl
