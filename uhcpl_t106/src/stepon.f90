

SUBROUTINE stepon(pmap,kdpmpf,gauw,cwava,kdpmph,lam,phi,dphi,sinlam,coslam,   &
                  lbasdy,lbasdz,lbassd,lbasiy,detam,detai,dlam,etamid,etaint, &
                  ub,vb,fb,kftype)

  ! Description:
  !
  ! Controls the time step.
  !
  ! Method:
  !
  ! This subroutine controls the structure of the scanning
  ! over the latitude lines and of the computations in spectral
  ! space. It also increments the time step and check for the
  ! completion of the run.
  !
  ! *stepon* is called from *drive*.
  !
  ! Externals:
  ! *scan1sl*   1st scans over gaussian latitudes.
  ! *scan2*     2nd scan over gaussaina latitudes.
  ! *hdiff*     horizontal diffusion.
  ! *scctp*     computations in spectral space for
  !             temperature and surface pressure equations.
  ! *sccd*      computations in spectral space for
  !             divergence equation.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, March 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics, nudging and nmi
  ! I. Kirchner, MPI, January 1999, add nmi
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! for more details see file AUTHORS
  !

  USE mo_memory_sp,     ONLY: stp
  USE mo_doctor,        ONLY: nout, nerr
  USE mo_exception,     ONLY: message, finish
  USE mo_control,       ONLY: dtime, l4ptime, labort, lg4x, lmidatm, lnmi,   &
                              lnudge, lnwp, lptime, lctime, ltdiag, lwtime, lcolumn, &
                              n4ptime, nk, nkp1, nlev, nresum, lamip2,       &
                              nptime, nctime, nstop, nsub, nsubint, nwtime,           &
                              lcouple, ldsst !fu++
  USE mo_hdiff,         ONLY: dampth, difd, dift, difvo
  USE mo_hyb,           ONLY: aktlrd, altrcp, rpr
  USE mo_start_dataset, ONLY: nstart, nstep, ntimeadj
  USE mo_constants,     ONLY: a, dayl
  USE mo_filename,      ONLY: nhm
  USE mo_decomposition, ONLY: dc=>local_decomposition, dcg=>global_decomposition !fu++ 
  USE mo_grid,          ONLY: pcnst,plat,platd,plon,plond,plev,plevp1
  USE mo_parallel_slt,  ONLY: nxpt_a
  USE mo_midatm,        ONLY: uspnge
  USE mo_nudging,       ONLY: Nudging, NudgingReadSST, NudgingOut, NudgingRerun, &
                              NDG_RERUN_WR
  USE mo_grib,          ONLY: outsp, set_output_time, open_grib_file, &
                              close_grib_file

  USE mo_diag_tendency, ONLY: DIAG_Init, DIAG_Write, ldinit,&
                              IDIAG_INI_GBUF, IDIAG_INI_PREV, IDIAG_RER_WR
  USE mo_nmi,           ONLY: NMI_Make, NMI_MAKE_NMI
  USE mo_mpi,           ONLY: p_pe, p_io
  USE mo_timeint,       ONLY: timeintold, timeint
  USE mo_time_control,  ONLY: time_control, TIME_INC_DAYS

  USE mo_couple,        ONLY: oce_run
  USE mo_dsst,          ONLY: sst_diurnal

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: cwava, dlam
  INTEGER :: pmap

  !  Array arguments 
  REAL :: coslam(plon), detai(plevp1), detam(plev), dphi(platd,2), &
          etaint(plevp1), etamid(plev), gauw(plat), lam(plond),    &
          lbasdy(4,2,platd,2), lbasdz(4,2,plev),                   &
          lbasiy(4,2,platd,2), lbassd(4,2,plevp1), phi(platd,2),   &
          sinlam(plon)

  REAL :: ub(plond,plev,platd,2)
  REAL :: vb(plond,plev,platd,2)
  REAL :: fb(plond,plev,pcnst,platd,2)

  INTEGER :: kdpmpf(pmap), kdpmph(pmap), kftype(pcnst)

  !  Local scalars: 
  REAL :: zutime, zstime, zrtime, zwtime
  INTEGER :: id1,id2,iday,im1,iy,iyo,idt, ihm, jlev, jsu, nstdamp !fu++
  LOGICAL :: lprint_tmean

  REAL :: tu0, tu1, ts0, ts1

  !  External functions 
  REAL, EXTERNAL    :: util_walltime
  INTEGER, EXTERNAL :: util_cputime

  !  External subroutines 
  EXTERNAL hdiff, helmo, scan1sl, scan2, sccd, scctp, subjob

  !  Intrinsic functions 
  INTRINSIC INT, MOD

  !  Executable statements

  ! In case of NWP forecast experiments set initial time for GRIB 1 output

  IF (lnwp) THEN
     CALL set_output_time
  END IF

  ! allocate nxpt_a

  ALLOCATE (nxpt_a(plev,platd,2))

  integration_loop: DO

  IF (util_cputime(tu0, ts0) == -1) THEN
     CALL message('stepon','Cannot determine used CPU time')
  END IF
 
!-- 1. 2-scan structure

  ! Set switches for postpro and history events

  IF (nstep > 0) THEN
    IF (ntimeadj >= 0) THEN
      IF (nwtime >= 0) THEN
        lwtime = time_control(nwtime, TIME_INC_DAYS, nstep+1)
      ELSE
        lwtime = time_control(nwtime, TIME_INC_DAYS, nstep+2)
      END IF
    ELSE
      lwtime = time_control(nwtime, TIME_INC_DAYS, nstep+1)
    END IF
  ELSE
    lwtime = .FALSE.
  END IF

  ! Set switch for postpro events

  IF (lcouple) lctime = time_control(nctime, TIME_INC_DAYS, nstep+1)

  lptime = time_control(nptime, TIME_INC_DAYS, nstep+1)
  IF (lg4x)    l4ptime = time_control(n4ptime, TIME_INC_DAYS, nstep+1)

  ! Calculates weighting factores for time interpolation
  IF (lamip2) THEN
    CALL timeint
  ELSE
    CALL timeintold
  END IF

  IF (lnudge) CALL NudgingReadSST

!-- 1.1 Grid point computations and direct *legendre transforms.

  ! First scans over gaussian latitudes

  ! if in NWP mode open new file 

  IF (lnwp .AND. lptime .AND. p_pe == p_io .AND..NOT. lcolumn ) &
    CALL open_grib_file 

  IF (ltdiag) CALL DIAG_Init(IDIAG_INI_GBUF)

  CALL scan1sl(pmap,kdpmpf,gauw,cwava,kdpmph,lam,phi,dphi,sinlam,coslam,   &
               lbasdy,lbasdz,lbassd,lbasiy,detam,detai,dlam,etamid,etaint, &
               kftype,ub,vb,fb)

!-- 1.2 Completion of divergence calculation

  CALL sccd

!-- 1.3 Completion of temperature and surface pressure equations.

  CALL scctp

!-- 1.4 Upper sponge layer 

  IF (lmidatm) CALL uspnge

!-- 1.5 Horizontal diffusion

  IF (lmidatm) THEN
    difvo = a*a/(nk*nkp1*3600.*dampth)
  ELSE
    nstdamp = INT(dayl/dtime)
    IF (nstep<=nstdamp) THEN
      difvo = a*a/(nk*nkp1*3600.*3.)
    ELSE
      difvo = a*a/(nk*nkp1*3600.*dampth)
    END IF
  END IF
  difd = 5.*difvo
  dift = 0.4*difvo

  CALL hdiff

  ! Call nudging after the horizontal diffusion

  IF (lnudge) THEN
     CALL Nudging
  ELSE IF (lnmi) THEN
     CALL NMI_Make(NMI_MAKE_NMI)
  END IF

  ! Postprocessing of spectral data, if necessary close output file

  IF (lptime .AND..NOT. lcolumn) then
    CALL outsp
    IF (lnwp .AND. p_pe == p_io) CALL close_grib_file
  END IF

  IF (ldsst) CALL sst_diurnal
  IF (lcouple) CALL oce_run

  ! Postprocessing of nudging data
  IF (lptime .AND. lnudge) CALL NudgingOut

  IF (ltdiag) THEN
     CALL DIAG_Init(IDIAG_INI_PREV) ! initialize total tendency history
     ldinit = .FALSE.               ! reset initial switch after first time step
     IF (lptime) CALL DIAG_Write    ! perform diagnostics output
  ENDIF

!-- 1.6 Inverse *legendre transforms

  CALL scan2

!-- 2. Continuation of run

!-- 2.1 Recompute matrix *cn* (used to compute divergence
!       in *sccd*) after the first time-step.

  IF (nstep == nstart) THEN
    idt = 2
    CALL helmo(idt)

!-- 2.2 Multiply by 2. arrays *aktlrd* and *altrcp* (used
!       by *conteq*) and *rpr* after the first time step.

    DO jlev = 1, nlev
      aktlrd(jlev) = aktlrd(jlev)*2.
      altrcp(jlev) = altrcp(jlev)*2.
    END DO
    rpr = rpr*2.

  END IF

!-- 2.3   Test for end of run

  IF (nstep == nstop) THEN

    WRITE (nout,*) 'Step ', nstep, ' completed.'
    WRITE (nout,*) 'Experiment finished.'
    IF (util_cputime(zutime, zstime) == -1) THEN
       CALL message('stepon','Cannot determine used CPU time')
    END IF
    zwtime = util_walltime()
    zrtime = (zutime+zstime)/zwtime
    WRITE (nout,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
    WRITE (nout,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
    WRITE (nout,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
    WRITE (nout,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'

    IF (labort) THEN

      WRITE (nout,*) 'Exit to interrupt rerun-chain.'
      CALL finish('stepon','Run terminated.')
    ELSE

      EXIT integration_loop

    END IF

  ELSE IF (lwtime) THEN

    IF (lnudge) CALL NudgingRerun(NDG_RERUN_WR)
    IF (ltdiag) CALL DIAG_Init(IDIAG_RER_WR)  ! store restart data

    WRITE (nout,*) 'Step ', nstep, ' completed.'
    WRITE (nout,*) 'Run finished.'
    IF (util_cputime(zutime, zstime) == -1) THEN
       CALL message('stepon','Cannot determine used CPU time')
    END IF
    zwtime = util_walltime()
    zrtime = (zutime+zstime)/zwtime
    WRITE (nout,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
    WRITE (nout,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
    WRITE (nout,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
    WRITE (nout,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'

    ! Submit jobs

    IF (p_pe == p_io) THEN
      DO jsu = 1, nsub
        IF (nsubint(jsu) < 0) THEN
          ihm = nhm - 1
          IF (ihm < 1) ihm = 12
          IF (MOD(ihm,-nsubint(jsu)) == 0) THEN
            CALL subjob(jsu)
          END IF
        ELSE
          IF (MOD(nhm,nsubint(jsu)) == 0) THEN
            CALL subjob(jsu)
          END IF
        END IF
      END DO
    END IF

    EXIT integration_loop

  END IF

  IF (util_cputime(tu1, ts1) == -1) THEN
    CALL message('stepon','Cannot determine used CPU time')
  END IF

  lprint_tmean = .FALSE.
  IF (dc%nsnm0 > 0) THEN
    IF ( dc%snn0(1) == 0 ) lprint_tmean = .TRUE.
  ENDIF

  IF (nstep <= nresum+NINT(dayl/dtime)) THEN
    IF (lprint_tmean) THEN  
      WRITE (nerr,'(a,i4,a,i7,f10.3," s",f22.13)') &
           ' PE',dc%pe,' stepon: ', nstep, (tu1+ts1)-(tu0+ts0), stp(nlev,1,1)
    ELSE
      WRITE (nerr,'(a,i4,a,i7,f10.3," s")') &
           ' PE',dc%pe,' stepon: ', nstep, (tu1+ts1)-(tu0+ts0)
    ENDIF
  ENDIF 

  ! 2.4   Increment time step

  nstep = nstep + 1

  END DO integration_loop

!ik  close(117)        !fu++

  DEALLOCATE (nxpt_a)

END SUBROUTINE stepon
