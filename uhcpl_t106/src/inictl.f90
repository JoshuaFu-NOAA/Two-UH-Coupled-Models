!+ preset constants in mo_control.
!+ $Id: inictl.f90,v 1.49 2000/09/19 07:57:50 m214003 Exp $ 

SUBROUTINE inictl

  ! Description:
  !
  ! Preset constants in mo_control.
  !
  ! Method:
  !
  ! Calculate space for memory manager
  !
  ! *inictl* is called from *initialise*.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, August 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! H.-S. Bauer, MPI, Jul 1998, changed
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! I. Kirchner, MPI, November 1998, modify nudging
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception,     ONLY: finish
  USE mo_parameters,    ONLY: jpnlev, jpg3xf
  USE mo_io_tables,     ONLY: ng4xl, ng4xp, ng3xl, ng3xp
  USE mo_post,          ONLY: l4xaccu, lxaccu, n4xpbits, nxpbits
  USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_io, p_bcast
  USE mo_doctor,        ONLY: nin, nout, nerr
  USE mo_control,       ONLY: dtime, eps, labort, lalai, lamip, lanalysis,    &
                              lavgrat, lcolumn, lcouple, ldebug, lg4x, lhg3x, &
                              lmidatm, lnmi, lnudge, lnwp, lpci, lrepro,      &
                              lsstadj, lsub, ltdiag, lvctch, maxrow, n4ptime, &
                              ncbase, ngl, nlev, nlevp1, nn, nproca, nprocb,  &
                              nptime, nctime, nresum, nrow, nscan, nstop, nsub, &
                              nsubint, ntbase, numfl1, numfl2, nvclev, nwlag, &
                              nwtime, twodt, vct, lamip2,imtt, lens, nens, &
                              lonudg,ldsst,lstratiform  !fu++
  USE mo_start_dataset, ONLY: ldebugio, ldebugmem, lg3force, lg3setzero,    &
                              lres, ly365, ndiahdf, nemi, nf1a, nfl1, nfl2, &
                              ng1a, ngl1a, ngribg, ngribs, ngribx, nhf1,    &
                              nhf3, nhf4, nhg1, nhg2, nhg3, nhg4, nhgl1,    &
                              nigp, nini, nisp, nist, njin, njout, nl1a,    &
                              nstart, nstep, ntimeadj, ntwodt, ndstart, ntstart
  USE mo_nudging,       ONLY: NudgingInit, NDG_INI_TIME
  USE mo_constants,     ONLY: inicon, dayl
  USE mo_machine,       ONLY: prec
  USE mo_io,            ONLY: IO_dt, IO_ng3xp, IO_xl, IO_init_dims
  USE m_alloc_mods,     ONLY: alloc_mods ! module subroutine
  USE mo_column,        ONLY: inicolumn

  IMPLICIT NONE

  !  Local scalars: 
  REAL    :: zdtold
  INTEGER :: istpold, js, jx

  !  Intrinsic functions 
  INTRINSIC NINT, ABS

  INCLUDE 'sdsctl.inc'
  INCLUDE 'runctl.inc'


  !  Executable statements 

!-- 1. Preset constants

  nproca = 1  ! assume one processor run as standard,
  nprocb = 1  ! so no decomposition at all ... 

  nscan  = 1

  maxrow = ngl
  labort = .TRUE.

  nstop  = 10
  dtime  = 0.
  eps    = 0.1

  nwlag  = 13
  numfl1 = 0
  numfl2 = 0

  nrow(1) = 1

  nptime  = 0
  nctime  = 0
  n4ptime = 0
  nwtime  = 0

  nsub = 0
  DO js = 1, 9
    nsubint(js) = 1
  END DO
  lcouple = .FALSE.  !fu++
  lonudg  = .FALSE.  !fu++; Ocean nudging experiments
  ldsst   = .FALSE.  !fu++; SST diurnal-cycle 
  lstratiform  = .FALSE.  !fu++; turn on prototype stratiform forcing 

  lamip   = .FALSE.
  lamip2  = .FALSE.
  lsstadj = .FALSE.
  lsub    = .TRUE.
  ldebug  = .FALSE.
  lrepro  = .TRUE.

  lmidatm = .FALSE.

  lpci    = .FALSE.

  lalai   = .FALSE.
  lavgrat = .FALSE.

  lcolumn = .FALSE.
  lvctch  = .FALSE.

  ! NWP forecast mode
  lnwp = .FALSE.

  ! Analysis mode
  lanalysis = .FALSE.

  ! nudging
  lnudge = .FALSE.

  ! tendency diagnostics and nmi
  lnmi   = .FALSE.
  ltdiag = .FALSE.

  ldebug = .FALSE.

  ng3xp = 0
  DO jx = 1, jpg3xf
    ng3xl(jx)   = 1
    lxaccu(jx)  = .FALSE.
    nxpbits(jx) = 16
  END DO
  lhg3x = .FALSE.

  ng4xp = 0
  lg4x  = .FALSE.
  DO jx = 1, jpg3xf
    ng4xl(jx)    = 1
    l4xaccu(jx)  = .FALSE.
    n4xpbits(jx) = 16
  END DO

  lens  = .false.
  nens  = 20

!-- 2. Read namelist runctl

  IF (p_parallel) THEN
     IF (p_parallel_io) THEN
        READ (nin,runctl)
     ENDIF
     CALL p_bcast (nproca, p_io)
     CALL p_bcast (nprocb, p_io)
     CALL p_bcast (nresum, p_io)
     CALL p_bcast (ncbase, p_io)
     CALL p_bcast (ntbase, p_io)
     CALL p_bcast (ng3xp, p_io)
     CALL p_bcast (ng3xl, p_io)
     CALL p_bcast (nxpbits, p_io)
     CALL p_bcast (ng4xp, p_io)
     CALL p_bcast (ng4xl, p_io)
     CALL p_bcast (n4xpbits, p_io)
     CALL p_bcast (n4ptime, p_io)
     CALL p_bcast (nstop, p_io)
     CALL p_bcast (nwtime, p_io)
     CALL p_bcast (nwlag, p_io)
     CALL p_bcast (nptime, p_io)
     CALL p_bcast (nctime, p_io)
     CALL p_bcast (nsub, p_io)
     CALL p_bcast (nsubint, p_io)
     CALL p_bcast (numfl1, p_io)
     CALL p_bcast (numfl2, p_io)
     CALL p_bcast (dtime, p_io)
     CALL p_bcast (lnwp, p_io)
     CALL p_bcast (lanalysis, p_io)
     CALL p_bcast (lnudge, p_io)
     CALL p_bcast (lmidatm, p_io)
     CALL p_bcast (lxaccu, p_io)
     CALL p_bcast (l4xaccu, p_io)
     CALL p_bcast (lcouple, p_io)   !fu++
     CALL p_bcast (lonudg, p_io)    !fu++
     CALL p_bcast (ldsst, p_io)     !fu++
     CALL p_bcast (lstratiform, p_io)     !fu++
     CALL p_bcast (lamip, p_io)
     CALL p_bcast (lamip2, p_io)
     CALL p_bcast (lsstadj, p_io)
     CALL p_bcast (labort, p_io)
     CALL p_bcast (lsub, p_io)
     CALL p_bcast (lhg3x, p_io)
     CALL p_bcast (ldebug, p_io)
     CALL p_bcast (lrepro, p_io)
     CALL p_bcast (lnmi, p_io)
     CALL p_bcast (ltdiag, p_io)
     CALL p_bcast (lpci, p_io)
     CALL p_bcast (lalai, p_io)
     CALL p_bcast (lavgrat, p_io)
     CALL p_bcast (lcolumn, p_io)
     CALL p_bcast (lens, p_io)
     CALL p_bcast (nens, p_io)
  ELSE
     READ (nin,runctl)
  ENDIF

  CALL IO_init_dims

  IF (lcolumn) CALL inicolumn (lcolumn, lvctch, nlev, nlevp1, nvclev, vct)

  CALL alloc_mods  ! previously in IO_init

  ! If in NWP mode ly365 must be true ...

  IF (lnwp) THEN
     ly365 = .TRUE.
  END IF

!-- 3. Initialise the input/output sub-system


!-- 3.2 

  IF (ABS(dtime) < prec) THEN
    IF (nn <= 21) THEN
      dtime = 2400.
    ELSE IF (nn <= 30) THEN
      dtime = 1800.
    ELSE IF (nn <= 42) THEN
      dtime = 1440.
    ELSE IF (nn <= 63) THEN
      dtime = 900.
    ELSE
      dtime = 720.
    END IF
  END IF


  ! Check range of SDS namelist parameter ntimeadj
  IF (ntimeadj >= NINT(dayl)) THEN
    WRITE (nerr,*) ' initsd : ntimeadj = ', ntimeadj
    CALL finish ('initsd', 'SDS namelist parameter ntimeadj out of range')
  END IF

  ! Check if timestep has been changed

  IF (lres) THEN
    zdtold = IO_dt()
    IF (ABS(zdtold-dtime) > 0.) THEN
      istpold = nstep
      nstep   = istpold*zdtold/dtime + 0.01

      WRITE (nout,*) 'Attention: length of timestep has been changed!'
      WRITE (nout,*) 'Old timestep was ', zdtold, ' seconds'
      WRITE (nout,*) 'New timestep is  ', dtime,  ' seconds'
      WRITE (nout,*) 'For consistency reasons number of previous '
      WRITE (nout,*) 'timesteps was automatically changed from ', istpold, &
                     ' to ', nstep, ' steps.'
    END IF
  END IF

  IF (nlev > jpnlev) THEN
    WRITE (nout,*) ' This version of the model does not support'
    WRITE (nout,*) ' more than ', jpnlev, ' vertical levels'
    CALL finish('inictl','Run terminated.')
  END IF

  ! Use g3x information from history files if requested

  IF (lhg3x .AND. lres) THEN
    WRITE (nout,*) 'G3X-information from history files is used'
    ng3xp = IO_ng3xp()
    
    DO jx = 1, ng3xp
      ng3xl(jx) = IO_xl(jx)
    END DO
  END IF

  IF (ng3xp > 0) THEN
    WRITE (nout,'(i2,a,/,a,/)') ng3xp, &
        ' Extra g3 fields have been created', &
        ' the number of levels in these extra fields are '
    WRITE (nout,'(16i5,/)') ng3xl(1:ng3xp)
  END IF

!-- 4. Modify mo_control and write namelist runctl

!-- 4.1 Modify constants in *mo_control*

  IF (nstop < 0) nstop = NINT(-nstop*dayl/dtime)

  twodt  = 2.*dtime
  ntwodt = twodt*1000. + 0.5

  ! moved from setdyn, needed for nudging setup
  CALL inicon

  ! initialise nudging reset initial date
  IF (lnudge) CALL NudgingInit(NDG_INI_TIME)

  ! History and postpro events

  IF (nptime == 0) THEN
    nptime = 12.*3600./dtime + 0.5
  END IF

  IF (nctime == 0) THEN
    nctime = 12.*3600./dtime + 0.5
  END IF

  IF (ng4xp > 0) lg4x = .TRUE.
  IF (n4ptime == 0) n4ptime = nptime

  nresum = nstep

  nens = max(0,nens)

#ifdef DEBUG
!-- 4.2 Write *runctl.*

  WRITE (nerr,sdsctl)
  WRITE (nerr,runctl)
#endif

  IF(p_parallel) THEN
     IF(lnmi)   CALL finish('inictl','NMI works not in parallel mode')
     IF(ltdiag) CALL finish('inictl','tendency diagnostics works not in parallel mode')
  END IF
END SUBROUTINE inictl
