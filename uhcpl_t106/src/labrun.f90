SUBROUTINE labrun

  ! Description:
  !
  ! Label a forecast run.
  !
  ! Method:
  !
  ! Write out details of a forecast run after the set-up is
  ! complete, just before computing the first timestep.
  !
  ! *labrun* has no parameters.
  !
  ! Various items are printed from modules, and the forecast
  ! namelists are written.
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, February 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

#ifdef NAG
  USE f90_unix
#endif

  USE mo_doctor,            ONLY: ylabel1, ylabel2, ylabel3, ylabel4, &
                                  ylabel5, ylabel6, ylabel7, ylabel8, nout
  USE mo_control,           ONLY: ncdata, ncbase, ntdata, ntbase, nlev, ngl, nlon,  &
                                  nstop, twodt, eps, lmidatm, lpci, vct, nvclev, lamip2
  USE mo_truncation,        ONLY: ntrm, ntrn, ntrk
  USE mo_start_dataset,     ONLY: lres
  USE mo_semi_impl,         ONLY: betadt, betazq, apr, tr
  USE mo_param_switches,    ONLY: ndiapfr, lgwdrag, lsurf, lcond, lvdiff, lconv, lice
  USE mo_rad_switches,      ONLY: nradfr, lrad
  USE mo_forecast_switches, ONLY: ndiadfr, ndiavfr
  USE mo_tracer,            ONLY: ntraca, ntrac, nhtrac, lhtrac
  USE mo_year,              ONLY: ic2ymd
  USE mo_filename,          ONLY: yomdn
  USE mo_mpi,               ONLY: p_pe, p_io, p_bcast
  USE mo_exception,         ONLY: finish

  IMPLICIT NONE

#if (! defined NAG) && (! defined CRAY)
  !  Local scalars: 
  INTEGER :: iic

  !  External Functions
  INTEGER, EXTERNAL :: getarg
#endif
#ifdef CRAY
  ! Local scalars:
  INTEGER :: ILEN, ierror

  ! External subroutines:
  EXTERNAL pxfgetarg
#endif


  !  Executable statements

  !-- 0.9  Print name of model

  IF (p_pe == p_io) THEN

#ifdef CRAY
     CALL pxfgetarg(0,yomdn,ILEN,ierror)
#else
#ifndef NAG
     iic = getarg (0, yomdn) 
#else
     CALL getarg (0, yomdn)
#endif
#endif

  END  IF

  CALL p_bcast (yomdn, p_io)

  IF (p_pe == p_io) THEN
     WRITE (nout, '(/,a,a,/)') ' Model: ',TRIM(yomdn) 
  END IF

  !-- 1. Type of run

  IF (p_pe == p_io) THEN

     IF (lres) THEN
        WRITE (nout,'(a)') ' Restarted run.'
        WRITE (nout,'(a)') ' Restart from history files.'
     ELSE
        IF (ncdata==ncbase .AND. ntdata==ntbase) THEN
           WRITE (nout,'(a)') ' Initial run.'
        ELSE
           WRITE (nout,'(/,a,i10,a,i10,a)')          &
                ' Initial data is a forecast from ', &
                ic2ymd(ncbase)*100+ntbase/3600,      &
                ', at forecast time ',               &
                ic2ymd(ncdata)*100+ntdata/3600,      &
                '.'
        END IF
     END IF

     WRITE (nout,'(/,a,/,4(a,i7,/),a,f6.1,/,4(a,i4,/),a,f7.3)') &
          ' General runtime parameter: ',                                        &
          '   number of vertical levels.                          (nlev)    = ', &
          nlev,                                                                  &
          '   number of gaussian latitudes.                       (ngl)     = ', &
          ngl,                                                                   &
          '   max number of points on each latitude line          (nlon)    = ', &
          nlon,                                                                  &
          '   last time step.                                     (nstop)   = ', &
          nstop,                                                                 &
          '   2*dtime.                                            (twodt)   = ', &
          twodt,                                                                 &
          '   frequency of full radiation computations            (nradfr)  = ', &
          nradfr,                                                                &
          '   frequency of dynamical diagnostics                  (ndiadfr) = ', &
          ndiadfr,                                                               &
          '   frequency of vertical dynamical diagnostics         (ndiavfr) = ', &
          ndiavfr,                                                               &
          '   frequency of physics budgets                        (ndiapfr) = ', &
          ndiapfr,                                                               &
          '   time filtering coefficient                          (eps)     = ', &
          eps             

     WRITE (nout,'(2(a,f4.1,/,a,/),2(a,e9.3,/))') &
          '   explicit scheme for d, t, alps (= 0.0)              (betadt)  = ', &
          betadt,                                                                &
          '   semi implicit scheme (= 1.0)                                    ', &
          '   explicit scheme for vo, q (= 0.0)                   (betazq)  = ', &
          betazq,                                                                &
          '   semi implicit scheme (= 1.0)                                    ', &
          '   reference surface pressure for semi-implicit scheme (apr)     = ', &
          apr,                                                                   &
          '   reference temperature for semi-implicit scheme      (tr)      = ', &
          tr

     WRITE (nout,'(a,/,10(a,l2,/))') ' Physical switches: ',                     &
          '   radiation                (lrad)    = ', lrad,                      &
          '   gravity wave drag        (lgwdrag) = ', lgwdrag,                   &
          '   surface exchanges        (lsurf)   = ', lsurf,                     &
          '   middle atmosphere        (lmidatm) = ', lmidatm,                   &
          '   amip2                    (lamip2)  = ', lamip2,                    &
          '   large scale condensation (lcond)   = ', lcond,                     &
          '   large scale cond. (pci)  (lpci)    = ', lpci,                      &
          '   vertical diffusion       (lvdiff)  = ', lvdiff,                    &
          '   convection               (lconv)   = ', lconv,                     &
          '   surface ice              (lice)    = ', lice

     WRITE (nout, '(a)') ' Vertical coordinate table (VCT): '
     WRITE (nout, '(10f7.0)') vct(1:nvclev)
     WRITE (nout, '(10f7.4)') vct(nvclev+1:2*nvclev)
     WRITE (nout, '(a)') ' Max zonal wave number (NTRM): '
     WRITE (nout, '(20i4)') ntrm(:)
     WRITE (nout, '(a)') ' Max meridional wave number for m=0 (NTRN): '
     WRITE (nout, '(20i4)') ntrn(1:nlev)
     WRITE (nout, '(a)') ' Max meridional wave number (NTRK): '  
     WRITE (nout, '(20i4)') ntrk(:)

     WRITE (nout,'(/,8(a,/),a)')    &
          ' -------------------------------------------------------', &
          TRIM(ylabel1), TRIM(ylabel2), TRIM(ylabel3), TRIM(ylabel4), &
          TRIM(ylabel5), TRIM(ylabel6), TRIM(ylabel7), TRIM(ylabel8)


     ! Print copyright for SLT scheme

     WRITE (nout,'(7(a,/))') &
          ' -------------------------------------------------------', &
          ' The semi Lagrangian transport scheme is based on:',       &
          ' NCAR Community Climate Model (CCM2)',                     &
          ' Version 2.1.2 [02/07/94]/, Copyright (C) 1993',           &
          ' University Corporation for Atmospheric Research',         &
          ' All Rights Reserved.',                                    &
          ' -------------------------------------------------------'

     ! Print tracer information

     WRITE (nout, *) &
          'Semi-lagrangian transport of specific humidity, '
     WRITE (nout, *) &
          'cloud water and ', ntrac, ' trace gase(s).'
     WRITE (nout,'(a)') &
          ' -------------------------------------------------------'

     IF (ntraca /= ntrac) THEN
        WRITE (nout,*) ntrac, ' tracers specified.'
     END IF

     IF (lres) THEN
        IF (ntrac>0 .OR. nhtrac>0) THEN
           WRITE (nout, *) ' Rerun: '
           WRITE (nout, *) ' Number of tracers in history files   : ', nhtrac
           WRITE (nout, *) ' Number of tracers in namelist tractl : ', ntrac
           WRITE (nout, *) '  lhtrac = ', lhtrac
           IF (lhtrac) THEN
              ntrac = nhtrac
              WRITE (nout, *) ' Number from history files used.'
              CALL finish ('labrun','lhtrac not implemented')
           ELSE
              WRITE (nout, *) ' Number from namelist tractl used.'
           END IF
        END IF
     END IF

  END IF

END SUBROUTINE labrun
