MODULE mo_sst

!+ $Id: mo_sst.f90,v 1.23 1999/11/30 13:52:27 m214003 Exp $

  IMPLICIT NONE

  REAL, ALLOCATABLE :: sst(:,:,:)  ! (nlon,ngl,0:13) in global coordinates
  REAL, ALLOCATABLE :: aice(:,:,:) ! (nlon,ngl,0:13)

CONTAINS

  SUBROUTINE readsst

    ! U. Schlese, DKRZ,  May 1993, original version
    ! U. Schulzweida, MPI, May 1999, netCDF version

    USE mo_start_dataset, ONLY: nist
    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: lamip
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io   
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_filename,      ONLY: nhy
    USE mo_io,            ONLY: sstnc0, sstnc1, sstnc2, io_read, io_open,  &
                                io_close, io_open_unit, io_get_var_double, &
                                io_get_vara_double, io_inq_varid

    REAL, ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL, POINTER :: gl_sst(:,:,:)

    CHARACTER (7) :: fn0, fn1, fn2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), count(3), nvarid
    INTEGER       :: i, nhy0, nhy1, nhy2

    nhy0 = nhy - 1
    nhy1 = nhy
    nhy2 = nhy + 1

    ! Allocate memory for sst per PE

    ALLOCATE (sst(lc%nglon, lc%nglat, 0:13))

    IF (p_pe == p_io) THEN

      IF (nhy < 100) THEN
        WRITE (fn0, '("sst",i2.2)') nhy0
        WRITE (fn1, '("sst",i2.2)') nhy1
        IF(nhy/= 99) THEN
          WRITE (fn2, '("sst",i2.2)') nhy2
        ELSE
          WRITE (fn2, '("sst",i3)') nhy2
        ENDIF
      ELSE IF (nhy< 1000) THEN
        IF (nhy/= 100) THEN
          WRITE (fn0, '("sst",i3)') nhy0
        ELSE
          WRITE (fn0, '("sst",i2.2)') nhy0
        ENDIF
          WRITE (fn1, '("sst",i3)') nhy1
        IF(nhy/= 999) THEN
          WRITE (fn2, '("sst",i3)') nhy2
        ELSE
          WRITE (fn2, '("sst",i4)') nhy2
        ENDIF
      ELSE
        IF(nhy/= 1000) THEN
          WRITE (fn0, '("sst",i4)') nhy0
        ELSE
          WRITE (fn0, '("sst",i3)') nhy0
        ENDIF
          WRITE (fn1, '("sst",i4)') nhy1
          WRITE (fn2, '("sst",i4)') nhy2
      ENDIF

      WRITE (nout, '(/)')

      ! Amip-type:

      IF(lamip) THEN
        WRITE (nout,*)  'This is an AMIP run (lamip = .true.).'
        INQUIRE (file=fn0, exist=lex0)
        INQUIRE (file=fn1, exist=lex1)
        INQUIRE (file=fn2, exist=lex2)
        IF (lex1) THEN
          CALL IO_open (fn1, sstnc1, IO_READ)
          WRITE (nout,*) 'Reading sst from files ',fn0, ', ',fn1,', ',fn2
          IF(lex0) THEN
            CALL IO_open (fn0, sstnc0, IO_READ)
          ELSE
            WRITE (nout,*) 'Could not open file <',fn0,'>'
            CALL finish ('readsst', 'run terminated.')
          ENDIF
          IF(lex2) THEN
            CALL IO_open (fn2, sstnc2, IO_READ)
          ELSE
            WRITE (nout,*) 'Could not open file <',fn2,'>'
            CALL finish ('readsst', 'run terminated.')
          ENDIF
        ELSE
          WRITE (nout,*) 'Could not open file <',fn1,'>'
          CALL finish ('readsst', 'run terminated.')
        ENDIF
      ELSE
        WRITE (nout,*)  'This is no AMIP run (lamip = .false.).'
        INQUIRE (nist, exist=lex)
        IF (lex) THEN
          CALL IO_open_unit (nist, sstnc1, IO_READ)
        ELSE
          WRITE (nout,*) 'Could not open sst file'
          CALL finish ('readsst', 'run terminated.')
        ENDIF
      ENDIF

      ! Read sst-file

      ! Allocate memory for sst global fields

      ALLOCATE (zin(lc%nlon, lc%nlat, 0:13))

      CALL IO_INQ_VARID (sstnc1%nc_file_id, 'sst', nvarid)
      CALL IO_GET_VAR_DOUBLE (sstnc1%nc_file_id, nvarid, zin(:,:,1:12))

      IF(.NOT.lamip) THEN
        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
      ELSE 
        CALL IO_INQ_VARID (sstnc0%nc_file_id, 'sst', nvarid)
        count(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 12 /)
        CALL IO_GET_VARA_DOUBLE (sstnc0%nc_file_id,nvarid,start,count,zin(1,1,0))

        CALL IO_INQ_VARID (sstnc2%nc_file_id, 'sst', nvarid)
        count(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 1 /)
        CALL IO_GET_VARA_DOUBLE (sstnc2%nc_file_id,nvarid,start,count,zin(1,1,13))
      END IF
    END IF

    NULLIFY (gl_sst)
    DO i = 0, 13
      IF (p_pe == p_io) gl_sst => zin(:,:,i:i)
      CALL scatter_gp (gl_sst, sst(:,:,i:i), global_decomposition)
    END DO

    IF (p_pe == p_io) THEN

      DEALLOCATE (zin)

      !    Close file(s)

      CALL IO_close(sstnc1)

      IF(lamip) THEN
        CALL IO_close(sstnc0)
        CALL IO_close(sstnc2)
      ENDIF

    ENDIF

  END SUBROUTINE readsst

  SUBROUTINE readice

    ! U. Schlese, DKRZ,  July 1998, original version (based on opensst)
    ! U. Schulzweida, MPI, Sep 1999, netCDF version

    USE mo_start_dataset, ONLY: nice
    USE mo_doctor,        ONLY: nout
    USE mo_control,       ONLY: lamip
    USE mo_exception,     ONLY: finish
    USE mo_mpi,           ONLY: p_pe, p_io   
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_filename,      ONLY: nhy
    USE mo_io,            ONLY: sicnc0, sicnc1, sicnc2, io_read, io_open,  &
                                io_close, io_open_unit, io_get_var_double, &
                                io_get_vara_double, io_inq_varid

    REAL, ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL, POINTER :: gl_ice(:,:,:)

    CHARACTER (7) :: fn0, fn1, fn2
    LOGICAL       :: lex, lex0, lex1, lex2
    INTEGER       :: start(3), count(3), nvarid
    INTEGER       :: i, nhy0, nhy1, nhy2

    nhy0  = nhy - 1
    nhy1  = nhy
    nhy2  = nhy + 1

    !     Allocate memory for ice per PE

    ALLOCATE (aice(lc%nglon, lc%nglat, 0:13))

    IF (p_pe == p_io) THEN

      IF(nhy < 100) THEN
        WRITE (fn0,'("ice",i2.2)') nhy0
        WRITE (fn1,'("ice",i2.2)') nhy1
        IF(nhy /= 99) THEN
          WRITE (fn2,'("ice",i2.2)') nhy2
        ELSE
          WRITE (fn2,'("ice",i3)') nhy2
        ENDIF
      ELSEIF(nhy < 1000) THEN
        IF(nhy /= 100) THEN
          WRITE (fn0,'("ice",i3)') nhy0
        ELSE
          WRITE (fn0,'("ice",i2.2)') nhy0
        ENDIF
        WRITE (fn1,'("ice",i3)') nhy1
        IF(nhy /= 999) THEN
          WRITE (fn2,'("ice",i3)') nhy2
        ELSE
          WRITE (fn2,'("ice",i4)') nhy2
        ENDIF
      ELSE
        IF(nhy /= 1000) THEN
          WRITE (fn0,'("ice",i4)') nhy0
        ELSE
          WRITE (fn0,'("ice",i3)') nhy0
        ENDIF
        WRITE (fn1,'("ice",i4)') nhy1
        WRITE (fn2,'("ice",i4)') nhy2
      ENDIF

      WRITE (nout, '(/)')

      ! Amip-type:

      IF(lamip) THEN
        ! WRITE (nout,*)  'This is an AMIP run (lamip = .true.).'
        INQUIRE(file=fn0,exist=lex0)
        INQUIRE(file=fn1,exist=lex1)
        INQUIRE(file=fn2,exist=lex2)
        IF (lex1) THEN
          CALL IO_open(fn1, sicnc1, IO_READ)
          WRITE (nout,*) 'Reading ice from files ',fn0, ', ',fn1,', ',fn2
          IF(lex0) THEN
            CALL IO_open(fn0, sicnc0, IO_READ)
          ELSE
            WRITE (nout,*) 'Could not open file <',fn0,'>'
            CALL finish ('readice', 'run terminated.')
          ENDIF
          IF(lex2) THEN
            CALL IO_open(fn2, sicnc2, IO_READ)
          ELSE
            WRITE (nout,*) 'Could not open file <',fn2,'>'
            CALL finish ('readice', 'run terminated.')
          ENDIF
        ELSE
          WRITE (nout,*) 'Could not open file <',fn1,'>'
          CALL finish ('readice', 'run terminated.')
        ENDIF
      ELSE
        ! WRITE (nout,*)  'This is no AMIP run (lamip = .false.).'
        INQUIRE(nice,exist=lex)
        IF (lex) THEN
          CALL IO_open_unit(nice, sicnc1, IO_READ)
        ELSE
          WRITE (nout,*) 'Could not open ice file'
          CALL finish ('readice', 'run terminated.')
        ENDIF
      ENDIF


      ! Read ice-file
 
      ! Allocate memory for sst global fields

      ALLOCATE (zin(lc%nlon, lc%nlat, 0:13))

      CALL IO_INQ_VARID (sicnc1%nc_file_id, 'sic', nvarid)
      CALL IO_GET_VAR_DOUBLE (sicnc1%nc_file_id, nvarid, zin(:,:,1:12))

      IF(.NOT.lamip) THEN
        zin(:,:,0)  = zin(:,:,12)
        zin(:,:,13) = zin(:,:,1)
      ELSE 
        CALL IO_INQ_VARID (sicnc0%nc_file_id, 'sic', nvarid)
        count(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 12 /)
        CALL IO_GET_VARA_DOUBLE (sicnc0%nc_file_id,nvarid,start,count,zin(1,1,0))

        CALL IO_INQ_VARID (sicnc2%nc_file_id, 'sic', nvarid)
        count(:) = (/ lc%nlon, lc%nlat, 1 /)
        start(:) = (/ 1, 1, 1 /)
        CALL IO_GET_VARA_DOUBLE (sicnc2%nc_file_id,nvarid,start,count,zin(1,1,13))
      ENDIF
    END IF

    nullify (gl_ice)
    DO i = 0, 13
      IF (p_pe == p_io) gl_ice => zin(:,:,i:i)
      CALL scatter_gp (gl_ice, aice(:,:,i:i), global_decomposition)
    END DO

    IF (p_pe == p_io) THEN

      DEALLOCATE (zin)

      !    Close file(s)

      CALL IO_close(sicnc1)
      IF(lamip) THEN
        CALL IO_close(sicnc0)
        CALL IO_close(sicnc2)
      ENDIF

    END IF

  END SUBROUTINE readice

  SUBROUTINE clsst

    ! Description:
    !
    ! Interpolates SST in time
    !
    ! Method:
    !
    ! This subroutine calculates the sea-surface-temperatures for
    ! each time step and updates tsm and tsm1m.
    !
    ! *clsst* is called from *gpc*.
    !
    ! Authors: 
    !
    ! U. Schlese, DKRZ, January 1993, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_memory_g3a,    ONLY: auxil1m, auxil2m, slmm, tsm, tsm1m, geospm
    USE mo_control,       ONLY: nrow, lsstadj
    USE mo_start_dataset, ONLY: nstart, nstep
    USE mo_physc2,        ONLY: ctfreez
    USE mo_constants,     ONLY: g
    USE mo_rad_switches,  ONLY: nmonth
    USE mo_decomposition, ONLY: dc=>local_decomposition
    USE mo_timeint,       ONLY: wgt1, wgt2, nmw1, nmw2

    !  Local scalars: 
    REAL    :: zcor, zdt, zgam, zts
    INTEGER :: im, jrow, jn, nglon

    !  Intrinsic functions 
    INTRINSIC MAX


    !  Executable Statements 

    jrow  = nrow(2)
    nglon = dc%nglon 

    ! Update temperatures

    ! Annual cycle

    IF (nmonth == 0) THEN

       ! Interpolation in time

       zdt = 0.01
       DO jn = 1, nglon
          zts = wgt1*sst(jn,jrow,nmw1) + wgt2*sst(jn,jrow,nmw2)
          IF (slmm(jn,jrow) <= 0.5) THEN
             IF (sst(jn,jrow,nmw1) <= ctfreez) THEN
                tsm(jn,jrow)   = MIN(zts,ctfreez-zdt)
                tsm1m(jn,jrow) = tsm(jn,jrow)
             ELSE
                tsm(jn,jrow)   = MAX(zts,ctfreez+zdt)
                tsm1m(jn,jrow) = tsm(jn,jrow)
             END IF
          END IF
       END DO

    ELSE

       ! Perpetual month

       im = nmonth
       DO jn = 1, nglon
          IF (slmm(jn,jrow) <= 0.5) THEN
             tsm(jn,jrow)   = sst(jn,jrow,im)
             tsm1m(jn,jrow) = tsm(jn,jrow)
          END IF
       END DO

    END IF

    ! Adjust sst to sea surface "orography".

    IF (lsstadj) THEN
       ! Lapse rate:
       zgam = 1./100.
       zcor = zgam/g

       DO jn = 1, nglon
          IF (slmm(jn,jrow) < 0.5 .AND. tsm(jn,jrow) > ctfreez) THEN
             tsm(jn,jrow)   = tsm(jn,jrow) - geospm(jn,jrow)*zcor
             tsm1m(jn,jrow) = tsm(jn,jrow)
          END IF
       END DO
    END IF

    ! Initialisation of seaice skin-temperature

    IF (nstep == nstart) THEN

       auxil1m(:,jrow) = tsm(:,jrow)
       auxil2m(:,jrow) = tsm(:,jrow)

    END IF

    RETURN
  END SUBROUTINE clsst

  SUBROUTINE clsst2

    ! Description:
    !
    ! Interpolates SST and sea-ice in time (for AMIP2)
    ! Including adjustment for new calendar
    !
    ! Method:
    !
    ! This subroutine calculates the sea-surface-temperatures and
    ! sea ice cover for each time step 
    ! and updates tsm, tsm1m and seaice, seaicem.
    !
    ! *clsst2* is called from *gpc*.
    !
    ! Authors: 
    !
    ! U. Schlese, DKRZ, July 1998, original source
    ! U. Schulzweida, MPI, October 1999, f90 rewrite
    !

    USE mo_memory_g3a,    ONLY: geospm, tsm, tsm1m, auxil1m, auxil2m, &
                                slmm, seaicem
    USE mo_memory_g3b,    ONLY: seaice
    USE mo_control,       ONLY: nrow, lsstadj
    USE mo_start_dataset, ONLY: nstart, nstep
    USE mo_physc2,        ONLY: ctfreez
    USE mo_constants,     ONLY: g
    USE mo_rad_switches,  ONLY: nmonth
    USE mo_decomposition, ONLY: dc=>local_decomposition
    USE mo_timeint,       ONLY: wgt1, wgt2, nmw1, nmw2

    ! Local scalars: 
    REAL    :: zcor, zdt, zgam, zic, zts
    INTEGER :: im, jrow, jn, nglon

    !  Intrinsic functions
    INTRINSIC MAX


    !  Executable Statements 

    jrow  = nrow(2)
    nglon = dc%nglon

    zdt   = 0.01

    ! Update temperatures and ice

    ! Annual cycle

    IF (nmonth == 0) THEN

       ! Interpolation in time

       DO jn = 1, nglon
          zts = wgt1* sst(jn,jrow,nmw1) + wgt2* sst(jn,jrow,nmw2)
          zic = wgt1*aice(jn,jrow,nmw1) + wgt2*aice(jn,jrow,nmw2)
          IF (slmm(jn,jrow) <= 0.5) THEN
             IF (zic < 50.) THEN
                seaice(jn,jrow)  = 0.
                seaicem(jn,jrow) = 0.
                tsm(jn,jrow)   = MAX(zts,ctfreez+zdt)
                tsm1m(jn,jrow) = tsm(jn,jrow)
             ELSE
                seaice(jn,jrow)  = 1.
                seaicem(jn,jrow) = 1.
                tsm(jn,jrow)   = ctfreez-zdt
                tsm1m(jn,jrow) = tsm(jn,jrow)
             END IF
          END IF
       END DO

    ELSE

       ! Perpetual month

       im = nmonth
       DO jn = 1, nglon
          IF (slmm(jn,jrow) <= 0.5) THEN
             IF (aice(jn,jrow,im) < 50.) THEN
                seaice(jn,jrow)  = 0.
                seaicem(jn,jrow) = 0.
                tsm(jn,jrow)     = MAX(sst(jn,jrow,im),ctfreez+zdt)
                tsm1m(jn,jrow)   = tsm(jn,jrow)
             ELSE
                seaice(jn,jrow)  = 1.
                seaicem(jn,jrow) = 1.
                tsm(jn,jrow)     = ctfreez-zdt
                tsm1m(jn,jrow)   = tsm(jn,jrow)
             ENDIF
          ENDIF
       END DO

    END IF

    ! Adjust SST to sea surface "orography".

    IF (lsstadj) THEN
       ! Lapse rate:
       zgam = 1./100.
       zcor = zgam/g

       DO jn = 1, nglon
          IF (slmm(jn,jrow) < 0.5 .AND. tsm(jn,jrow) > ctfreez) THEN
             tsm(jn,jrow)   = tsm(jn,jrow) - geospm(jn,jrow)*zcor
             tsm1m(jn,jrow) = tsm(jn,jrow)
          END IF
       END DO
    END IF

    ! Initialisation of seaice skin-temperature

    IF (nstep == nstart) THEN

       auxil1m(:,jrow) = tsm(:,jrow)
       auxil2m(:,jrow) = tsm(:,jrow)

    END IF

    RETURN
  END SUBROUTINE clsst2

  SUBROUTINE cplsst

    ! Description:
    !
    ! Interpolates SST in time + coupled SST
    !
    ! Method:
    !
    ! This subroutine calculates the sea-surface-temperatures for
    ! each time step and updates tsm and tsm1m.
    !
    ! *clsst* is called from *gpc*.
    !
    ! Authors: 
    !
    ! U. Schlese, DKRZ, January 1993, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! X. Fu, IPRC, Sep. 2003 
    ! I. Kirchner, IfM, Sep. 2005, ocean coupling

    ! for more details see file AUTHORS
    !

    USE mo_memory_g3a,    ONLY: auxil1m, auxil2m, slmm, tsm, tsm1m, geospm
    USE mo_control,       ONLY: nrow, lsstadj, ldsst
    USE mo_start_dataset, ONLY: nstart, nstep
    USE mo_physc2,        ONLY: ctfreez
    USE mo_constants,     ONLY: g
    USE mo_rad_switches,  ONLY: nmonth
    USE mo_decomposition, ONLY: dc=>local_decomposition
    USE mo_timeint,       ONLY: wgt1, wgt2, nmw1, nmw2

    USE mo_couple,        ONLY: ssto, bzo
    USE mo_dsst,          ONLY: dssto

    !  Local scalars: 
    REAL    :: zcor, zdt, zgam, zts
    INTEGER :: im, jrow, jn, nglon

    !  Intrinsic functions 
    INTRINSIC MAX


    !  Executable Statements 

    jrow  = nrow(2)
    nglon = dc%nglon 

    ! Update temperatures

    ! Annual cycle

    IF (nmonth == 0) THEN

       ! Interpolation in time

       zdt = 0.01
       DO jn = 1, nglon
          zts = wgt1*sst(jn,jrow,nmw1) + wgt2*sst(jn,jrow,nmw2)
          IF (slmm(jn,jrow) <= 0.5) THEN
             IF (sst(jn,jrow,nmw1) <= ctfreez) THEN
                tsm(jn,jrow)   = MIN(zts,ctfreez-zdt)
                tsm1m(jn,jrow) = tsm(jn,jrow)
             ELSE
                tsm(jn,jrow)   = MAX(zts,ctfreez+zdt)
                tsm1m(jn,jrow) = tsm(jn,jrow)
             END IF
          END IF
       END DO

    ELSE

       ! Perpetual month

       im = nmonth
       DO jn = 1, nglon
          IF (slmm(jn,jrow) <= 0.5) THEN
             tsm(jn,jrow)   = sst(jn,jrow,im)
             tsm1m(jn,jrow) = tsm(jn,jrow)
          END IF
       END DO

    END IF

    ! Adjust sst to sea surface "orography".

    IF (lsstadj) THEN
       ! Lapse rate:
       zgam = 1./100.
       zcor = zgam/g

       DO jn = 1, nglon
          IF (slmm(jn,jrow) < 0.5 .AND. tsm(jn,jrow) > ctfreez) THEN
             tsm(jn,jrow)   = tsm(jn,jrow) - geospm(jn,jrow)*zcor
             tsm1m(jn,jrow) = tsm(jn,jrow)
          END IF
       END DO
    END IF

    ! merge with ocean-model SST for bzo=1 fu++
     do jn=1,nglon
     if(bzo(jn,jrow,0).gt.0.8.AND.&
              bzo(jn,jrow,0).lt.1.5) then   !over ocean
     tsm(jn,jrow)=ssto(jn,jrow,0)+273.16

!fu++ turn on diurnal SST 
     IF(ldsst) tsm(jn,jrow)=tsm(jn,jrow)+dssto(jn,jrow,0)

     tsm1m(jn,jrow) = tsm(jn,jrow)
     end if
     end do
!fu++

    ! Initialisation of seaice skin-temperature

    IF (nstep == nstart) THEN

       auxil1m(:,jrow) = tsm(:,jrow)
       auxil2m(:,jrow) = tsm(:,jrow)

    END IF

    RETURN
  END SUBROUTINE cplsst

END MODULE mo_sst
