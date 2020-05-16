PROGRAM createtab

  ! Input: ECHAM data in service format
  ! Output: ASCII table for some codes

  IMPLICIT NONE

  LOGICAL, PARAMETER :: lprint = .FALSE.
  INTEGER, PARAMETER :: maxlev=99
  TYPE srv_header
    INTEGER :: code,level,date,time,nlon,nlat,disp1,disp2
  endtype srv_header

  TYPE dat_header
    INTEGER :: code,levels
    REAL :: xmin(0:maxlev),xmax(0:maxlev),xmean(0:maxlev)
  endtype dat_header

  INTEGER,PARAMETER :: nconst=2
  INTEGER,PARAMETER :: nvari=1
  INTEGER,PARAMETER :: nclat=12
  LOGICAL :: ltime,lnewstep, lsamecole
  TYPE (srv_header)   :: head
  TYPE (dat_header)   :: dcinfo(nconst)
  TYPE (dat_header)   :: dvinfo(nvari)
  CHARACTER(128)  :: ifile,ofile,dummy
  INTEGER :: i,j,ic,nlon,nlat,nscan,nstep,date,time,code,level, code0,level0
  INTEGER :: const(nconst)=(/129, 172/)
  INTEGER :: vari(nvari)=(/130/)
  !  INTEGER :: vari(nvari)=(/130, 142, 143, 156/)
  REAL, ALLOCATABLE   :: field(:,:)
  REAL :: xmin,xmax,xmean
  REAL, ALLOCATABLE   :: glat(:)
  INTEGER :: nclati(nclat)=(/-80,-60,-45,-30,-20,-10,10,20,30,45,60,80/)
  INTEGER :: nvlat(nclat)
  REAL :: zvals(nclat)


  CALL GETFILENAMES

  OPEN (unit=20,file=ifile,form='unformatted',status='old')
  OPEN (unit=30,file=ofile,form='formatted',status='unknown')

  ! nscan=1 : constant codes
  ! nscan=2 : variable codes
  ! nscan=3 : zonal average for -80 -70 ... 70 80
  nscan=1
  DO
    READ(20,END=777) head
    nlon=head%nlon
    nlat=head%nlat
    date=head%date
    time=head%time
    code0=head%code
    level0=head%level
    IF (.NOT. ALLOCATED (field)) ALLOCATE (field(nlon,nlat))
    REWIND(20)
    IF (nscan.EQ.1) THEN
      CALL PUTHEADINFO
      GO TO 777
      WRITE(*,*) 'check constant codes'
      WRITE(30,*) 'constant codes:',nconst
    END IF
    IF (nscan.EQ.2) THEN
      GO TO 777
      WRITE(*,*) 'check variable codes'
      WRITE(30,*) 'variable codes:',nvari
    END IF
    IF (nscan.EQ.3) THEN
      IF (.NOT. ALLOCATED (glat)) ALLOCATE(glat(nlat))
      CALL calclat
      CALL fnlat
      DO i=1,nclat
        IF (lprint) WRITE(*,*) i,nclati(i),nvlat(i)
      ENDDO
      WRITE(*,*) 'make zonal average'
      WRITE(30,'(A)') '#3'
      WRITE(30,'(A)') '# zonal average: '
      WRITE(30,'(A)') '# code level -80    -60    -45    -30    -20    -10     10     20     30     45     60     80'
    END IF

    ltime=.TRUE.
    lnewstep=.TRUE.
    nstep=0

    DO
      READ(20,END=777) head
      code=head%code
      level=head%level
      IF (nlon.NE.head%nlon.OR.nlat.NE.head%nlat) THEN
        WRITE(*,*) nlon,head%nlon,nlat,head%nlat
        STOP 'error in dimensions'
      ENDIF
      lsamecole = code.EQ.code0.AND.level.EQ.level0
      IF (date.NE.head%date.OR.time.NE.head%time.OR.lsamecole) THEN
        date=head%date
        time=head%time
        nstep=nstep+1
        lnewstep=.TRUE.
        IF (lprint) WRITE(*,*) date,head%date,time,head%time
      ENDIF
      READ(20) field
      xmin=MINVAL(field)
      xmax=MAXVAL(field)
      xmean=SUM(field)/(nlon*nlat)
      IF (nscan.EQ.1) THEN
        DO i=1,nconst
          IF (code.EQ.const(i)) EXIT
        END DO
        IF (i.GT.nconst) CYCLE

        IF (nstep.EQ.1) THEN
          IF (lprint) WRITE(*,*) 'code: ',code, 'level: ',level
          CALL PUTCDATINFO
        ELSE
          IF (xmin.NE.dcinfo(i)%xmin(level) .OR. &
               &          xmax.NE.dcinfo(i)%xmax(level) .OR. &
               &          xmean.NE.dcinfo(i)%xmean(level)) THEN
            WRITE(*,*) 'Warning: ',code,level,xmin,dcinfo(i)%xmin(level), &
                 & xmax,dcinfo(i)%xmax(level),xmean,dcinfo(i)%xmean(level)
          END IF
        ENDIF
      ENDIF
      IF (nscan.EQ.2) THEN
        DO i=1,nvari
          IF (code.EQ.vari(i)) EXIT
        END DO
        IF (i.GT.nvari) CYCLE
        CALL PUTVDATINFO
        IF (lprint) WRITE(*,*) 'code: ',code, 'level: ',level
      ENDIF
      IF (nscan.EQ.3) THEN
        DO i=1,nvari
          IF (code.EQ.vari(i)) EXIT
        END DO
        IF (i.GT.nvari) CYCLE
        CALL PUTZDATINFO
        IF (lprint) WRITE(*,*) 'code: ',code, 'level: ',level
      ENDIF

    ENDDO

777 CONTINUE
    WRITE(*,*) 'nstep: ',nstep
    IF (nscan.EQ.3) EXIT
    nscan=nscan+1
    REWIND(20)
  ENDDO


  STOP

CONTAINS

  SUBROUTINE GETFILENAMES

    WRITE(6,*) 'Enter Inputfilename:'
    READ '(a128)',ifile

    DO ic=128,1,-1
      IF (ifile(ic:ic).NE.' ') EXIT
    END DO

    dummy=ifile(1:ic)
    dummy(ic+1:ic+4)='.tab'
    WRITE(*,*) 'Enter Outputfilename (default: ',dummy(1:ic+4),')'
    READ '(a128)',ofile

    DO i=128,1,-1
      IF (ofile(i:i).NE.' ') EXIT
    END DO

    IF (i.EQ.0) ofile = dummy

  END SUBROUTINE GETFILENAMES

  SUBROUTINE PUTHEADINFO
    WRITE(*,'(A)') '#0'
    WRITE(*,'(A)') '# nlon  nlat head8 filename'
    WRITE(*,'(3I6,1X,A)') nlon,nlat,head%disp2,ifile(1:ic)
    WRITE(30,'(A)') '#0'
    WRITE(30,'(A)') '# nlon  nlat head8 filename'
    WRITE(30,'(3I6,1X,A)') nlon,nlat,head%disp2,ifile(1:ic)
  END SUBROUTINE PUTHEADINFO

  SUBROUTINE PUTCDATINFO
    dcinfo(i)%code = code
    dcinfo(i)%levels = dcinfo(i)%levels + 1
    dcinfo(i)%xmin(level)  = xmin
    dcinfo(i)%xmax(level)  = xmax
    dcinfo(i)%xmean(level) = xmean
    WRITE(30,*) code,level,xmin,xmax,xmean
  END SUBROUTINE PUTCDATINFO

  SUBROUTINE PUTVDATINFO
    IF (nstep.EQ.1) THEN
      dvinfo(i)%code = code
      dvinfo(i)%levels = dcinfo(i)%levels + 1
      dvinfo(i)%xmin(level)  = xmin
      dvinfo(i)%xmax(level)  = xmax
      dvinfo(i)%xmean(level) = xmean
    ENDIF
    IF (lnewstep) THEN
      WRITE(*,*) nstep, date, time
      WRITE(30,*) nstep, date, time
      lnewstep = .FALSE.
    ENDIF
    WRITE(30,*) code,level,xmin,xmax,xmean
  END SUBROUTINE PUTVDATINFO

  SUBROUTINE PUTZDATINFO
    IF (nstep.EQ.1) THEN
      dvinfo(i)%code = code
      dvinfo(i)%levels = dcinfo(i)%levels + 1
      dvinfo(i)%xmin(level)  = xmin
      dvinfo(i)%xmax(level)  = xmax
      dvinfo(i)%xmean(level) = xmean
    ENDIF
    IF (lnewstep) THEN
      WRITE(*,*) nstep, date, time
      WRITE(30,*) nstep, date, time
      lnewstep = .FALSE.
    ENDIF
    DO j = 1, nclat
      zvals(j) = -999.
      DO i = 1, nlat
        IF (i.EQ.nvlat(j)) THEN
          zvals(j) = SUM(field(:,i))/nlon
          EXIT
        END IF
      END DO
      IF (zvals(j).EQ. -999.) THEN
        WRITE(*,*) ' putzdatinfo: error',j
        STOP
      END IF
    END DO
    WRITE(30,'(I5,I5,12F7.2)') code,level,zvals
  END SUBROUTINE PUTZDATINFO

  SUBROUTINE calclat
    DO i = 1, nlat
      glat(i) = 90. - (i-1)*(180./nlat) - 90./nlat
    END DO
  END SUBROUTINE calclat

  SUBROUTINE fnlat
    REAL :: zabs, mabs
    INTEGER :: izabs

    DO j = 1, nclat
      zabs  = 90.
      izabs = 0
      DO i = 1, nlat
        mabs = MIN(zabs,ABS(glat(i)-nclati(j)))
        IF (mabs.LT.zabs)  THEN
          zabs=mabs
          izabs=i
        ENDIF
      END DO
      IF (izabs.GT.0) THEN
        nvlat(j)=izabs
      ELSE
        WRITE(*,*) 'fnlat: error'
        STOP
      ENDIF
    END DO

  END SUBROUTINE fnlat
END PROGRAM createtab
