PROGRAM tabinfo
  ! Input: Outputs from createtab

  IMPLICIT NONE

  INTEGER, PARAMETER :: mstep=999
  INTEGER, PARAMETER :: maxlev=99
  INTEGER :: nlev
  TYPE header
     INTEGER :: nlon,nlat,nexp
     CHARACTER(128) :: exp
  END TYPE header
  TYPE dateheader
     INTEGER :: tstep,date,time
  END TYPE dateheader

  TYPE (header)      :: head
  TYPE (dateheader)  :: dhead(mstep,2)
  CHARACTER(128)     :: table1
  INTEGER :: nlon, nlat, i, is, nres, nstep
  INTEGER :: nunit(2)=(/20,21/)
  INTEGER :: steps(2)=(/0,0/)
  INTEGER :: nlevs(2)=(/maxlev,maxlev/)
  REAL :: zfield(12,maxlev,mstep,2)

  CALL GETFILENAME

  OPEN (unit=nunit(1),file=table1,form='formatted',status='old')
  OPEN (unit=77,file='tinfo',form='formatted')

  CALL RSEC0(nunit(1))
  nlon = head%nlon
  nlat = head%nlat

  CALL RSEC3(1)
  nstep = steps(1)

  nlev  = nlevs(1)

  IF (nlat == 32) THEN
     nres = 21
  ELSE IF (nlat == 48) THEN
     nres = 30
  ELSE IF (nlat == 64) THEN
     nres = 42
  ELSE IF (nlat == 96) THEN
     nres = 63
  ELSE IF (nlat == 160) THEN
     nres = 106
  ELSE
     WRITE(*,*) nlon, nlat, nlev, nstep
     STOP 'unsupported resolution'
  ENDIF

  WRITE(77,'(I4,'':'',I4,'':'',I4)') nres, nlev, nstep

CONTAINS

  SUBROUTINE GETFILENAME
    WRITE(6,*) 'Enter table:'
    READ '(a128)',table1
  END SUBROUTINE GETFILENAME

  SUBROUTINE RSEC0(IUNIT)
    INTEGER :: iunit
    CHARACTER(256)  :: dummy

    DO
       READ(iunit,'(A)',END=777) dummy
       IF (dummy(1:2).EQ.'#0') THEN
          DO
             READ(iunit,'(A)') dummy        
             IF (dummy(1:1).NE.'#') THEN
                READ(dummy,*) head%nlon,head%nlat,head%nexp,head%exp
                RETURN
             ENDIF
          ENDDO
       ENDIF
    ENDDO

777 CONTINUE
    WRITE(*,*) 'RSEC0 not found:',IUNIT
  END SUBROUTINE RSEC0

  SUBROUTINE RSEC3(IU)
    INTEGER :: iunit, iu, code, level, klev
    CHARACTER(256)  :: dummy
    LOGICAL :: lread

    iunit=nunit(iu)
    is = 1
    klev = maxlev
    lread = .TRUE.
    DO
       READ(iunit,'(A)',END=777) dummy
       IF (dummy(1:2).EQ.'#3') THEN
          DO
             IF (lread) READ(iunit,'(A)',END=888) dummy
             lread = .TRUE.      
             IF (dummy(1:1).NE.'#') THEN
                READ(dummy,*) dhead(is,iu)%tstep,dhead(is,iu)%date,dhead(is,iu)%time

                DO i=1,klev
                   READ(iunit,'(A)',END=888) dummy
                   IF (is .EQ. 1) THEN     
                      READ(dummy,*) code, level
                      IF (code.NE.130.OR.level.NE.i) THEN
                         nlevs(iu) = i-1
                         klev  = nlevs(iu)
                         lread = .FALSE.      
                         EXIT
                      ELSE
                         nlevs(iu) = i
                         steps(iu) = 1
                      ENDIF
                   ENDIF
                   READ(dummy,*) code, level,zfield(1:12,i,is,iu)
                   IF (code.NE.130.OR.level.NE.i) THEN
                      WRITE(*,*) 'RSEC3: wrong format of table file'
                      STOP
                   ENDIF
                END DO
                steps(iu) = is
                is = is + 1
             ENDIF
          ENDDO
       ENDIF
    ENDDO

777 CONTINUE
    WRITE(*,*) 'RSEC3 not found:',IUNIT
888 CONTINUE
  END SUBROUTINE RSEC3

END PROGRAM tabinfo
