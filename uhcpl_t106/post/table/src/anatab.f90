PROGRAM anatab
  ! Input: Outputs from createtab

  IMPLICIT NONE

  LOGICAL, PARAMETER :: lprint=.FALSE.
  INTEGER, PARAMETER :: mstep=999
  INTEGER, PARAMETER :: maxlev=99
  INTEGER :: nlev
  INTEGER, PARAMETER :: nclat=12
  INTEGER :: nclati(nclat)=(/-80,-60,-45,-30,-20,-10,10,20,30,45,60,80/)
  TYPE header
     INTEGER :: nlon,nlat,nexp
     CHARACTER(128) :: exp
  END TYPE header
  TYPE dateheader
     INTEGER :: tstep,date,time
  END TYPE dateheader

  TYPE (header)   :: head1, head2, head
  TYPE (dateheader)  :: dhead(mstep,2)
  CHARACTER(128)  :: table1, table2
  INTEGER :: nlon, nlat, i, is, ir, ilat, ilev,j
  INTEGER :: nunit(2)=(/20,21/)
  INTEGER :: steps(2)=(/0,0/)
  INTEGER :: nlevs(2)=(/maxlev,maxlev/)
  REAL :: zfield(12,maxlev,mstep,2), vabs

  CALL GETFILENAMES

  OPEN (unit=nunit(1),file=table1,form='formatted',status='old')
  OPEN (unit=nunit(2),file=table2,form='formatted',status='old')

  ! SEC0
  CALL RSEC0(nunit(1))
  head1=head
  nlon = head1%nlon
  nlat = head1%nlat

  CALL RSEC0(nunit(2))
  head2=head

  IF (nlon.NE.head2%nlon.OR.nlat.NE.head2%nlat) THEN
     WRITE(*,*) nlon,head2%nlon,nlat,head2%nlat
     STOP 'error in dimensions'
  ENDIF

  !  CALL COMPARESEC0

  CALL RSEC3(1)
  WRITE(*,*) 'steps1:',steps(1)
  DO is=1,steps(1)
     DO i=1,nlevs(1)
        IF (lprint) WRITE(*,'(I5,12F7.2)') i,zfield(1:12,i,is,1)
     END DO
  END DO

  CALL RSEC3(2)
  WRITE(*,*) 'steps2:',steps(2)
  DO is=1,steps(2)
     DO i=1,nlevs(2)
        IF (lprint) WRITE(*,'(I5,12F7.2)') i,zfield(1:12,i,is,1)
     END DO
  END DO

  IF (nlevs(1).NE.nlevs(2)) THEN
     WRITE(*,*) nlevs(1), nlevs(2)
     STOP 'error in levels'
  ENDIF

  nlev = nlevs(1)
  WRITE(*,*) 'levels:', nlev

  DO is=1,steps(2)
     DO ir=1,steps(1)
        IF (dhead(is,2)%date.EQ.dhead(ir,1)%date.AND. &
             &  dhead(is,2)%time.EQ.dhead(ir,1)%time) THEN
           IF (is.GT.1) THEN
              IF (dhead(is,2)%date.EQ.dhead(is-1,2)%date.AND. &
                   &  dhead(is,2)%time.EQ.dhead(is-1,2)%time.AND. &
                   &  is.NE.ir) CYCLE
           ENDIF
           WRITE(*,'(A, I4, I10, I4)') 'compare table I: ', &
                        dhead(ir,1)%tstep,dhead(ir,1)%date,dhead(ir,1)%time
           WRITE(*,'(A, I4, I10, I4)') 'with table II  : ', &
                        dhead(is,2)%tstep,dhead(is,2)%date,dhead(is,2)%time
           ilat = 0
           ilev = 0
           vabs = 0.
           DO i=1,nlev
              DO j=1,nclat
                 IF (ABS(zfield(j,i,is,2)-zfield(j,i,ir,1)).GT.vabs) THEN
                    vabs = ABS(zfield(j,i,is,2)-zfield(j,i,ir,1))
                    ilat = j
                    ilev = i
                 END IF
              END DO
           END DO
           IF (vabs.GT.0) THEN
              WRITE(*,'(A,I4,A,I4,A,F7.2)') '  T [K] maximum difference at lat', &
                   & nclati(ilat),' and level ',ilev, ' : ',vabs
           ELSE
              WRITE(*,'(A)') '  T [K] no difference at this timestep'
           ENDIF
           EXIT
        ENDIF
     END DO
  END DO


CONTAINS
  SUBROUTINE GETFILENAMES
    WRITE(6,*) 'Enter table I:'
    READ '(a128)',table1
    WRITE(6,*) 'Enter table II:'
    READ '(a128)',table2
  END SUBROUTINE GETFILENAMES

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

END PROGRAM anatab
