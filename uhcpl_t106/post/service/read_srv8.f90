PROGRAM read_srv8

  ! Input:  afterburner data service format 8 byte
  ! Output: contents of the file

  IMPLICIT NONE

  TYPE srv_header
     INTEGER*8 :: code,level,date,time,nlon,nlat,disp1,disp2
  endtype srv_header
  TYPE (srv_header)   :: head

  CHARACTER(128)  :: ifile
  INTEGER :: nlon,nlat,nstep,date,time,code,level
  REAL*8, ALLOCATABLE   :: field(:,:)
  REAL :: xmin,xmax,xmean


  WRITE(6,*) 'Enter Inputfilename:'
  READ '(a128)',ifile

  OPEN (unit=20,file=ifile,form='unformatted',status='old')

  nstep=0

  DO
     READ(20,END=777) head
     IF (nstep .EQ. 0) THEN
        nlon=head%nlon
        nlat=head%nlat
        IF (.NOT. ALLOCATED (field)) ALLOCATE (field(nlon,nlat))
     END IF
     code=head%code
     level=head%level
     date=head%date
     time=head%time
     IF (nlon.NE.head%nlon.OR.nlat.NE.head%nlat) THEN
        WRITE(*,*) nlon,head%nlon,nlat,head%nlat
        STOP 'error in dimensions'
     ENDIF
     READ(20) field
     nstep = nstep + 1
     xmin=MINVAL(field)
     xmax=MAXVAL(field)
     xmean=SUM(field)/(nlon*nlat)
     WRITE (*,'(I5,I4,I7,2I9,2I5,3G12.5)') &
          nstep, code, level, date, time, nlon, nlat, xmin, xmean, xmax
  END DO

777 CONTINUE

  STOP

END PROGRAM read_srv8
