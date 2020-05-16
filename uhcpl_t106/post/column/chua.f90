program chua
  !-------------
  ! modules used
  !-------------
  use mo_system                  ! getarg
  use mo_1d_io, only: wp,&       ! working precision
                      readfile, &!
                      readfile2
  implicit none
  !----------
  ! variables
  !----------
  character(len=128)    :: ufile, afile        ! file names
  logical               :: exist_u, exist_a    ! input files exist
  logical               :: ifmt, ofmt          ! flags for formatted i/o
  integer               :: ios                 ! iostat return value
  integer               :: k                   ! number of commandline args
  integer               :: index(2), lonlat(2)
  character(len=16)     :: name                ! name of variable read
  real(wp) ,pointer     :: x1(:)               ! fieldt input/output
  real(wp) ,pointer     :: x3(:,:,:)           ! fieldt input/output
  !--------------
  ! get arguments
  !--------------
  k = iargc()
print *,k
  if(k/=1) goto 90
  call getarg (1,ufile)
  afile = trim(ufile)//'.asc'
  !---------------------------
  ! mark arguments not present
  !---------------------------
  if(ufile=='') ufile = '?'
  !-----------------------
  ! do input files exist ?
  !-----------------------
  inquire (file=ufile,exist=exist_u)
  inquire (file=afile,exist=exist_a)
  !----------------
  ! print arguments
  !----------------
  print *
  print *,'convert 1d forcing file:'
  if(exist_u.and..not.exist_a) then
    print *,'  from : ',trim(ufile)
    print *,'  to   : ',trim(afile)
    ifmt = .false.
    ofmt = .true.
  else if (exist_a.and..not.exist_u) then
    print *,'  from : ',trim(afile)
    print *,'  to   : ',trim(ufile)
    ifmt = .true.
    ofmt = .false.
  else
    if (.not.exist_u) print *,'         ',trim(ufile),' DOES NOT EXIST'
    if (.not.exist_a) print *,'         ',trim(afile),' DOES NOT EXIST'
    if (     exist_u) print *,'         ',trim(ufile),' EXISTS'
    if (     exist_a) print *,'         ',trim(afile),' EXISTS'
    goto 90
  endif
  !-----------
  ! open files
  !-----------
  if (     ifmt) open(1,file=afile,status='old',form=  'formatted',iostat=ios)
  if (.not.ifmt) open(1,file=ufile,status='old',form='unformatted',iostat=ios)
  if (ios/=0) goto 91
  if (     ofmt) open(2,file=afile,status='new',form=  'formatted',iostat=ios)
  if (.not.ofmt) open(2,file=ufile,status='new',form='unformatted',iostat=ios)
  if (ios/=0) goto 91
  !----------------------------------------
  ! loop: read column - interpolate - write
  !----------------------------------------
  nullify (x1)
  nullify (x3)
  index  = -1
  lonlat = (/0,0/)
  do
    !------------
    ! read column
    !------------
!   call readfile (x1, name, 1, index, lonlat, ifmt=ifmt, ounit=2 ,ofmt=ofmt)
!   if(.not.associated(x1)) goto 92
    call readfile2 (x3, name, 1, ifmt=ifmt, ounit=2 ,ofmt=ofmt)
    if(.not.associated(x3)) goto 92
  end do
  close (1)
  close (2)
  stop
  !-----------------
  ! error conditions
  !-----------------
90 continue
  print *
  !-----------------------------
  ! print error message and exit
  !-----------------------------
  print *
  print *,'USAGE: chua unformatted_file_name'
  stop
91 continue
  print *
  print *,'ERROR'
  stop
92 continue
  print *
  print *,'EOF'
  stop
end program chua
