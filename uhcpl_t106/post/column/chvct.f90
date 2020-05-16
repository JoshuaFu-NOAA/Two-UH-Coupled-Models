program chvct
  !-------------
  ! modules used
  !-------------
  use mo_system                ! getarg
  use mo_1d_io, only: wp,&     ! working precision
                      readfile !
  implicit none
  !----------
  ! variables
  !----------
  character(len=128)    :: ifile, ofile, table ! file names
  logical               :: exist_i, exist_t    ! input files exist
  integer               :: ios                 ! iostat return value
  character(len=128)    :: comment             ! vct file comment line
  integer               :: nlevp1              ! number of coefficients
  real(wp) ,allocatable :: ao (:)              ! coefficients A on output
  real(wp) ,allocatable :: bo (:)              ! coefficients B on output
  real(wp) ,allocatable :: ai (:)              ! coefficients A on input
  real(wp) ,allocatable :: bi (:)              ! coefficients B on input
  integer               :: i, k                ! loop index
  character(len=16)     :: name                ! name of variable read
  integer               :: nlevi, nlevo
  integer               :: rank, shap(3), index(2), lonlat(2)
  real(wp) ,pointer     :: xi(:), xo(:)        ! fieldt input/output
  real(wp) ,allocatable :: pi(:), po(:)        ! pressure input/output
  real(wp) ,allocatable :: pmi(:), pmo(:)      ! pressure at model levels
  integer               :: nstep
  real(wp)              :: ps
  !--------------
  ! get arguments
  !--------------
  call getarg (1,ifile)
  call getarg (2,ofile)
  call getarg (3,table)
  !---------------------------
  ! mark arguments not present
  !---------------------------
  if(ifile=='') ifile = '?'
  if(ofile=='') ofile = '?'
  if(table=='') table = '?'
  !-----------------------
  ! do input files exist ?
  !-----------------------
  inquire (file=ifile,exist=exist_i)
  inquire (file=table,exist=exist_t)
  !----------------
  ! print arguments
  !----------------
  print *
  print *,'convert 1d forcing file:'
  print *,'  from : ',trim(ifile)
  print *,'  to   : ',trim(ofile)
  print *,'  table: ',trim(table)
  print *
  !-----------------------------
  ! print error message and exit
  !-----------------------------
  if (.not.exist_i) print *,'         ',trim(ifile),' DOES NOT EXIST'
  if (.not.exist_t) print *,'         ',trim(table),' DOES NOT EXIST'
  if (.not.exist_i .or. .not.exist_t .or. &
      ifile=='?' .or. ofile=='?' .or. table=='?') then
    print *
    print *,'USAGE: chvct ifile ofile table'
    stop
  endif
  !-----------------------------
  ! read hybrid coordinate table
  !-----------------------------
  print *
  open(3,file=table,status='old',iostat=ios); if (ios/=0) goto 93
  read (3,'(a)',iostat=ios) comment;          if (ios/=0) goto 93
  print *,trim(comment)
  read (3,*,iostat=ios) nlevp1;               if (ios/=0) goto 93
  print *,nlevp1,' coefficients'
  if                                         (nlevp1 < 2) goto 93
  allocate (ao (nlevp1),  bo  (nlevp1))
  allocate (po (nlevp1),  pmo (nlevp1-1))
  do i=1,nlevp1
    read (3,*,iostat=ios) k, ao(i), bo(i);    if (ios/=0) goto 93
    write (*,'(i4,2f12.4)') k, ao(i), bo(i)
  end do
  close (3)
  !-----------
  ! open files
  !-----------
  open(1,file=ifile,status='old',form='unformatted',iostat=ios)
                                              if (ios/=0) goto 91
  open(2,file=ofile             ,form='unformatted',iostat=ios)
                                              if (ios/=0) goto 92
  !----------------------------------------
  ! loop: read column - interpolate - write
  !----------------------------------------
  nullify (xi)
  index  = -1
  lonlat = (/0,0/)
  nlevi  = 0
  nlevo  = nlevp1 - 1
  nstep  = 0
  ps     = 100000._wp  
  do
    !------------
    ! read column
    !------------
    call readfile (xi, name, 1, index, lonlat, rank=rank, shap=shap)
    if(.not.associated(xi)) goto 91
    select case (name)
    case ('AK')
      !---------------------
      ! process coefficients
      !---------------------
      nlevi  = size(xi) -1
      allocate (ai (size(xi)), xo (nlevp1))
      ai = xi
      xo = ao
      shap(2) = nlevp1
      print *
      print *, trim(ifile),' : ',size(xi),' coefficients AK read'
    case ('BK')
      nlevi  = size(xi) -1
      allocate (bi (size(xi)), xo (nlevp1))
      allocate (pi (nlevi+1),   pmi(nlevi))
      bi = xi
      xo = bo
      shap(2) = nlevp1
      print *, trim(ifile),' : ',size(xi),' coefficients BK read'
      print *
      print *,'  Old Levels'
      print *,nlevi+1,' coefficients'
      do i=1,nlevi+1
        write (*,'(i4,2f12.4)') i, ai(i), bi(i)
      end do
    print *
    case default 
      if(size(xi)==nlevi .and. nlevi == shap(2)) then
        !------------------------------
        ! interpolate
        !   source : xi (1:nlevi)
        !   dest   : xo (1:nlevo)
        ! coefficients
        !   source : ai, bi (1:nlevi+1)
        !   dest   : ao, bo (1:nlevo+1)
        ! pressure in between model levels
        !   source : pi (1:nlevi+1)
        !   dest   : po (1:nlevo+1)
        ! pressure at model levels
        !   source : pmi (1:nlevi)
        !   dest   : pmo (1:nlevo)
        !------------------------------
!print *, 'interpolate ', name, size(xi), size(xo), ps
        allocate (xo (nlevo))
        shap(2) = nlevo
        pi  = ai + bi * ps
        po  = ao + bo * ps
        pmi = (pi(1:nlevi)+pi(2:nlevi+1))/2._wp
        pmo = (po(1:nlevo)+po(2:nlevo+1))/2._wp

        do k=1,nlevo
          if (pmi(1) > pmo(k)) then
            xo(k)=xi(1)
!i=2
          elseif (pmi(nlevi) < pmo(k)) then
            xo(k)=xi(nlevi)
!i=nlevi+1
          else
            i=1
            do
              if (pmi(i) <= pmo(k) .and. i < nlevi) then
                i=i+1
              else
                exit
              endif
            end do
            if (i > nlevi) then
              print *,'Mistake in interpolation'
              stop
            else
             xo(k) = xi(i)*(pmo(k)-pmi(i-1))/(pmi(i)-pmi(i-1)) &
                   + xi(i-1)*(pmi(i)-pmo(k))/(pmi(i)-pmi(i-1))
            endif
          endif
!print *,xo(k),pmo(k),ao(k), bo(k) * ps !,xi(i),pmi(i)
!print *,xo(k),pmo(k),i-1,xi(i-1),pmi(i-1)
        enddo
!stop

      else
        !----------------
        ! leave unchanged
        !----------------
        allocate (xo (size(xi)))
        xo = xi
!print *, 'copy        ', name, size(xi), size(xo)
        !---------------------------------
        ! special tretment for nstep, alps
        !---------------------------------
        select case (name)
        case ('NSTEP')
          nstep = xi(1)
          ps    = 100000._wp
        case ('ALPS')
          ps    = exp(xi(1))
        end select
      endif      
    end select
    !------
    ! write
    !------
    write (2) name, rank, shap
    write (2) xo
    deallocate (xo)
  end do
  close (1)
  close (2)
  stop
  !-----------------
  ! error conditions
  !-----------------
91 continue
  print *
  print *,'EOF reading ',trim(ifile)
  stop
92 continue
  print *
  print *,'ERROR writing ',trim(ofile)
  stop
93 continue
  print *
  print *,'ERROR reading table ',trim(table)
contains
  subroutine interpolate
  end subroutine interpolate
end program chvct
