module mo_1d_io
  !-------------
  ! modules used
  !-------------
  use mo_system       ! abort
  use ask_module, only: ask
  implicit none
  !----------------
  ! public entities
  !----------------
  private
  public :: wp        ! working precision kind parameter
  public :: sp        ! single precision
  public :: readfile  ! subroutine to read/copy column model files (1 column)
  public :: readfile2 ! subroutine to read/copy whole column model file 
  !----------------
  ! Kind parameters
  !----------------
  integer, parameter  :: dp = selected_real_kind(13) ! Double Precision
  integer, parameter  :: sp = selected_real_kind(6)  ! Single Precision
  integer, parameter  :: wp = dp
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  subroutine readfile (y ,name ,iunit ,index ,lonlat &
                      ,ounit ,ifmt ,ofmt, rank, shap)
  !-------------------------------------------
  ! subroutine to read/copy column model files
  !-------------------------------------------
  real(wp)          ,pointer       :: y (:)     ! field (column) to return
  character(len=16) ,intent(out)   :: name      ! name of field read
  integer           ,intent(in)    :: iunit     ! Fortran unit to read from
  integer           ,intent(inout) :: index(:)  ! index of column file
  integer           ,intent(inout) :: lonlat(2) ! lon,lat index in echam
  integer ,optional ,intent(in)    :: ounit     ! Fortran unit to copy to
  logical ,optional ,intent(in)    :: ifmt      ! true: formatted inputfile
  logical ,optional ,intent(in)    :: ofmt      ! true: formatted outputfile
  integer ,optional ,intent(out)   :: rank      ! rank of array read
  integer ,optional ,intent(out)   :: shap(3)   ! shape of array read
    real(wp) ,allocatable :: z (:,:)
    integer               :: rnk, shp(3), siz, i
    logical               :: lifmt, lofmt
    integer               :: ios
    lifmt = .false. ;if(present (ifmt)) lifmt = ifmt
    lofmt = .false. ;if(present (ofmt)) lofmt = ofmt
    if(present (rank)) rank = 0
    if(present (shap)) shap = 0
    if(associated(y)) deallocate (y)
    if(lifmt) then
      read(iunit,*,iostat=ios) name, rnk, shp
    else
      read(iunit,iostat=ios) name, rnk, shp
    endif
    if(ios/=0) return
    siz = product (shp(2:rnk))
    allocate(y(siz))
    allocate(z(shp(1),siz))
    if(lifmt) then
      read(iunit,*)z
    else
      read(iunit)z
    endif
    if (shp(1)==1) then
      y = z (1,:)
    else
      if (index(iunit) < 1) then
        if(name/='ILON_JLAT') call abort ('readfile: undefined index')
lola:   do
          call ask(lonlat,'enter lon.,lat. indices')
          do i = 1, shp(1)
            write(0,'(2f5.0)') z(i,:)
            if (all(lonlat==z(i,:))) then
              index(iunit)=i
              exit lola
            endif
          end do
        end do lola
      endif
      if (index(iunit) > shp(1)) call abort ('readfile: index too large')
      y = z (index(iunit),:)
    endif
    !-----------------------------------
    ! rewrite file if 'ounit' is present
    !-----------------------------------
    if (present(ounit)) then
      if (lofmt) then
        write (ounit,*)         name, rnk, 1, shp(2:3)
        if (shp(1)==1) then
          write (ounit,'(g22.15)') z (1,:)
        else
          write (ounit,'(g22.15)') z (index(iunit),:)
        endif
      else
        write (ounit) name, rnk, 1, shp(2:3)
        if (shp(1)==1) then
          write (ounit) z (1,:)
        else
          write (ounit) z (index(iunit),:)
        endif
      endif
    endif
    if(present (rank)) rank = rnk
    if(present (shap)) shap = (/1,shp(2),shp(3)/)
  end subroutine readfile
!------------------------------------------------------------------------------
  subroutine readfile2 (y ,name ,iunit ,ounit ,ifmt ,ofmt, rank, shap)
  !-------------------------------------------
  ! subroutine to read/copy column model files
  !-------------------------------------------
  real(wp)          ,pointer       :: y (:,:,:) ! field (column) to return
  character(len=16) ,intent(out)   :: name      ! name of field read
  integer           ,intent(in)    :: iunit     ! Fortran unit to read from
  integer ,optional ,intent(in)    :: ounit     ! Fortran unit to copy to
  logical ,optional ,intent(in)    :: ifmt      ! true: formatted inputfile
  logical ,optional ,intent(in)    :: ofmt      ! true: formatted outputfile
  integer ,optional ,intent(out)   :: rank      ! rank of array read
  integer ,optional ,intent(out)   :: shap(3)   ! shape of array read
    integer               :: rnk, shp(3)
    logical               :: lifmt, lofmt
    integer               :: ios
    lifmt = .false. ;if(present (ifmt)) lifmt = ifmt
    lofmt = .false. ;if(present (ofmt)) lofmt = ofmt
    if(present (rank)) rank = 0
    if(present (shap)) shap = 0
    if(associated(y)) deallocate (y)
    if(lifmt) then
      read(iunit,*,iostat=ios) name, rnk, shp
    else
      read(iunit,iostat=ios) name, rnk, shp
    endif
    if(ios/=0) return

    select case (rnk)
    case (1)
      allocate (y(shp(1),1,1))
    case (2)
      allocate (y(shp(1),shp(2),1))
    case (3)
      allocate (y(shp(1),shp(2),shp(3)))
    end select
    if(lifmt) then
      read(iunit,*) y
    else
      read(iunit) y
    endif
    !-----------------------------------
    ! rewrite file if 'ounit' is present
    !-----------------------------------
    if (present(ounit)) then
      if (lofmt) then
        write (ounit,*)         name, rnk, shp
        write (ounit,'(g22.15)') y
      else
        write (ounit) name, rnk, shp
        write (ounit) y
      endif
    endif
    if(present (rank)) rank = rnk
    if(present (shap)) shap = shp
  end subroutine readfile2
!------------------------------------------------------------------------------
end module mo_1d_io
