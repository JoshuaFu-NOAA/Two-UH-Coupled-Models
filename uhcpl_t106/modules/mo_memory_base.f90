MODULE mo_memory_base

  ! Definition of data type 'memory_type' holding pointers to
  ! fields and metadata required for I/O.
  !
  ! This module only defines the data type and basic operations.
  ! A specific instance should be declared in another module.

  USE mo_kind,        ONLY: dp 
  USE mo_doctor,      ONLY: nout
  USE mo_netCDF,      ONLY: IO_get_varindx
  USE mo_linked_list

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: gptr
  TYPE gptr                         
     REAL(dp), POINTER   :: x(:,:,:)  
  END TYPE gptr

  PUBLIC :: create_list         ! construct empty list 
  PUBLIC :: delete_list         ! destruct list

  PUBLIC :: print_memory_table  ! print list/memory information
  PUBLIC :: print_memory_use    ! print used memory
  PUBLIC :: print_sinfo         ! print short information
  
  PUBLIC :: ass_entry           ! assign a new list entry
  PUBLIC :: new_entry           ! create/allocate a new list entry
  PUBLIC :: get_entry           ! obtain reference to existing list entry

  PUBLIC :: get_info            ! obtain meta data of list entry
  PUBLIC :: set_info_restart    ! set restart flag

  PUBLIC :: memory_info         ! meta data structure

  INTERFACE ass_entry
     MODULE PROCEDURE ass_list_entry_4d ! assign a new list entry
     MODULE PROCEDURE ass_list_entry_3d 
     MODULE PROCEDURE ass_list_entry_2d 
     MODULE PROCEDURE ass_list_entry_1d 
  END INTERFACE

  INTERFACE new_entry
     MODULE PROCEDURE new_list_entry_4d ! create a new list entry
     MODULE PROCEDURE new_list_entry_3d 
     MODULE PROCEDURE new_list_entry_2d 
     MODULE PROCEDURE new_list_entry_1d 
  END INTERFACE
  
  INTERFACE get_entry
     MODULE PROCEDURE get_list_entry_4d ! obtain reference to a list entry
     MODULE PROCEDURE get_list_entry_3d
     MODULE PROCEDURE get_list_entry_2d
     MODULE PROCEDURE get_list_entry_1d
  END INTERFACE
  
CONTAINS

  SUBROUTINE create_list (this_list)

    TYPE (list) :: this_list

    CALL construct_list (this_list)

  END SUBROUTINE create_list

  SUBROUTINE delete_list (this_list)

    TYPE (list) :: this_list

    CALL destruct_list (this_list)

  END SUBROUTINE delete_list

  SUBROUTINE get_info (this_list, name, info)
    TYPE (list)       , INTENT(in)  :: this_list     ! list
    CHARACTER (*)     , INTENT(in)  :: name          ! name of variable
    TYPE (memory_info), INTENT(out) :: info          ! variable meta data

    TYPE (list_element), POINTER :: requested_list_element

    ! obtain info

    requested_list_element => find_list_element (this_list, name)
  
    info = requested_list_element%field%info

  END SUBROUTINE get_info

  SUBROUTINE set_info_restart (this_list, name, restart)
    TYPE (list)       , INTENT(in)  :: this_list     ! list
    CHARACTER (*)     , INTENT(in)  :: name          ! name of variable
    LOGICAL           , INTENT(in)  :: restart 

    TYPE (list_element), POINTER :: requested_list_element

    ! obtain info

    requested_list_element => find_list_element (this_list, name)

    ! set restart flag
  
    requested_list_element%field%info%restart = restart

  END SUBROUTINE set_info_restart

  SUBROUTINE new_list_entry_4d (this_list, name, ptr, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 4      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:,:,:,:)  ! reference to allocated field
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names

    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 4d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name

    ALLOCATE (new_list_element%field%ptr(ldims(1),ldims(2),ldims(3),ldims(4)))

    new_list_element%field%ptr = 0.

    this_list%memory_used = this_list%memory_used &
                            + 8*SIZE(new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = ldims(2)
    new_list_element%field%info%dim_3 = ldims(3)
    new_list_element%field%info%dim_4 = ldims(4)

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = gdims(2)
    new_list_element%field%info%gdim_3 = gdims(3)
    new_list_element%field%info%gdim_4 = gdims(4)

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE new_list_entry_4d

  SUBROUTINE new_list_entry_3d (this_list, name, ptr, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 3      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:,:,:)    ! reference to allocated field
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names

    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 3d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name

    ALLOCATE (new_list_element%field%ptr(ldims(1),ldims(2),ldims(3),1))

    new_list_element%field%ptr = 0.

    this_list%memory_used = this_list%memory_used &
                            + 8*SIZE(new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = ldims(2)
    new_list_element%field%info%dim_3 = ldims(3)
    new_list_element%field%info%dim_4 = 1

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = gdims(2)
    new_list_element%field%info%gdim_3 = gdims(3)
    new_list_element%field%info%gdim_4 = 1

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr(:,:,:,1)

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE new_list_entry_3d

  SUBROUTINE new_list_entry_2d (this_list, name, ptr, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 2      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:,:)      ! reference to allocated field
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names


    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 2d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name

    ALLOCATE (new_list_element%field%ptr(ldims(1),ldims(2),1,1))

    new_list_element%field%ptr = 0.

    this_list%memory_used = this_list%memory_used &
                          + 8*SIZE(new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = ldims(2)
    new_list_element%field%info%dim_3 = 1
    new_list_element%field%info%dim_4 = 1

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = gdims(2)
    new_list_element%field%info%gdim_3 = 1
    new_list_element%field%info%gdim_4 = 1

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr(:,:,1,1)

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE new_list_entry_2d

  SUBROUTINE new_list_entry_1d (this_list, name, ptr, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 1      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:)        ! reference to allocated field
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names

    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 1d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name
    ALLOCATE (new_list_element%field%ptr(ldims(1),1,1,1))

    new_list_element%field%ptr = 0.

    this_list%memory_used = this_list%memory_used &
                            + 8*SIZE(new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = 1
    new_list_element%field%info%dim_3 = 1
    new_list_element%field%info%dim_4 = 1

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = 1
    new_list_element%field%info%gdim_3 = 1
    new_list_element%field%info%gdim_4 = 1

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr(:,1,1,1)

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE new_list_entry_1d

  SUBROUTINE get_list_entry_4d (this_list, name, ptr)

    TYPE (list)               :: this_list    ! list
    CHARACTER (*), INTENT(in) :: name         ! name of variable
    REAL(dp)     , POINTER    :: ptr(:,:,:,:) ! reference to allocated field

    TYPE (list_element), POINTER :: requested_list_element

    ! obtain pointer to 4d-field

    requested_list_element => find_list_element (this_list, name)
  
    ptr => requested_list_element%field%ptr

  END SUBROUTINE get_list_entry_4d

  SUBROUTINE get_list_entry_3d (this_list, name, ptr)

    TYPE (list)               :: this_list    ! list
    CHARACTER (*), INTENT(in) :: name         ! name of variable
    REAL(dp)     , POINTER    :: ptr(:,:,:)   ! reference to allocated field

    TYPE (list_element), POINTER :: requested_list_element

    ! obtain pointer to 3d-field

    requested_list_element => find_list_element (this_list, name)

    ptr => requested_list_element%field%ptr(:,:,:,1)

  END SUBROUTINE get_list_entry_3d

  SUBROUTINE get_list_entry_2d (this_list, name, ptr)

    TYPE (list)       , INTENT(in) :: this_list    ! list
    CHARACTER (*)     , INTENT(in) :: name         ! name of variable
    REAL(dp)          , POINTER    :: ptr(:,:)   ! reference to allocated field

    TYPE (list_element), POINTER :: requested_list_element

    ! obtain pointer to 2d-field

    requested_list_element => find_list_element (this_list, name)

    ptr => requested_list_element%field%ptr(:,:,1,1)

  END SUBROUTINE get_list_entry_2d

  SUBROUTINE get_list_entry_1d (this_list, name, ptr)

    TYPE (list)  , INTENT(in) :: this_list    ! list
    CHARACTER (*), INTENT(in) :: name         ! name of variable
    REAL(dp)     , POINTER    :: ptr(:)       ! reference to allocated field

    TYPE (list_element), POINTER :: requested_list_element

    ! obtain pointer to 1d-field

    requested_list_element => find_list_element (this_list, name)

    ptr => requested_list_element%field%ptr(:,1,1,1)

  END SUBROUTINE get_list_entry_1d

  SUBROUTINE print_memory_use (this_list)

    TYPE (list)       , INTENT(in) :: this_list    ! list

    WRITE (nout,'(a,i10,a,i4,a)') &
         'Memory in use: ', this_list%memory_used, ' bytes in ', &
         this_list%list_elements, ' fields.'

  END SUBROUTINE print_memory_use

  SUBROUTINE print_memory_table (this_list)

    TYPE (list),  INTENT(in) :: this_list ! list

    ! print current memory table 
    
    WRITE (nout,'(/,/,a,/)') &
         'Status of base memory:'     

    CALL print_linked_list (this_list)
    
  END SUBROUTINE print_memory_table

  SUBROUTINE print_sinfo (this_list)

    TYPE (list),  INTENT(in) :: this_list ! list

    ! print current stat table 
    
    WRITE (nout,'(/,/,a,/)') &
         'Statistic of base memory:'     

    CALL print_sinfo_list (this_list)
    
  END SUBROUTINE print_sinfo

  SUBROUTINE ass_list_entry_4d (this_list, name, ptr, ass_list, ass_name, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 4      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:,:,:,:)  ! reference to allocated field
    TYPE (list)       , INTENT(in)    :: ass_list      ! assign list
    CHARACTER (*)     , INTENT(in)    :: ass_name      ! name of assign variable
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names


    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 4d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name

    CALL get_entry (ass_list, ass_name, new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = ldims(2)
    new_list_element%field%info%dim_3 = ldims(3)
    new_list_element%field%info%dim_4 = ldims(4)

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = gdims(2)
    new_list_element%field%info%gdim_3 = gdims(3)
    new_list_element%field%info%gdim_4 = gdims(4)

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr(:,:,:,:)

    new_list_element%field%info%assign = .TRUE.

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE ass_list_entry_4d


  SUBROUTINE ass_list_entry_3d (this_list, name, ptr, ass_list, ass_name, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 3      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:,:,:)    ! reference to allocated field
    TYPE (list)       , INTENT(in)    :: ass_list      ! assign list
    CHARACTER (*)     , INTENT(in)    :: ass_name      ! name of assign variable
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names


    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 3d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name

    CALL get_entry (ass_list, ass_name, new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = ldims(2)
    new_list_element%field%info%dim_3 = ldims(3)
    new_list_element%field%info%dim_4 = 1

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = gdims(2)
    new_list_element%field%info%gdim_3 = gdims(3)
    new_list_element%field%info%gdim_4 = 1

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr(:,:,:,1)

    new_list_element%field%info%assign = .TRUE.

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE ass_list_entry_3d

  SUBROUTINE ass_list_entry_2d (this_list, name, ptr, ass_list, ass_name, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 2      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:,:)      ! reference to allocated field
    TYPE (list)       , INTENT(in)    :: ass_list      ! assign list
    CHARACTER (*)     , INTENT(in)    :: ass_name      ! name of assign variable
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names


    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 2d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name

    CALL get_entry (ass_list, ass_name, new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = ldims(2)
    new_list_element%field%info%dim_3 = 1
    new_list_element%field%info%dim_4 = 1

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = gdims(2)
    new_list_element%field%info%gdim_3 = 1
    new_list_element%field%info%gdim_4 = 1

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr(:,:,1,1)

    new_list_element%field%info%assign = .TRUE.

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE ass_list_entry_2d

  SUBROUTINE ass_list_entry_1d (this_list, name, ptr, ass_list, ass_name, ldims, gdims, &
       gribcode, gribtable, outint, accumulate, restart, dimnames)

    INTEGER           , PARAMETER     :: dims = 1      ! dimensions
    TYPE (list)       , INTENT(inout) :: this_list     ! list
    CHARACTER (*)     , INTENT(in)    :: name          ! name of variable
    REAL(dp)          , POINTER       :: ptr(:)        ! reference to allocated field
    TYPE (list)       , INTENT(in)    :: ass_list      ! assign list
    CHARACTER (*)     , INTENT(in)    :: ass_name      ! name of assign variable
    INTEGER           , INTENT(in)    :: ldims(dims)   ! shape of array to allocate
    INTEGER           , INTENT(in)    :: gdims(dims)   ! global size of field

    CHARACTER (*), OPTIONAL  , INTENT(in)  :: dimnames(dims)   ! dimension names


    INTEGER,  OPTIONAL, INTENT(in)    :: gribcode      ! gribcode number
    INTEGER,  OPTIONAL, INTENT(in)    :: gribtable     ! gribcode table number

    INTEGER,  OPTIONAL, INTENT(in)    :: outint        ! output interval (time steps)
    LOGICAL,  OPTIONAL, INTENT(in)    :: accumulate    ! accumulation flag
    LOGICAL,  OPTIONAL, INTENT(in)    :: restart       ! restart file flag

    INTEGER :: idim

    ! create (allocate) a new table entry
    ! optionally obtain pointer to 1d-field
    ! optionally overwrite default meta data 

    TYPE (list_element), POINTER :: new_list_element

    ! add list entry

    CALL add_list_element (this_list, new_list_element)

    ! and set meta data
      
    new_list_element%field%info%name = name

    CALL get_entry (ass_list, ass_name, new_list_element%field%ptr)

    new_list_element%field%info%dim_1 = ldims(1)
    new_list_element%field%info%dim_2 = 1
    new_list_element%field%info%dim_3 = 1
    new_list_element%field%info%dim_4 = 1

    new_list_element%field%info%gdim_1 = gdims(1)
    new_list_element%field%info%gdim_2 = 1
    new_list_element%field%info%gdim_3 = 1
    new_list_element%field%info%gdim_4 = 1

    new_list_element%field%info%ndim = dims

    ptr => new_list_element%field%ptr(:,1,1,1)

    new_list_element%field%info%assign = .TRUE.

    ! pass optional arguments
          
    IF (PRESENT(dimnames)) THEN
       DO idim = 1, dims
          new_list_element%field%info%IO_var_indx(idim) = IO_get_varindx(dimnames(idim))
       END DO
    END IF

    IF (PRESENT(gribtable))  new_list_element%field%info%gribtable  = gribtable 
    IF (PRESENT(gribcode))   new_list_element%field%info%gribcode   = gribcode
    IF (PRESENT(outint))     new_list_element%field%info%outint     = outint
    IF (PRESENT(accumulate)) new_list_element%field%info%accumulate = accumulate
    IF (PRESENT(restart))    new_list_element%field%info%restart    = restart

  END SUBROUTINE ass_list_entry_1d

END MODULE mo_memory_base


