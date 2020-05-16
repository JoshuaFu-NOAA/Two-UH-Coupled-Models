MODULE mo_linked_list

  USE mo_kind, ONLY : dp, i8
  USE mo_doctor, ONLY: nout

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: memory_type         ! data type to hold field list entry
  PUBLIC :: memory_info         ! meta data type
  PUBLIC :: empty_info          ! empty list entry

  PUBLIC :: list                ! anchor for a whole list
  PUBLIC :: list_element        

  PUBLIC :: construct_list      
  PUBLIC :: destruct_list
  
  PUBLIC :: add_list_element
  PUBLIC :: find_list_element

  PUBLIC :: print_linked_list
  PUBLIC :: print_sinfo_list

  ! The following definitions are the base for the linked list and may be
  ! replaced by other definitions or extended. Don't forget to change the
  ! default constructor according to the changes in TYPE memory info.

  TYPE memory_info                      ! meta data type
     CHARACTER (32)    :: name          ! variable name 
     INTEGER           :: dim_1         ! local dimensions of variable
     INTEGER           :: dim_2
     INTEGER           :: dim_3
     INTEGER           :: dim_4
     INTEGER           :: gdim_1        ! global dimaensions of variable
     INTEGER           :: gdim_2        ! global dimaensions of variable
     INTEGER           :: gdim_3        ! global dimaensions of variable
     INTEGER           :: gdim_4        ! global dimaensions of variable
     INTEGER           :: ndim
     INTEGER           :: gribtable     ! gribcode table number
     INTEGER           :: gribcode      ! gribcode number
     INTEGER           :: IO_var_indx(4)  ! NETCDF id for internal use
     INTEGER           :: IO_var_id     ! NETCDF id for internal use
     CHARACTER (128)   :: IO_name       ! name specifier for NETCDF file
     CHARACTER (32)    :: IO_unit       ! unit specifier for NETCDF file
     INTEGER           :: outint        ! output interval (in time steps)
     LOGICAL           :: accumulate    ! accumulation flag
     LOGICAL           :: assign        ! assign flag
     LOGICAL           :: restart       ! restart file flag
  END TYPE memory_info
  
  TYPE memory_type                         ! linked list entry type
     REAL(dp), POINTER   :: ptr (:,:,:,:)  ! pointer to 3D-field
     TYPE (memory_info)  :: info           ! meta data for this entry
  END TYPE memory_type

  ! default (empty) meta data entry

  TYPE (memory_info), PARAMETER :: empty_info = &
       memory_info('        '      , &
                  0,   1,   1,    1, &     ! 0 size, but allows 1/2/3d assign  
                  0,   0,   0,    0, &     ! global dimensions unknown 
                  0,                 &     !
                  128              , &     ! gribcode table number
                  0                , &     ! gribcode number
                  (/ 0, 0, 0, 0 /) , &     ! NETCDF indx for internal use
                  0,                 &     ! NETCDF id for internal use
                  'undefined'      , &     ! name specifier for NETCDF file
                  'undefined'      , &     ! unit specifier for NETCDF file
                  0                , &     ! output interval
                  .FALSE.          , &     ! accumulation flag
                  .FALSE.          , &     ! assign flag
                  .FALSE.)                 ! restart file flag

  TYPE list_element
     TYPE (memory_type) :: field
     TYPE(list_element), POINTER :: next_list_element
  END TYPE list_element

  TYPE list
     TYPE (list_element), POINTER :: first_list_element 
     INTEGER (i8) :: memory_used 
     INTEGER :: list_elements
  END TYPE list

CONTAINS

  SUBROUTINE construct_list (this_list)

    TYPE (list) :: this_list

    NULLIFY (this_list%first_list_element)
    this_list%memory_used = 0
    this_list%list_elements = 0
    
  END SUBROUTINE construct_list

  SUBROUTINE destruct_list (this_list)

    TYPE (list) :: this_list

    CALL destruct_list_element (this_list, this_list%first_list_element)

    nullify (this_list%first_list_element)
    IF (this_list%memory_used /= 0) THEN
       WRITE (nout,*) &
            'List destructor didn''t work proper (memory counter) ...'
       STOP
    ENDIF

    IF (this_list%list_elements /= 0) THEN
       WRITE (nout,*) &
            'List destructor didn''t work proper (element counter) ...'
       STOP
    ENDIF

  END SUBROUTINE destruct_list

  RECURSIVE SUBROUTINE destruct_list_element (this_list, this_list_element)

    TYPE (list) :: this_list
    TYPE (list_element), POINTER :: this_list_element

    IF (ASSOCIATED(this_list_element)) THEN

       CALL destruct_list_element (this_list, &
            this_list_element%next_list_element)

       ! 8 as constant has to be adjusted with information from mo_machine
       ! the variable to be used is mp_real8

       this_list%memory_used = this_list%memory_used &
            -8*SIZE(this_list_element%field%ptr)

       DEALLOCATE (this_list_element%field%ptr)
    
       this_list%list_elements = this_list%list_elements-1
    
       DEALLOCATE (this_list_element)
 
    ENDIF
  
  END SUBROUTINE destruct_list_element

  SUBROUTINE create_list_element (this_list, current_list_element)

    TYPE (list) :: this_list
    TYPE (list_element), POINTER :: current_list_element
    
    ALLOCATE (current_list_element)

    this_list%list_elements = this_list%list_elements+1
    
    IF (.NOT. ASSOCIATED (current_list_element)) THEN
       WRITE (nout,*) 'Cannot add element to linked list ...'
       STOP
    ENDIF

    NULLIFY (current_list_element%next_list_element)
    NULLIFY (current_list_element%field%ptr)
    current_list_element%field%info = empty_info

  END SUBROUTINE create_list_element

  SUBROUTINE add_list_element (this_list, new_list_element)

    TYPE (list) :: this_list
    TYPE (list_element), POINTER :: new_list_element

    TYPE (list_element), POINTER :: current_list_element

   IF (.NOT. ASSOCIATED (this_list%first_list_element)) THEN
      CALL create_list_element (this_list, this_list%first_list_element)
      new_list_element => this_list%first_list_element
      RETURN
   ENDIF   

   current_list_element => this_list%first_list_element
   DO WHILE (ASSOCIATED(current_list_element%next_list_element)) 
      current_list_element => current_list_element%next_list_element
   ENDDO

   CALL create_list_element (this_list, new_list_element)
   current_list_element%next_list_element => new_list_element 

  END SUBROUTINE add_list_element

  SUBROUTINE delete_list_element (this_list, delete_this_list_element)

    TYPE (list) :: this_list
    TYPE (list_element), POINTER :: delete_this_list_element

    TYPE (list_element), POINTER :: current_list_element

    IF (ASSOCIATED(delete_this_list_element, &
         this_list%first_list_element)) THEN
       this_list%first_list_element &
            => delete_this_list_element%next_list_element
    ELSE
       current_list_element => this_list%first_list_element
       DO WHILE ((ASSOCIATED(current_list_element)) &
            .AND. (.NOT. ASSOCIATED(current_list_element%next_list_element, &
                                    delete_this_list_element)))
          current_list_element => current_list_element%next_list_element
       ENDDO
       IF (.NOT. ASSOCIATED(current_list_element)) THEN
          WRITE (nout,*) 'Cannot find element to be deleted ...'
          RETURN
       ENDIF
       current_list_element%next_list_element &
            => current_list_element%next_list_element%next_list_element
    ENDIF

    this_list%memory_used = this_list%memory_used &
         -8*SIZE(delete_this_list_element%field%ptr)

    DEALLOCATE (delete_this_list_element%field%ptr)

    this_list%list_elements = this_list%list_elements-1

    DEALLOCATE (delete_this_list_element)

  END SUBROUTINE delete_list_element

  ! Should be overloaded to be able to search for the different information 
  ! In the proposed structure for the linked list, in the example only
  ! A character string is used so it is straight forward only one find

  FUNCTION find_list_element (this_list, name) RESULT (this_list_element)

    TYPE (list) :: this_list
    CHARACTER (*), INTENT(in) :: name

    TYPE (list_element), POINTER :: this_list_element

    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))
       IF (name == this_list_element%field%info%name) THEN
          RETURN
       ENDIF
       this_list_element => this_list_element%next_list_element
    ENDDO

    WRITE (nout,*) 'find_list_element: element ', name, ' not available ...'

    NULLIFY (this_list_element)

  END FUNCTION find_list_element

  SUBROUTINE print_linked_list (this_list)
    
    TYPE (list) :: this_list
    TYPE (list_element), POINTER :: this_list_element

    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))

       ! print ....
       IF (this_list_element%field%info%name /= '') THEN

          WRITE (nout,'(a,a)')       &
               'Table entry name      : ', &
               TRIM(this_list_element%field%info%name)
#ifdef DEBUG
          WRITE (nout,'(a,i10)') 'Address of data field : ', &
               LOC(this_list_element%field%ptr)
#endif
          IF (ASSOCIATED(this_list_element%field%ptr)) THEN
             WRITE (nout,'(a)')      &
                  'Pointer status        : in use.'

             IF (SIZE(this_list_element%field%ptr,4) == 1 & 
                  .AND. SIZE(this_list_element%field%ptr,3) == 1 &
                  .AND. SIZE(this_list_element%field%ptr,2) == 1) THEN
                WRITE (nout,'(a,1(i4,a))') &
                     'Local field dimensions      : (',  &
                     SIZE(this_list_element%field%ptr,1), ')' 
             ELSE IF (SIZE(this_list_element%field%ptr,4) == 1 &
                  .AND. SIZE(this_list_element%field%ptr,3) == 1) THEN
                WRITE (nout,'(a,2(i4,a))') &
                     'Local field dimensions      : (',  &
                     SIZE(this_list_element%field%ptr,1), ',', &
                     SIZE(this_list_element%field%ptr,2), ')' 
             ELSE IF (SIZE(this_list_element%field%ptr,4) == 1) THEN
                WRITE (nout,'(a,3(i4,a))') &
                     'Local field dimensions      : (',  &
                     SIZE(this_list_element%field%ptr,1), ',', &
                     SIZE(this_list_element%field%ptr,2), ',', &
                     SIZE(this_list_element%field%ptr,3), ')'       
             ELSE
                WRITE (nout,'(a,4(i4,a))') &
                     'Local field dimensions      : (',  &
                     SIZE(this_list_element%field%ptr,1), ',', &
                     SIZE(this_list_element%field%ptr,2), ',', &
                     SIZE(this_list_element%field%ptr,3), ',', &
                     SIZE(this_list_element%field%ptr,4), ')'       
             ENDIF
          ELSE
             WRITE (nout,'(a)')      &
                  'Pointer status       : not in use.'
          ENDIF
          
          IF (this_list_element%field%info%gdim_4 /= 0 & 
               .AND. this_list_element%field%info%gdim_3 /= 0 &
               .AND. this_list_element%field%info%gdim_2 /= 0) THEN
             WRITE (nout,'(a,1(i4,a))') &
                  'Global field dimensions      : (',  &
                  this_list_element%field%info%gdim_1, ')' 
          ELSE IF (this_list_element%field%info%gdim_4 /= 0 &
               .AND. this_list_element%field%info%gdim_3 /= 0) THEN
             WRITE (nout,'(a,2(i4,a))') &
                  'Global field dimensions      : (',  &
                  this_list_element%field%info%gdim_1, ',', &
                  this_list_element%field%info%gdim_2, ')' 
          ELSE IF (this_list_element%field%info%gdim_4 /= 0) THEN
             WRITE (nout,'(a,3(i4,a))') &
                  'Global field dimensions      : (',  &
                  this_list_element%field%info%gdim_1, ',', &
                  this_list_element%field%info%gdim_2, ',', &
                  this_list_element%field%info%gdim_3, ')'       
          ELSE
             WRITE (nout,'(a,4(i4,a))') &
                  'Global field dimensions      : (',  &
                  this_list_element%field%info%gdim_1, ',', &
                  this_list_element%field%info%gdim_2, ',', &
                  this_list_element%field%info%gdim_3, ',', &
                  this_list_element%field%info%gdim_4, ')'       
          ENDIF

          WRITE (nout,'(a,i3,/,a,i3)') &
               'Assigned GRIB table   : ', &
               this_list_element%field%info%gribtable, &
               '         GRIB code    : ', &
               this_list_element%field%info%gribcode

          WRITE (nout,'(a,i6,/,a,a,/,a,a)') &
               'IO id                 : ', &
               this_list_element%field%info%IO_var_id, &
               '   name               : ', &
               TRIM(this_list_element%field%info%IO_name), &
               '   unit               : ', &
               TRIM(this_list_element%field%info%IO_unit)

          WRITE (nout,'(a,i4,a)')      &
               'Output intervall      : ', &
               this_list_element%field%info%outint, ' timesteps'
          IF (this_list_element%field%info%accumulate) THEN
             WRITE (nout,'(a)')      &
                  'Accumulation          : on.'
          ELSE
             WRITE (nout,'(a)')      &
                  'Accumulation          : off.'
          ENDIF
          IF (this_list_element%field%info%restart) THEN
             WRITE (nout,'(a)')      &
                  'Restart table         : added.'
          ELSE
             WRITE (nout,'(a)')      &
                  'Restart table         : unused.'
          ENDIF
          WRITE (nout,'(/)')
       ENDIF
       
       ! select next element in linked list 

       this_list_element => this_list_element%next_list_element
    ENDDO

  END SUBROUTINE print_linked_list

  SUBROUTINE print_sinfo_list (this_list)
    
    TYPE (list) :: this_list
    TYPE (list_element), POINTER :: this_list_element
    CHARACTER (80) :: cout

!                  123456789+123456789+123456789+123456789+123456789+123456789+
    WRITE(nout,*) '   Name   Local dimension   Tab Code Outint Accu  Restart'
    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))

       ! print ....
       IF (this_list_element%field%info%name /= '') THEN

          WRITE(cout(1:10),'(a10)') TRIM(this_list_element%field%info%name)

          IF (ASSOCIATED(this_list_element%field%ptr)) THEN
             IF (SIZE(this_list_element%field%ptr,4) == 1 & 
                  .AND. SIZE(this_list_element%field%ptr,3) == 1 &
                  .AND. SIZE(this_list_element%field%ptr,2) == 1) THEN
                WRITE (cout(12:30),'(a,1(i3,a))') &
                     ' (',  &
                     SIZE(this_list_element%field%ptr,1), ')' 
             ELSE IF (SIZE(this_list_element%field%ptr,4) == 1 &
                  .AND. SIZE(this_list_element%field%ptr,3) == 1) THEN
                WRITE (cout(12:30),'(a,2(i3,a))') &
                     ' (',  &
                     SIZE(this_list_element%field%ptr,1), ',', &
                     SIZE(this_list_element%field%ptr,2), ')' 
             ELSE IF (SIZE(this_list_element%field%ptr,4) == 1) THEN
                WRITE (cout(12:30),'(a,3(i3,a))') &
                     ' (',  &
                     SIZE(this_list_element%field%ptr,1), ',', &
                     SIZE(this_list_element%field%ptr,2), ',', &
                     SIZE(this_list_element%field%ptr,3), ')'       
             ELSE
                WRITE (cout(12:30),'(a,4(i3,a))') &
                     ' (',  &
                     SIZE(this_list_element%field%ptr,1), ',', &
                     SIZE(this_list_element%field%ptr,2), ',', &
                     SIZE(this_list_element%field%ptr,3), ',', &
                     SIZE(this_list_element%field%ptr,4), ')'       
             ENDIF
          ELSE
             WRITE (cout(12:30),'(a)')      &
                  '    not in use '
          ENDIF
          WRITE (cout(31:33),'(i3)') this_list_element%field%info%gribtable
          WRITE (cout(34:38),'(i5)') this_list_element%field%info%gribcode

          WRITE (cout(40:43),'(i4)') this_list_element%field%info%outint
          IF (this_list_element%field%info%accumulate) THEN
             WRITE (cout(47:52),'(a)') ' on  '
          ELSE
             WRITE (cout(47:52),'(a)') ' off '
          ENDIF
          IF (this_list_element%field%info%restart) THEN
             WRITE (cout(53:60),'(a)') ' added '
          ELSE
             WRITE (cout(53:60),'(a)') ' unused'
          ENDIF
       ENDIF

       WRITE(nout,'(a)') cout
       
       ! select next element in linked list 

       this_list_element => this_list_element%next_list_element
    ENDDO

  END SUBROUTINE print_sinfo_list

END MODULE mo_linked_list
