MODULE mo_memory_sp

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base, ONLY: create_list, delete_list, new_entry, get_entry,&
                            print_memory_table, print_memory_use, print_sinfo,&
                            get_info, memory_info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_sp ! construct the sp table
  PUBLIC :: destruct_sp  ! destruct  the sp table

  PUBLIC :: new_entry
  PUBLIC :: get_entry
  PUBLIC :: get_info

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  PUBLIC :: memory_info            ! meta data

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC  :: svo(:,:,:)
  REAL(dp), POINTER, PUBLIC  :: sd(:,:,:) 
  REAL(dp), POINTER, PUBLIC  :: stp(:,:,:)
  REAL(dp), POINTER, PUBLIC  :: su0(:,:) 

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: sp

CONTAINS

  SUBROUTINE construct_sp (lnlev, lnnp1, lnsp, nlev, nnp1, nsp)

    INTEGER, INTENT (in) :: lnlev, lnnp1, lnsp
    INTEGER, INTENT (in) :: nlev, nnp1, nsp

    ! construct the sp table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (sp)

    ! assign pointers

    CALL new_entry (sp,'SVO',svo, (/ lnlev, 2, lnsp /),   (/ nlev, 2, nsp /),   138)
    CALL new_entry (sp,'SD', sd,  (/ lnlev, 2, lnsp /),   (/ nlev, 2, nsp /),   155)
    CALL new_entry (sp,'STP',stp, (/ lnlev+1, 2, lnsp /), (/ nlev+1, 2, nsp /), 130)

    CALL new_entry (sp,'SU0',su0, (/ lnlev, lnnp1 /),     (/ nlev, nnp1 /))

  END SUBROUTINE construct_sp

  SUBROUTINE destruct_sp

    CALL delete_list (sp)

  END SUBROUTINE destruct_sp

END MODULE mo_memory_sp
