MODULE mo_memory_ls

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base, ONLY: create_list, delete_list, new_entry, get_entry,   &
                            print_memory_table, print_memory_use, print_sinfo,&
                            get_info, memory_info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_ls ! construct the ls table
  PUBLIC :: destruct_ls  ! destruct  the ls table

  PUBLIC :: new_entry
  PUBLIC :: get_entry
  PUBLIC :: get_info

  PUBLIC :: print_memory_table     ! print information on ls table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  PUBLIC :: memory_info            ! meta data

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC  :: lvo (:,:,:)
  REAL(dp), POINTER, PUBLIC  :: ld  (:,:,:) 
  REAL(dp), POINTER, PUBLIC  :: ltp (:,:,:)
  REAL(dp), POINTER, PUBLIC  :: lu0 (:,:) 

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: ls

CONTAINS

  SUBROUTINE construct_ls (lnlev, lnlevp1, lnnp1, lnsp, nlev, nnp1, nsp)

    INTEGER, INTENT (in) :: lnlev, lnlevp1, lnnp1, lnsp ! local array bounds
    INTEGER, INTENT (in) :: nlev, nnp1, nsp             ! global array bounds

    ! construct the ls table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (ls)

    ! assign pointers

    CALL new_entry (ls,'SVO',lvo, (/ lnlev,   2, lnsp /), (/ nlev, 2, nsp /),   138)
    CALL new_entry (ls,'SD', ld,  (/ lnlev,   2, lnsp /), (/ nlev, 2, nsp /),   155)
    CALL new_entry (ls,'STP',ltp, (/ lnlevp1, 2, lnsp /), (/ nlev+1, 2, nsp /), 130)

    CALL new_entry (ls,'SU0',lu0, (/ lnlev,   lnnp1 /),   (/ nlev, nnp1 /))

  END SUBROUTINE construct_ls

  SUBROUTINE destruct_ls

    CALL delete_list (ls)

  END SUBROUTINE destruct_ls

END MODULE mo_memory_ls
