MODULE mo_memory_g2b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g2b ! construct the g2b table
  PUBLIC :: destruct_g2b  ! destruct  the g2b table

  PUBLIC :: new_entry
  PUBLIC :: get_entry

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: uf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: vf(:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: g2b

CONTAINS

  SUBROUTINE construct_g2b (lnlon, nglpx, lnlev, lngl, nlon, nlev, ngl)

    INTEGER, INTENT (in) :: lnlon, nglpx, lnlev, lngl
    INTEGER, INTENT (in) ::  nlon,         nlev,  ngl

    INTEGER :: nlp2
    INTEGER :: dim1(3), dim1p(3)

    ! construct the g2b table
    !
    ! all information specific to this table is set in this subroutine

    nlp2   = nlon  + 2

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (g2b)

    ! assign pointers

    dim1p = (/ nglpx, lnlev, lngl /)
    dim1  = (/  nlp2,  nlev,  ngl /)

    CALL new_entry (g2b, 'UF', uf, dim1p, dim1)
    CALL new_entry (g2b, 'VF', vf, dim1p, dim1)

  END SUBROUTINE construct_g2b

  SUBROUTINE destruct_g2b

    CALL delete_list (g2b)

  END SUBROUTINE destruct_g2b

END MODULE mo_memory_g2b
