MODULE mo_memory_g2a

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base 
  USE mo_netCDF,      ONLY: max_dim_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g2a ! construct the g2a table
  PUBLIC :: destruct_g2a  ! destruct  the g2a table

  PUBLIC :: new_entry
  PUBLIC :: get_entry

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: um1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: vm1(:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: g2a

CONTAINS

  SUBROUTINE construct_g2a (lnlon, nglpx, lnlev, lngl, nlon, nlev, ngl)

    INTEGER, INTENT (in) :: lnlon, nglpx, lnlev, lngl
    INTEGER, INTENT (in) ::  nlon,         nlev,  ngl

    INTEGER :: nlp2
    INTEGER :: dim1(3), dim1p(3)
    CHARACTER (max_dim_name) :: dim1n(3)

    ! construct the g2a table
    !
    ! all information specific to this table is set in this subroutine

    nlp2   = nlon  + 2

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (g2a)

    ! assign pointers

    dim1p = (/ nglpx, lnlev, lngl  /)
    dim1  = (/  nlp2,  nlev,  ngl  /)
    dim1n = (/ "nlp2","nlev","ngl "/)

    CALL new_entry (g2a, 'UM1' ,um1, dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g2a, 'VM1' ,vm1, dim1p, dim1, dimnames=dim1n, restart=.true.)

  END SUBROUTINE construct_g2a

  SUBROUTINE destruct_g2a

    CALL delete_list (g2a)

  END SUBROUTINE destruct_g2a

END MODULE mo_memory_g2a
