MODULE mo_memory_g1b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g1b ! construct the g1b table
  PUBLIC :: destruct_g1b  ! destruct  the g1b table

  PUBLIC :: new_entry
  PUBLIC :: get_entry

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: vof(:,:,:)
  REAL(dp), POINTER, PUBLIC :: df(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: alpsf(:,:)
  REAL(dp), POINTER, PUBLIC :: qf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xtf(:,:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: g1b

CONTAINS

  SUBROUTINE construct_g1b (lnlon, nglpx, lnlev, lntrac, lngl, &
                             nlon,         nlev,  ntrac,  ngl)

    INTEGER, INTENT (in) :: lnlon, nglpx, lnlev, lntrac, lngl
    INTEGER, INTENT (in) ::  nlon,         nlev,  ntrac,  ngl

    INTEGER :: nlp2, nlevp1
    INTEGER :: dim1(3), dim1p(3)
    INTEGER :: dim2(2), dim2p(2)
    INTEGER :: dim3(4), dim3p(4)

    ! construct the g1b table
    !
    ! all information specific to this table is set in this subroutine

    nlp2   = nlon  + 2
    nlevp1 = nlev  + 1

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (g1b)

    ! assign pointers

    dim1p = (/ nglpx, lnlev, lngl /)
    dim1  = (/  nlp2,  nlev,  ngl /)

    dim2p = (/ nglpx, lngl /)
    dim2  = (/  nlp2,  ngl /)

    dim3p = (/ lnlon, lnlev,  lntrac, lngl /)
    dim3  = (/  nlon,  nlev,   ntrac,  ngl /)

    CALL new_entry (g1b, 'VOF',   vof,   dim1p, dim1)
    CALL new_entry (g1b, 'DF',    df,    dim1p, dim1)
    CALL new_entry (g1b, 'TF',    tf,    dim1p, dim1)
    CALL new_entry (g1b, 'ALPSF', alpsf, dim2p, dim2)
    CALL new_entry (g1b, 'QF',    qf,    dim1p, dim1)
    CALL new_entry (g1b, 'XF',    xf,    dim1p, dim1)
    CALL new_entry (g1b, 'XTF',   xtf,   dim3p, dim3)

  END SUBROUTINE construct_g1b

  SUBROUTINE destruct_g1b

    CALL delete_list (g1b)

  END SUBROUTINE destruct_g1b

END MODULE mo_memory_g1b
