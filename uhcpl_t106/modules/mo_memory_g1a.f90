MODULE mo_memory_g1a

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base 
  USE mo_netCDF,      ONLY: max_dim_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g1a ! construct the g1a table
  PUBLIC :: destruct_g1a  ! destruct  the g1a table

  PUBLIC :: new_entry
  PUBLIC :: get_entry

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: vom1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: dm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: alpsm1(:,:)
  REAL(dp), POINTER, PUBLIC :: qm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xm1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xtm1(:,:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: g1a

CONTAINS

  SUBROUTINE construct_g1a (lnlon, nglpx, lnlev, lntrac, lngl, &
                             nlon,         nlev,  ntrac,  ngl)

    INTEGER, INTENT (in) :: lnlon, nglpx, lnlev, lntrac, lngl
    INTEGER, INTENT (in) ::  nlon,         nlev,  ntrac,  ngl

    INTEGER :: nlp2, nlevp1
    INTEGER :: dim1(3), dim1p(3)
    INTEGER :: dim2(2), dim2p(2)
    INTEGER :: dim3(4), dim3p(4)
    CHARACTER (max_dim_name) :: dim1n(3), dim2n(2), dim3n(4)

    ! construct the g1a table
    !
    ! all information specific to this table is set in this subroutine

    nlp2   = nlon  + 2
    nlevp1 = nlev  + 1

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (g1a)

    ! assign pointers

    dim1p = (/ nglpx, lnlev, lngl  /)
    dim1  = (/  nlp2,  nlev,  ngl  /)
    dim1n = (/ "nlp2","nlev","ngl "/)

    dim2p = (/ nglpx, lngl  /)
    dim2  = (/  nlp2,  ngl  /)
    dim2n = (/ "nlp2","ngl "/)

    dim3p = (/ lnlon,   lnlev,   lntrac,  lngl    /)
    dim3  = (/  nlon,    nlev,    ntrac,   ngl    /)
    dim3n = (/ "nlon  ","nlev  ","nhtrac","ngl   "/)

    CALL new_entry (g1a, 'VOM1',   vom1,   dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g1a, 'DM1',    dm1,    dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g1a, 'TM1',    tm1,    dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g1a, 'ALPSM1', alpsm1, dim2p, dim2, dimnames=dim2n, restart=.true.)
    CALL new_entry (g1a, 'QM1',    qm1,    dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g1a, 'XM1',    xm1,    dim1p, dim1, dimnames=dim1n, restart=.true.)
    IF (ntrac > 0) THEN
       CALL new_entry (g1a, 'XTM1',   xtm1,   dim3p, dim3, dimnames=dim3n, restart=.true.)
    ELSE
       CALL new_entry (g1a, 'XTM1',   xtm1,   dim3p, dim3)
    END IF

  END SUBROUTINE construct_g1a

  SUBROUTINE destruct_g1a

    CALL delete_list (g1a)

  END SUBROUTINE destruct_g1a

END MODULE mo_memory_g1a
