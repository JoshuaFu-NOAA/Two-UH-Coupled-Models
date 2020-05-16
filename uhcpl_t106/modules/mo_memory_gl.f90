MODULE mo_memory_gl

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base, ONLY: create_list, delete_list, new_entry, get_entry,&
                            print_memory_table, print_memory_use, print_sinfo,&
                            get_info, memory_info 
  USE mo_netCDF,      ONLY: max_dim_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_gl ! construct the gl table
  PUBLIC :: destruct_gl  ! destruct  the gl table

  PUBLIC :: new_entry
  PUBLIC :: get_entry
  PUBLIC :: get_info

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  PUBLIC :: memory_info            ! meta data

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: q(:,:,:)
  REAL(dp), POINTER, PUBLIC :: x(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xt(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: lammp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: phimp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sigmp(:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: gl

CONTAINS

  SUBROUTINE construct_gl (lnlon, nglpx, lnlev, lntrac, lngl, &
                           nlon,         nlev,  ntrac,  ngl)

    INTEGER, INTENT (in) :: lnlon, nglpx, lnlev, lntrac, lngl
    INTEGER, INTENT (in) ::  nlon,         nlev,  ntrac,  ngl

    INTEGER :: nlp2
    INTEGER :: dim1(3), dim1p(3)
    INTEGER :: dim2(4), dim2p(4)
    INTEGER :: dim3(3), dim3p(3)
    CHARACTER (max_dim_name) :: dim1n(3), dim2n(4), dim3n(3)

    ! construct the gl table
    !
    ! all information specific to this table is set in this subroutine

    nlp2  = nlon  + 2

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (gl)

    ! assign pointers

    dim1p = (/  nglpx, lnlev, lngl  /)
    dim1  = (/  nlp2,   nlev,  ngl  /)
    dim1n = (/ "nlp2", "nlev","ngl "/)

    dim2p = (/ lnlon,   lnlev,   lntrac,  lngl    /)
    dim2  = (/  nlon,    nlev,    ntrac,   ngl    /)
    dim2n = (/ "nlon  ","nlev  ","nhtrac","ngl   "/)

    dim3p = (/ lnlon, lnlev, lngl  /)
    dim3  = (/  nlon,  nlev,  ngl  /)
    dim3n = (/ "nlon","nlev","ngl "/)

    CALL new_entry (gl, 'Q',  q, dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (gl, 'X',  x, dim1p, dim1, dimnames=dim1n, restart=.true.)
    IF (ntrac > 0) THEN
       CALL new_entry (gl, 'XT', xt, dim2p, dim2, dimnames=dim2n, restart=.true.)
    ELSE
       CALL new_entry (gl, 'XT', xt, dim2p, dim2)
    END IF
    CALL new_entry (gl, 'LAMMP',  lammp, dim3p, dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (gl, 'PHIMP',  phimp, dim3p, dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (gl, 'SIGMP',  sigmp, dim3p, dim3, dimnames=dim3n, restart=.true.)

  END SUBROUTINE construct_gl

  SUBROUTINE destruct_gl

    CALL delete_list (gl)

  END SUBROUTINE destruct_gl

END MODULE mo_memory_gl
