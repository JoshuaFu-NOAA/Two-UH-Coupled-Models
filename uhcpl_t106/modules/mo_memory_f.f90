MODULE mo_memory_f

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base 
  USE mo_netCDF,      ONLY: max_dim_name

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: construct_f ! construct the f table
  PUBLIC :: destruct_f  ! destruct  the f table

  PUBLIC :: new_entry
  PUBLIC :: get_entry

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  ! declaration of predefined fields within this module 

  ! used in inverse Legendre transformation

  ! pointers for *f1* space.
  REAL(dp), POINTER, PUBLIC  :: fsvo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: favo(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsu(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fau(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsv(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fav(:,:,:,:)
  ! pointers for *f3* space.
  REAL(dp), POINTER, PUBLIC  :: fsd(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fad(:,:,:,:)
  ! pointers for *f4* space.
  REAL(dp), POINTER, PUBLIC  :: fstp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fatp(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fstpm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fatpm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC  :: fsu0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fau0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fsdu0(:,:)
  REAL(dp), POINTER, PUBLIC  :: fadu0(:,:)

  ! used in direct Legendre transformation

  ! pointers for *f5* space.
  REAL(dp), POINTER, PUBLIC :: fszl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fazl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fszm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fazm(:,:,:,:)
  ! pointers for *f7* space.
  REAL(dp), POINTER, PUBLIC :: fsdl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fadl(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsdm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fadm(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsr(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: far(:,:,:,:)
  ! pointers for *f8* space.
  REAL(dp), POINTER, PUBLIC :: fstp1(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fatp1(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: fsul(:,:)
  REAL(dp), POINTER, PUBLIC :: faul(:,:)

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: f

CONTAINS

  SUBROUTINE construct_f (lnlev, lnlevp1, lnmp1, lnhgl, nlev, nmp1, nhgl)

    INTEGER, INTENT (in) :: lnlev, lnlevp1, lnmp1, lnhgl
    INTEGER, INTENT (in) ::  nlev,           nmp1,  nhgl
    INTEGER :: dim1(4), dim1p(4)
    INTEGER :: dim2(4), dim2p(4)
    INTEGER :: dim3(2), dim3p(2)
    CHARACTER (max_dim_name) :: dim1n(4), dim2n(4), dim3n(2)


    ! construct the f table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (f)

    ! assign pointers

    dim1p = (/ lnlev,   2,   lnmp1, lnhgl /)
    dim1  = (/  nlev,   2,    nmp1,  nhgl /)
    dim1n = (/ "nlev","n2  ","nmp1","nhgl"/)

    dim2p = (/ lnlevp1,   2,     lnmp1,   lnhgl   /)
    dim2  = (/  nlev+1,   2,      nmp1,    nhgl   /)
    dim2n = (/ "nlevp1","n2    ","nmp1  ","nhgl  "/)

    dim3p = (/ lnlev, lnhgl /)
    dim3  = (/  nlev,  nhgl /)
    dim3n = (/ "nlev","nhgl"/)

    ! Arrays used by inverse transform

    CALL new_entry (f, 'FSVO',  fsvo,  dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FAVO',  favo,  dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FSU',   fsu,   dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FAU',   fau,   dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FSV',   fsv,   dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FAV',   fav,   dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FSD',   fsd,   dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FAD',   fad,   dim1p, dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (f, 'FSTP',  fstp,  dim2p, dim2, dimnames=dim2n, restart=.true.)
    CALL new_entry (f, 'FATP',  fatp,  dim2p, dim2, dimnames=dim2n, restart=.true.)
    CALL new_entry (f, 'FSTPM', fstpm, dim2p, dim2, dimnames=dim2n, restart=.true.)
    CALL new_entry (f, 'FATPM', fatpm, dim2p, dim2, dimnames=dim2n, restart=.true.)
    CALL new_entry (f, 'FSU0',  fsu0,  dim3p, dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (f, 'FAU0',  fau0,  dim3p, dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (f, 'FSDU0', fsdu0, dim3p, dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (f, 'FADU0', fadu0, dim3p, dim3, dimnames=dim3n, restart=.true.)

    ! Arrays used by direct transform (not yet in memory buffer)

    ALLOCATE (fsdl (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fadm (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fsr  (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fszl (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fazm (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fstp1(lnlevp1,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fsul (lnlev            ,lnhgl))
    ALLOCATE (fadl (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fsdm (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (far  (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fazl (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fszm (lnlev  ,2 ,lnmp1 ,lnhgl))
    ALLOCATE (fatp1(lnlevp1,2 ,lnmp1 ,lnhgl))
    ALLOCATE (faul (lnlev            ,lnhgl))

  END SUBROUTINE construct_f

  SUBROUTINE destruct_f

    CALL delete_list (f)

  END SUBROUTINE destruct_f

END MODULE mo_memory_f
