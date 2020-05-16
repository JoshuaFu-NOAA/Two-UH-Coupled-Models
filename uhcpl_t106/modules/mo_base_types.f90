MODULE mo_base_types
  TYPE base_pointer
     REAL, POINTER :: vp(:)
     CHARACTER (8) :: cp
     INTEGER :: np
  END TYPE base_pointer

  TYPE pointer_set
     TYPE (base_pointer), POINTER :: set(:)
     INTEGER :: nsets
  END TYPE pointer_set

END MODULE mo_base_types
