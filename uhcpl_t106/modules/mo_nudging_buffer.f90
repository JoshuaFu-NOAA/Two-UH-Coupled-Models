MODULE mo_nudging_buffer
  REAL, ALLOCATABLE :: sdobs0(:,:,:) 
  REAL, ALLOCATABLE :: svobs0(:,:,:)
  REAL, ALLOCATABLE :: stobs0(:,:,:)

  REAL, ALLOCATABLE :: sdobs1(:,:,:) 
  REAL, ALLOCATABLE :: svobs1(:,:,:)
  REAL, ALLOCATABLE :: stobs1(:,:,:)

  REAL, ALLOCATABLE :: sdobs2(:,:,:)
  REAL, ALLOCATABLE :: svobs2(:,:,:)
  REAL, ALLOCATABLE :: stobs2(:,:,:)

  REAL, ALLOCATABLE :: sdobs3(:,:,:) 
  REAL, ALLOCATABLE :: svobs3(:,:,:)
  REAL, ALLOCATABLE :: stobs3(:,:,:)

  REAL, ALLOCATABLE :: sdobs(:,:,:)
  REAL, ALLOCATABLE :: svobs(:,:,:)
  REAL, ALLOCATABLE :: stobs(:,:,:)

  ! buffer for transpostions
  REAL, POINTER :: sdobs3_global(:,:,:) 
  REAL, POINTER :: svobs3_global(:,:,:)
  REAL, POINTER :: stobs3_global(:,:,:)

  REAL, POINTER :: worksp_global(:,:,:)

  ! fast mode temporary arrays
  logical :: lfill_a, lfill_b
  integer :: ifast_accu

  real, allocatable :: sdfast_a(:,:,:)
  real, allocatable :: svfast_a(:,:,:)
  real, allocatable :: stfast_a(:,:,:)

  real, allocatable :: sdfast_b(:,:,:)
  real, allocatable :: svfast_b(:,:,:)
  real, allocatable :: stfast_b(:,:,:)

  real, allocatable :: sdfast_accu(:,:,:)
  real, allocatable :: svfast_accu(:,:,:)
  real, allocatable :: stfast_accu(:,:,:)

  ! residue correction term
  LOGICAL :: lrescor_a, lrescor_b

  real, allocatable :: sdres_a(:,:,:)
  real, allocatable :: svres_a(:,:,:)
  real, allocatable :: stres_a(:,:,:)

  real, allocatable :: sdres_b(:,:,:)
  real, allocatable :: svres_b(:,:,:)
  real, allocatable :: stres_b(:,:,:)

  ! SITE calculation
  real,allocatable :: sd_o_n0(:,:,:)
  real,allocatable :: sv_o_n0(:,:,:)
  real,allocatable :: st_o_n0(:,:,:)
  real,allocatable :: sd_o_n1(:,:,:)
  real,allocatable :: sv_o_n1(:,:,:)
  real,allocatable :: st_o_n1(:,:,:)

  real,allocatable :: sd_m_n0(:,:,:)
  real,allocatable :: sv_m_n0(:,:,:)
  real,allocatable :: st_m_n0(:,:,:)
  real,allocatable :: sd_m_n1(:,:,:)
  real,allocatable :: sv_m_n1(:,:,:)
  real,allocatable :: st_m_n1(:,:,:)

END MODULE mo_nudging_buffer
