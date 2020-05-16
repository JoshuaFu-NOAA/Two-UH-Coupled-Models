MODULE mo_f5
  ! pointers for *f5* space.
  REAL, TARGET,  ALLOCATABLE :: fszl(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fazl(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fszm(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fazm(:,:,:)
END MODULE mo_f5
MODULE mo_f7
  ! pointers for *f7* space.
  REAL, TARGET,  ALLOCATABLE :: fsdl(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fadl(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fsdm(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fadm(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fsr(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: far(:,:,:)
END MODULE mo_f7
MODULE mo_f8
  ! pointers for *f8* space.
  REAL, TARGET,  ALLOCATABLE :: fstp1(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fatp1(:,:,:)
  REAL, TARGET,  ALLOCATABLE :: fsul(:)
  REAL, TARGET,  ALLOCATABLE :: faul(:)
END MODULE mo_f8

MODULE mo_f
  USE mo_f5
  USE mo_f7
  USE mo_f8
END MODULE mo_f
