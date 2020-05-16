MODULE mo_grid

  IMPLICIT NONE

  !
  !  basic grid point resolution parameters
  !

!  INTEGER, PARAMETER :: nxpt=1   ! no.of points outside active domain for 
!                                 ! interpolant
!  INTEGER, PARAMETER :: jintmx=1 ! number of extra latitudes in polar region

! choose the max. value for nxpt depending on resolution
!
! recommended: t21      : 18
!              t106     : 60
!              nprocb=1 : 2       ! 1 PE in E-W direction
!
  INTEGER, PARAMETER :: nxpt=60
! INTEGER, PARAMETER :: nxpt=18
! INTEGER, PARAMETER :: nxpt=2

  INTEGER, PARAMETER :: joverlap=3
  INTEGER, PARAMETER :: jintmx=joverlap-nxpt

  !  integer pointers to physical data starting locations in 3-d arrays
  INTEGER, PARAMETER :: i1  = 1 + nxpt           ! model starting longitude
  INTEGER, PARAMETER :: j1  = 1 + nxpt + jintmx  ! model starting latitude

  INTEGER ::   plon    ! number of local longitudes
  INTEGER ::   plev    ! number of vertical levels
  INTEGER ::   plat    ! number of local latitudes
  INTEGER ::   plato2  ! number of local latitudes/2
  INTEGER ::   pcnst   ! number of constituents (including water vapor)
  INTEGER ::   plevp1  ! plev + 1
  INTEGER ::   plond   ! local slt extended domain longitude
  INTEGER ::   platd   ! local slt extended domain latitude per hemisphere
  INTEGER ::   plevd   ! fold plev,pcnst indices into one
  INTEGER ::   pgls    ! length of local latitude slice
  INTEGER ::   jfirst  ! first index to be computed
  INTEGER ::   jlast   ! last  index to be computed
  INTEGER ::   i2pi    ! start of eastern long. extension
  INTEGER ::   plono2  ! plon/2
  INTEGER ::   istart  ! index to start computation
  INTEGER ::   istop   ! index to stop  computation
  INTEGER ::   js      ! index of southernmost model lat
  INTEGER ::   jn      ! index of northernmost model lat
  INTEGER ::   jstart  ! index for first model lat.
  INTEGER ::   jstop   ! index for last  model lat.
  INTEGER ::   pbpts   ! (length of local latitude slice)*fields
  INTEGER ::   plevm1
  REAL    ::   dphibr  ! reciprocal of maximum del phi
  
END MODULE mo_grid

