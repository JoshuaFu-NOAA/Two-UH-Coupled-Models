MODULE mo_filename

  IMPLICIT NONE

  ! -----------------------------------------------------------------
  !
  ! module *mo_filename* - quantities needed for file names etc.
  !
  ! -----------------------------------------------------------------
  !
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  !

  INTEGER, PARAMETER :: path_limit=256, name_limit=19
  INTEGER, PARAMETER :: dir_limit = path_limit-name_limit

  INTEGER :: nhm
  INTEGER :: nhy
  INTEGER :: nhd

  CHARACTER (dir_limit) :: yomdn
  CHARACTER (dir_limit) :: ypath
  CHARACTER (9)         :: yexp

  CHARACTER (path_limit) :: standard_grib_file
  CHARACTER (path_limit) :: extended_grib_file

END MODULE mo_filename
