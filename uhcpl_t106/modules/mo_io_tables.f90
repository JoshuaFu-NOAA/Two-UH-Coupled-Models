MODULE mo_io_tables

  USE mo_parameters

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  !
  ! *mo_io_tables* - module containing input/output tables.
  !
  ! j. k. gibson     e.c.m.w.f.     17/11/81.
  !
  ! ------------------------------------------------------------------

  INTEGER :: ng3xp            !   no. of extra *g3 fields
  INTEGER :: ng3xl(jpg3xf)    !   no. of levels of each extra *g3 field
  INTEGER :: ng4xp            !   no. of extra *g4 fields
  INTEGER :: ng4xl(jpg3xf)    !   no. of levels of each extra *g4 field

END MODULE mo_io_tables
