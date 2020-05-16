!---------------------------------------------------------!
!  IMAS, COPYRIGHT (C) 1996 BY GKSS                       !
!  THIS PROGRAM IS EXPERIMENTAL AND IS PROVIDED "AS IS"   !
!  WITHOUT REPRESENTATION OF WARRANTY OF ANY KIND, EITHER !
!  EXPRESS OR IMPLIED. THE ENTIRE RISK AS TO THE QUALITY  !
!  AND PERFORMANCE OF THE PROGRAM IS WITH THE USER.       !
!---------------------------------------------------------!
 
module unit_number_module

!----------------------------------------------------------------------------
! This module keeps a list of unit numbers used in the application.
! In order to avoid collisions the application should ask for
! unit numbers by:
!     get_unit_number or
!     reserve_unit_number (unit) (only if a certain unit number is required)
! and return unit numbers by:
!     return_unit_number (unit)
!
! Cray f90 Compiler:
!   0, 5, 6           : associated to stderr, stdin, stdout, 
!                       but may be reassigned
!   0 - 99, 105 - 299 : available for user specification
!   100 - 104         : reserved for system use
!
! NAG f90 Compiler:
!   0, 5, 6           : preconnected to stderr, stdin, stdout
!   63                : maximum unit number
!
! This module:
!   1 - 4, 7 - 63     : available for user specification
!   5, 6              : reserved
!----------------------------------------------------------------------------
use mo_system ! only : abort
implicit none
private

public :: get_unit_number,     &
          return_unit_number,  &
          reserve_unit_number, &
          print_used_units


  integer, parameter :: highest_unit_number = 63
  logical            :: initialized         = .false.

  type unit_number_list_type
    logical, dimension (highest_unit_number) :: free
  end type unit_number_list_type

  type (unit_number_list_type) :: unit_number_list = &
        unit_number_list_type (.false.)

!==============================================================================
contains
!==============================================================================
  subroutine init_unit_number_list
    !-------------------------------
    ! Initialize unit_number_list:
    ! all units besides 5,6 are free
    !-------------------------------
    unit_number_list% free    = .true.
    unit_number_list% free(5) = .false.
    unit_number_list% free(6) = .false.
    initialized               = .true.
  end subroutine init_unit_number_list
!------------------------------------------------------------------------------
  function reserve_unit_number (unit) result (ok)
  !-------------------------------------------------
  ! Reserves unit number 'unit' for the application
  ! Returns .true. if 'unit' is not in use
  !-------------------------------------------------
  integer :: unit
  logical :: ok
    if (.not. initialized) call init_unit_number_list
    select case (unit)
    case (1:highest_unit_number)
      ok = unit_number_list% free (unit)
      unit_number_list% free (unit) = .false.
    case default
      ok = .false.
    end select
  end function reserve_unit_number
!------------------------------------------------------------------------------
  subroutine return_unit_number (unit)
  !-------------------------------------------------
  ! Returns unit number 'unit' if not used further.
  ! aborts if 'unit' has never been required before.
  !-------------------------------------------------
  integer, intent (inout) :: unit
    if (.not. initialized) call init_unit_number_list
    select case (unit)
    case (1:highest_unit_number)
      if (unit_number_list% free (unit) ) then
        print *,'ERROR in return_unit_number'
        print *,'  unit number',unit,'was never used'
        call abort ('return_unit_number')
      endif
      unit_number_list% free (unit) = .true.
    case default
      print *,'ERROR in return_unit_number'
      print *,'  wrong unit number',unit
      call abort ('return_unit_number')
    end select
    unit = 0
  end subroutine return_unit_number
!------------------------------------------------------------------------------
  function get_unit_number () result (unit)
  !--------------------------------------------
  ! Returns a free unit number
  ! Returns '0' if all unit numbers are in use
  !--------------------------------------------
  integer :: unit
    integer :: i
    if (.not. initialized) call init_unit_number_list
    unit = 0
    do i = highest_unit_number, 1, -1
      if (unit_number_list% free (i)) then
        unit = i
        unit_number_list% free (i) = .false.
        exit
      endif
    end do
  end function get_unit_number
!------------------------------------------------------------------------------
  subroutine print_used_units
    integer :: i
    if (.not. initialized) call init_unit_number_list
    print *,'used unit numbers:'
    do i = 1, highest_unit_number
      if (.not. unit_number_list% free(i)) write(*,'(i4)',advance='NO') i
    end do
    print *
  end subroutine print_used_units
!------------------------------------------------------------------------------
end module unit_number_module

