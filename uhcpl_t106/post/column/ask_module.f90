!---------------------------------------------------------!
!  IMAS, COPYRIGHT (C) 1996 BY GKSS                       !
!  THIS PROGRAM IS EXPERIMENTAL AND IS PROVIDED "AS IS"   !
!  WITHOUT REPRESENTATION OF WARRANTY OF ANY KIND, EITHER !
!  EXPRESS OR IMPLIED. THE ENTIRE RISK AS TO THE QUALITY  !
!  AND PERFORMANCE OF THE PROGRAM IS WITH THE USER.       !
!---------------------------------------------------------!

module ask_module
use mo_system         ! only : abort, flush
use unit_number_module, only : get_unit_number, return_unit_number

  implicit none

  private

  public :: ring_bell,           &! ring the bell 
            ask,                 &! ask for anything
            ask_y_n,             &! ask for a logical (answer 'y' or 'n')
            ask_character,       &! ask for a character
            ask_integer,         &! ask for an integer
            ask_real,            &! ask for a real
            ask_random_seed,     &! ask for seed of system clock
            ask_menue             ! enter menue

  public :: log_on,              &! start writing to logfile
            read_log,            &! start reading from logfile
            log_off,             &! stop reading or writing logfile
            get_logunit,         &! get unit_number of logfile
            get_log               ! get logfile status

  public :: menu_item_length,    &! length of menu item character string
            max_menu_items,      &! max. number of menu items
            stdout_unit           ! unit for standard output (default=6)

!------------------------------------------------------------------------------
  !----------------
  ! Kind parameters
  !----------------
  integer, parameter :: dp = selected_real_kind(13) ! Double Precision
  integer, parameter :: sp = selected_real_kind(6)  ! Single Precision
  integer, parameter :: wp = dp
  !--------------------------------------
  ! Define non-printable ascii characters
  !--------------------------------------
  !                       mnemonic  decimal ! octal
  character, parameter :: HT  = achar ( 09) ! 011
  character, parameter :: BEL = achar ( 07) ! 007  
!------------------------------------------------------------------------------
! host variables

  integer,          private            :: iostat
  character(len=*), private, parameter :: deli =' '//HT
  integer,          private            :: logunit           = 0
  logical,          private            :: readlog           = .false.
  logical,          private            :: writelog          = .false.

  integer                       :: stdout_unit       = 6
  integer,            parameter :: menu_item_length  = 8
  integer,            parameter :: max_menu_items    = 30
  character (len=*),  parameter :: menu_item_fmt     ='(1x,a8,": ",a)'
!------------------------------------------------------------------------------
! Interface

  interface ask
    module procedure subr_ask_string,         &! ask for a string
                     subr_ask_integer,        &! ask for an integer
                     subr_ask_real_dp,        &! ask for a real(dp) number
                     subr_ask_real_sp,        &! ask for a real(sp) number
                     subr_ask_logical_vector, &! ask for a vector of logicals
                     subr_ask_integer_vector, &! ask for a vector of integer
                     subr_ask_real_vector_dp, &! ask for a vector of reals
                     subr_ask_real_vector_sp   ! ask for a vector of reals
  end interface

  interface ask_y_n
    module procedure fct_ask_y_n
  end interface

  interface ask_character
    module procedure fct_ask_character
  end interface

  interface ask_integer
    module procedure fct_ask_integer
  end interface

  interface ask_real
    module procedure fct_ask_real
  end interface



!---------
 contains
!---------

!==============================================================================
!                      INTERACTIVE INPUT (SUBROUTINES)
!==============================================================================
  subroutine pline(unit ,c, text)
  integer, intent(in) :: unit
  character,         optional, intent(in) :: c
  character (len=*) ,optional, intent(in) :: text
    character (len=72) :: l
    l = line (c, text)
    write (unit,'(a72)') l
  end subroutine pline
!------------------------------------------------------------------------------
  function line (c, text)
  character,         optional, intent(in) :: c
  character (len=*) ,optional, intent(in) :: text
  character (len=72)                      :: line
  !----------------------------------------------------
  ! returns a line with unique character '-'
  ! if 'c' is present the character is 'c'
  ! if text is present the text is centered on the line
  !----------------------------------------------------
    character (len=72) :: string
    integer            :: l
    if (present(c)) then
      line = repeat( c ,72)
    else
      line = repeat('-',72)
    endif
    if(present(text)) then
      string = adjustl(text)
      l = len_trim (string)
      line(36-l/2  :35-l/2+l) = string
      line(35-l/2  :35-l/2)   = ' '
      line(36-l/2+l:36-l/2+l) = ' '
    endif
  end function line
!------------------------------------------------------------------------------
  subroutine ring_bell (times)
  integer, optional :: times
    integer              :: i
    if (present(times)) then
      do i = 1, times
        write (stdout_unit,'(a1)',advance='NO') BEL
      end do
    else
      write (stdout_unit,'(a1)',advance='NO') BEL
    endif
  end subroutine ring_bell
!------------------------------------------------------------------------------
  subroutine subr_ask_string (a, q, opt1, opt2, opt3, opt4, opt5 ,  &
                                    opt6, opt7, opt8, opt9, opt10)
  character (len=*), intent (out)          :: a
  character (len=*), intent (in)           :: q
  character (len=*), intent (in), optional :: opt1, opt2, opt3, opt4, opt5,  &
                                              opt6, opt7, opt8, opt9, opt10
    integer :: i
     do
       a = ' '
       if (readlog) then
         write (stdout_unit,*) q
         write (stdout_unit,'(">")',advance='NO')
         read  (logunit,'(a)',iostat=iostat) a
         if (iostat/=0)then
           call read_error
           cycle
         else
           write (stdout_unit,      '(a)') a
         endif
       else
         write (stdout_unit,*) q
         if (present(opt1 )) write (stdout_unit,*) opt1
         if (present(opt2 )) write (stdout_unit,*) opt2
         if (present(opt3 )) write (stdout_unit,*) opt3
         if (present(opt4 )) write (stdout_unit,*) opt4
         if (present(opt5 )) write (stdout_unit,*) opt5
         if (present(opt6 )) write (stdout_unit,*) opt6
         if (present(opt7 )) write (stdout_unit,*) opt7
         if (present(opt8 )) write (stdout_unit,*) opt8
         if (present(opt9 )) write (stdout_unit,*) opt9
         if (present(opt10)) write (stdout_unit,*) opt10
         write (stdout_unit,'(">")',advance='NO')
         read (5,'(a)',iostat=iostat) a
         if (iostat /= 0) then
            write (stdout_unit,*) BEL,'invalid input, please type again!'
            cycle
         end if
       endif
       exit
     end do
     if (writelog) then
       write (logunit,'(a,a,a)') a, HT//'! ', q
       call flush (logunit)
     endif
     i = scan(a,' '//HT)
     if (i/=0) a(i:) = ' '
  end subroutine subr_ask_string
!------------------------------------------------------------------------------
  subroutine subr_ask_logical_vector (a, q)
  logical,   intent(out) :: a(:)
  character, intent(in)  :: q
    do
      write (stdout_unit,*) q
      write (stdout_unit,'(">")',advance='NO')
      if (readlog) then
        read  (logunit,*,iostat=iostat) a
        if (iostat/=0) then
          call read_error
          cycle
        else
          write (stdout_unit,      *) a
        endif
      else
        read (5,*,iostat=iostat) a
        if (iostat /= 0) then
          write (stdout_unit,*) BEL,'invalid input, please type again!'
          cycle
        end if
      endif
      exit
    end do
    if (writelog) then
      write (logunit,*) a, HT//'!', q
      call flush (logunit)
    endif
  end subroutine subr_ask_logical_vector
!------------------------------------------------------------------------------
  subroutine subr_ask_integer (a, q, opt1, opt2, opt3, opt4, opt5,  &
                                     opt6, opt7, opt8, opt9, opt10)
  integer,           intent(out)          :: a
  character (len=*), intent(in)           :: q
  character (len=*), intent(in), optional :: opt1, opt2, opt3, opt4, opt5,  &
                                             opt6, opt7, opt8, opt9, opt10
      do
        write (stdout_unit,*) q
        if (present(opt1 )) write (stdout_unit,*) opt1
        if (present(opt2 )) write (stdout_unit,*) opt2
        if (present(opt3 )) write (stdout_unit,*) opt3
        if (present(opt4 )) write (stdout_unit,*) opt4
        if (present(opt5 )) write (stdout_unit,*) opt5
        if (present(opt6 )) write (stdout_unit,*) opt6
        if (present(opt7 )) write (stdout_unit,*) opt7
        if (present(opt8 )) write (stdout_unit,*) opt8
        if (present(opt9 )) write (stdout_unit,*) opt9
        if (present(opt10)) write (stdout_unit,*) opt10
        write (stdout_unit,'(">")',advance='NO')
        if (readlog) then
          read  (logunit,*,iostat=iostat) a
          if (iostat/=0) then
            call read_error
            cycle
          else
            write (stdout_unit,      *) a
          endif
        else
          read (5,*,iostat=iostat) a
          if (iostat /= 0) then
             write (stdout_unit,*) BEL,'invalid input, please type again!'
             cycle
          end if
        endif
        exit
      end do
      if (writelog) then
        write (logunit,*) a, HT//'!', q
        call flush (logunit)
      endif
  end subroutine subr_ask_integer
!------------------------------------------------------------------------------
  subroutine subr_ask_integer_vector (a, q, write_log)
  integer,           intent(out)          :: a(:)
  character (len=*), intent(in)           :: q
  logical,           intent(in), optional :: write_log
      do
        write (stdout_unit,*) q
        write (stdout_unit,'(">")',advance='NO')
        if (readlog) then
          read  (logunit,*,iostat=iostat) a
          if (iostat/=0) then
            call read_error
            cycle
          else
            write (stdout_unit,      *) a
          endif
        else
          read (5,*,iostat=iostat) a
          if (iostat /= 0) then
             write (stdout_unit,*) BEL,'invalid input, please type again!'
             cycle
          end if
        endif
        exit
      end do
      if (writelog) then
        if (present(write_log)) then
          if (write_log) then
            write (logunit,*) a, HT//'!', q
            call flush (logunit)
          endif
        else
          write (logunit,*) a, HT//'!', q
          call flush (logunit)
        endif
      endif
  end subroutine subr_ask_integer_vector
!------------------------------------------------------------------------------
  subroutine subr_ask_real_vector_dp (a, q)
  real(DP),          intent(out) :: a(:)
  character (len=*), intent(in)  :: q
      do
        write (stdout_unit,*) q
        write (stdout_unit,'(">")',advance='NO')
        if (readlog) then
          read  (logunit,*,iostat=iostat) a
          if (iostat/=0) then
            call read_error
            cycle
          else  
            write (stdout_unit,      *) a
          endif
        else
          read (5,*,iostat=iostat) a
          if (iostat /= 0) then
             write (stdout_unit,*) BEL,'invalid input, please type again!'
             cycle
          end if
        endif
        exit
      end do
      if (writelog) then
        write (logunit,*) a, HT//'!', q
        call flush (logunit)
      endif
  end subroutine subr_ask_real_vector_dp
!------------------------------------------------------------------------------
  subroutine subr_ask_real_dp (a, q, opt1, opt2, opt3, opt4, opt5,  &
                                  opt6, opt7, opt8, opt9, opt10)
  real(DP),          intent(out) :: a   
  character (len=*), intent(in)  :: q
  character (len=*), intent(in), optional :: opt1, opt2, opt3, opt4, opt5,  &
                                             opt6, opt7, opt8, opt9, opt10
      do!e
        write (stdout_unit,*) q
        if (present(opt1 )) write (stdout_unit,*) opt1
        if (present(opt2 )) write (stdout_unit,*) opt2
        if (present(opt3 )) write (stdout_unit,*) opt3
        if (present(opt4 )) write (stdout_unit,*) opt4
        if (present(opt5 )) write (stdout_unit,*) opt5
        if (present(opt6 )) write (stdout_unit,*) opt6
        if (present(opt7 )) write (stdout_unit,*) opt7
        if (present(opt8 )) write (stdout_unit,*) opt8
        if (present(opt9 )) write (stdout_unit,*) opt9
        if (present(opt10)) write (stdout_unit,*) opt10
        write (stdout_unit,'(">")',advance='NO')
        if (readlog) then
          read  (logunit,*,iostat=iostat) a
          if (iostat/=0) then
            call read_error
            cycle
          else
            write (stdout_unit,      *) a
          endif
        else
          read (5,*,iostat=iostat) a
          if (iostat /= 0) then
             write (stdout_unit,*) BEL,'invalid input, please type again!'
             cycle
          end if
        endif
        exit
      end do
      if (writelog) then
        write (logunit,*) a, HT//'!', q
        call flush (logunit)
      endif
  end subroutine subr_ask_real_dp
!------------------------------------------------------------------------------
  subroutine subr_ask_real_vector_sp (a, q)
  real(SP),          intent(out) :: a(:)
  character (len=*), intent(in)  :: q
      do
        write (stdout_unit,*) q
        write (stdout_unit,'(">")',advance='NO')
        if (readlog) then
          read  (logunit,*,iostat=iostat) a
          if (iostat/=0) then
            call read_error
            cycle
          else  
            write (stdout_unit,      *) a
          endif
        else
          read (5,*,iostat=iostat) a
          if (iostat /= 0) then
             write (stdout_unit,*) BEL,'invalid input, please type again!'
             cycle
          end if
        endif
        exit
      end do
      if (writelog) then
        write (logunit,*) a, HT//'!', q
        call flush (logunit)
      endif
  end subroutine subr_ask_real_vector_sp
!------------------------------------------------------------------------------
  subroutine subr_ask_real_sp (a, q, opt1, opt2, opt3, opt4, opt5,  &
                                  opt6, opt7, opt8, opt9, opt10)
  real(SP),          intent(out) :: a   
  character (len=*), intent(in)  :: q
  character (len=*), intent(in), optional :: opt1, opt2, opt3, opt4, opt5,  &
                                             opt6, opt7, opt8, opt9, opt10
      do!e
        write (stdout_unit,*) q
        if (present(opt1 )) write (stdout_unit,*) opt1
        if (present(opt2 )) write (stdout_unit,*) opt2
        if (present(opt3 )) write (stdout_unit,*) opt3
        if (present(opt4 )) write (stdout_unit,*) opt4
        if (present(opt5 )) write (stdout_unit,*) opt5
        if (present(opt6 )) write (stdout_unit,*) opt6
        if (present(opt7 )) write (stdout_unit,*) opt7
        if (present(opt8 )) write (stdout_unit,*) opt8
        if (present(opt9 )) write (stdout_unit,*) opt9
        if (present(opt10)) write (stdout_unit,*) opt10
        write (stdout_unit,'(">")',advance='NO')
        if (readlog) then
          read  (logunit,*,iostat=iostat) a
          if (iostat/=0) then
            call read_error
            cycle
          else
            write (stdout_unit,      *) a
          endif
        else
          read (5,*,iostat=iostat) a
          if (iostat /= 0) then
             write (stdout_unit,*) BEL,'invalid input, please type again!'
             cycle
          end if
        endif
        exit
      end do
      if (writelog) then
        write (logunit,*) a, HT//'!', q
        call flush (logunit)
      endif
  end subroutine subr_ask_real_sp
!------------------------------------------------------------------------------
  subroutine ask_random_seed
  !----------------------------------------------------------------
  ! sets random generator seed to the same value in subsequent runs
  !----------------------------------------------------------------
    integer                            :: siz
    integer, dimension(:), allocatable :: seed
    logical                            :: answer
    call random_seed (size=siz)
    allocate (seed(siz))
    call random_seed (get=seed)
    write (stdout_unit,*) 'actual random_seed = ',seed
    !------------------------------------
    ! change seed if reading from stdin ?
    !------------------------------------
    if (.not. readlog) then
      if (ask_y_n('change seed?',write_log=.false.)) then
        call ask(seed, 'enter new seed',write_log=.false.)
        call random_seed (put=seed)
        write (stdout_unit,*) 'new random_seed = ',seed
      endif
    endif
    !----------------------
    ! write seed to logfile
    !----------------------
    if (writelog) then
      write (logunit,*) seed, HT//'! random_seed'
      call flush (logunit)
      answer = ask_y_n &
        ('set random generator seed to random_seed in subsequent runs?')
    elseif (readlog) then
    !-----------------------
    ! read seed from logfile
    !-----------------------
      call ask (seed,' ')
      if (ask_y_n(' ')) then
        call random_seed (put=seed)
        write (stdout_unit,*) 'restored random_seed = ',seed
      endif
    endif
    deallocate (seed)    
  end subroutine ask_random_seed
!==============================================================================
!                      INTERACTIVE INPUT (FUNCTIONS)
!==============================================================================
  function fct_ask_y_n (q,write_log) result (yes)
  character (len=*) :: q
  logical, optional :: write_log
  logical           :: yes
    character :: c
    integer   :: iostat
    do
      write (stdout_unit,*) q
      write (stdout_unit,'(">")',advance='NO')
      if (readlog) then
        read  (logunit,*,iostat=iostat) c
        if (iostat/=0) then
          call read_error
          cycle
        else
          write (stdout_unit,      *) c
        endif
      else
        read (*,*,iostat=iostat) c
        if (iostat /= 0) c = 'x'
      endif
      select case (c)
        case ('y','Y')
          yes = .true.
          exit
        case ('n','N')
          yes = .false.
          exit
        case default
          write (stdout_unit,*) BEL,'invalid input, please type again! (y,Y,n,N)'
          if (readlog) call read_error
          cycle
      end select
    end do
      if (writelog) then
        if (present(write_log)) then
          if (write_log) then
            write (logunit,*) c, HT//'!', q
            call flush (logunit)
          endif
        else
          write (logunit,*) c, HT//'!', q
          call flush (logunit)
        endif
      endif
  end function fct_ask_y_n
!------------------------------------------------------------------------------
  function fct_ask_character (q, opt1, opt2, opt3, opt4, opt5 ,  &
                                 opt6, opt7, opt8, opt9, opt10)    result (a)
  character (len=*)           :: q
  character (len=*), optional :: opt1, opt2, opt3, opt4, opt5,  &
                                 opt6, opt7, opt8, opt9, opt10
  character (len=1)           :: a

        call ask (a, q, opt1, opt2, opt3, opt4, opt5 ,    &
                        opt6, opt7, opt8, opt9, opt10) 

  end function fct_ask_character
!------------------------------------------------------------------------------
  function fct_ask_integer (q, opt1, opt2, opt3, opt4, opt5,  &
                               opt6, opt7, opt8, opt9, opt10)   result (a)
  character (len=*)           :: q
  character (len=*), optional :: opt1, opt2, opt3, opt4, opt5,  &
                                 opt6, opt7, opt8, opt9, opt10
  integer                     :: a

         call ask (a, q, opt1, opt2, opt3, opt4, opt5,    &
                         opt6, opt7, opt8, opt9, opt10)

  end function fct_ask_integer
!------------------------------------------------------------------------------
  function fct_ask_real (q, opt1, opt2, opt3, opt4, opt5,  &
                            opt6, opt7, opt8, opt9, opt10)   result (a)
  character (len=*)           :: q
  character (len=*), optional :: opt1, opt2, opt3, opt4, opt5,  &
                                 opt6, opt7, opt8, opt9, opt10
  real(WP)                    :: a   

            call ask (a, q, opt1, opt2, opt3, opt4, opt5,  &
                            opt6, opt7, opt8, opt9, opt10)

  end function fct_ask_real
!------------------------------------------------------------------------------
  function ask_menue (q,c01,c02,c03,c04,c05,c06,c07,c08,c09,c10,&
                        c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,&
                        c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,&
                        title)
  !-------------------------------------------------------------
  ! Displays a menu.
  ! Returns a character string containing the first word of the
  ! menue item choosen 
  !-------------------------------------------------------------
  character(len=menu_item_length)      :: ask_menue
  character(len=*),intent(in)          :: q
  character(len=*),intent(in)          :: c01
  character(len=*),intent(in),optional :: c02  
  character(len=*),intent(in),optional :: c03
  character(len=*),intent(in),optional :: c04
  character(len=*),intent(in),optional :: c05
  character(len=*),intent(in),optional :: c06
  character(len=*),intent(in),optional :: c07
  character(len=*),intent(in),optional :: c08
  character(len=*),intent(in),optional :: c09
  character(len=*),intent(in),optional :: c10
  character(len=*),intent(in),optional :: c11
  character(len=*),intent(in),optional :: c12
  character(len=*),intent(in),optional :: c13
  character(len=*),intent(in),optional :: c14
  character(len=*),intent(in),optional :: c15
  character(len=*),intent(in),optional :: c16
  character(len=*),intent(in),optional :: c17
  character(len=*),intent(in),optional :: c18
  character(len=*),intent(in),optional :: c19
  character(len=*),intent(in),optional :: c20
  character(len=*),intent(in),optional :: c21
  character(len=*),intent(in),optional :: c22
  character(len=*),intent(in),optional :: c23
  character(len=*),intent(in),optional :: c24
  character(len=*),intent(in),optional :: c25
  character(len=*),intent(in),optional :: c26
  character(len=*),intent(in),optional :: c27
  character(len=*),intent(in),optional :: c28
  character(len=*),intent(in),optional :: c29
  character(len=*),intent(in),optional :: c30
  character(len=*),intent(in),optional :: title
    !----------------
    ! local variables
    !----------------
    character(len=menu_item_length)      :: answer
    integer           :: found
    integer           :: iostat
    integer           :: i
    do
      call print_menue
!..........................................................................
!      call prompt !    INLINED internal subroutine (sgi bug workaround)
!..........................................................................
!    subroutine prompt
      write (stdout_unit,'(a,"> ")',advance='NO') trim(q)
      if (readlog) then
        read  (logunit,'(a)',iostat=iostat) answer
        if (iostat/=0)then
          call read_error
        else
          write (stdout_unit,      '(a)') answer
        endif
      else
        read  (5,'(a)',iostat=iostat) answer
      endif
      i = scan(answer,' '//HT)
      if (i/=0) answer(i:) = ' '
!    end subroutine prompt
!..........................................................................
      if(iostat==0) then
        call search
        if (found == 1) exit
      endif
      write(stdout_unit,*) 'Invalid input, type again!'
    end do
    if (writelog) then
      write (logunit,'(a,a,a)') answer, HT//'! ', q
      call flush (logunit)
    endif
    ask_menue = answer
  contains
    !..........................................................................
    subroutine print_menue
      call pline(stdout_unit,'=',title)
      call print_item (c01)
      call print_item (c02)
      call print_item (c03)
      call print_item (c04)
      call print_item (c05)
      call print_item (c06)
      call print_item (c07)
      call print_item (c08)
      call print_item (c09)
      call print_item (c10)
      call print_item (c11)
      call print_item (c12)
      call print_item (c13)
      call print_item (c14)
      call print_item (c15)
      call print_item (c16)
      call print_item (c17)
      call print_item (c18)
      call print_item (c19)
      call print_item (c20)
      call print_item (c21)
      call print_item (c22)
      call print_item (c23)
      call print_item (c24)
      call print_item (c25)
      call print_item (c26)
      call print_item (c27)
      call print_item (c28)
      call print_item (c29)
      call print_item (c30)
    end subroutine print_menue
    !..........................................................................
    subroutine print_item (c)
    character(len=*), intent(in),optional :: c
      integer :: i1,i2
      if (present(c)) then
        i1 = index (c       ,' ') - 1
        i2 = verify(c(i1+2:),' ') + i1
        if (i1>0) then
          write (stdout_unit,menu_item_fmt) c(:i1),c(i2:)
        else if (c/=' ') then
          call pline(stdout_unit,'-',c)
        endif
      endif
    end subroutine print_item
    !..........................................................................
    subroutine search
      found = 0
      call search_item (c01)
      call search_item (c02)
      call search_item (c03)
      call search_item (c04)
      call search_item (c05)
      call search_item (c06)
      call search_item (c07)
      call search_item (c08)
      call search_item (c09)
      call search_item (c10)
      call search_item (c11)
      call search_item (c12)
      call search_item (c13)
      call search_item (c14)
      call search_item (c15)
      call search_item (c16)
      call search_item (c17)
      call search_item (c18)
      call search_item (c19)
      call search_item (c20)
      call search_item (c21)
      call search_item (c22)
      call search_item (c23)
      call search_item (c24)
      call search_item (c25)
      call search_item (c26)
      call search_item (c27)
      call search_item (c28)
      call search_item (c29)
      call search_item (c30)
      if (found>1) &
        call abort('ask_menue: more than one entries == '//answer)
    end subroutine search
    !..........................................................................
    subroutine search_item (c)
    character(len=*), intent(in),optional :: c
      integer :: i1
      if (present(c)) then
        i1 = index (c       ,' ')-1
        if (i1>0 .and. answer==c(1:i1)) found = found + 1
      endif
    end subroutine search_item
    !..........................................................................
  end function ask_menue
!==============================================================================
! Log - File
!==============================================================================
  subroutine log_on (position)
  character (len=*), optional :: position
    logunit = get_unit_number()
    if (present(position)) then
      open (logunit, file='logfile', action='WRITE', position=position)
    else
      open (logunit, file='logfile', action='WRITE')
    endif
    writelog = .true.
  end subroutine log_on
!------------------------------------------------------------------------------
  subroutine read_error
    character :: c
    integer   :: iostat
    write (stdout_unit,*)' '
    write (stdout_unit,*)'ERROR while reading logfile !'
    write (stdout_unit,*)'continue reading from stdin'
    write (stdout_unit,"('>')",advance = 'NO')
    read (*,*,iostat=iostat) c
    if(iostat /= 0) call abort ('read_error')
    select case (c)
    case ('y','Y')
      call log_off
      if (ask_y_n('continue writing to logunit')) then
        call log_on ('APPEND')
        write (logunit,'("!------------------------------------------------")')
        write (logunit,'("! appended dialog:")')
      endif
    case default
      call abort ('read_error')
    end select
  end subroutine read_error
!------------------------------------------------------------------------------
  subroutine log_off
    close (logunit)
    call return_unit_number(logunit)
    logunit  = 0
    readlog  = .false.
    writelog = .false.
  end subroutine log_off
!------------------------------------------------------------------------------
  subroutine read_log
    logunit = get_unit_number()
    open (logunit, file='logfile', action='READ')
    readlog = .true.
  end subroutine read_log
!------------------------------------------------------------------------------
  function get_logunit ()
  integer :: get_logunit
    get_logunit = logunit
  end function get_logunit
!------------------------------------------------------------------------------
  function get_log ()
  character :: get_log
    get_log               = ' '
    if (writelog) get_log = 'w' 
    if (readlog)  get_log = 'r'
  end function get_log
!------------------------------------------------------------------------------
  subroutine error_in_logfile
    write (stdout_unit,*)'ERROR while reading logfile, now reading stdin.'
    readlog  = .false.
    close (logunit)
    call return_unit_number (logunit)
    logunit  = 0 
  end subroutine error_in_logfile
!------------------------------------------------------------------------------
end module ask_module
