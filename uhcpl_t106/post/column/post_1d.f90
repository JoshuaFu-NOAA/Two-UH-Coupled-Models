program post_1d
!
!+ $Id: post_1d.f90,v 1.13 2000/05/18 13:18:06 m214030 Exp $
!
! Authors:
!
! A. Rhodin, MPI, November 1999, original source
!
  !-------------
  ! modules used
  !-------------
  use mo_system !only: abort, flush
  use mo_1d_io,  only: readfile, & ! subroutine to read column model files
                       sp, wp      ! single precision, working precision
  use ask_module
  implicit none
  !--------------------------------
  ! data type for list of variables
  !--------------------------------
  type t_list
    character(len=16) :: name  ! name of variable
    integer           :: size  ! size of array variable
    integer           :: n     ! number of time steps available
    integer           :: temp
  end type t_list
  !-------------------
  ! variables declared
  !-------------------
  character(len=16)   :: file, name
  character(len=128)  :: epsfile
  character           :: c
  integer             :: ios, iunit
  real(wp) ,pointer   :: z(:), x(:,:)
  type(t_list)        :: l(100)
  integer             :: i,j
  integer             :: nvar
  integer             :: gpu
  character(len=16)   :: range
  integer             :: lonlat(2)
  integer             :: index(4)
  logical             :: p_asc, p_bak, opened
  character(len=*) ,parameter :: ccol(2)=(/'color solid','mono dashed'/)
  integer             :: nlev
  !----------------
  ! plot parameters
  !----------------
  logical             :: lmin  = .false. ! subtract mean
  integer             :: t0    = 1       ! first time step plotted
  integer             :: t9    = 9999    ! last  time step plotted
  integer             :: t     = 1       ! time step
  integer             :: k     = 1       ! level 
  integer             :: n     = 0       ! variable index
  character           :: mode  = 'l'     ! plot mode l,t=constant level,time 
  real(sp)            :: mean  = 0._sp
  integer             :: icol  = 1
  character(len=7)    :: m(5)  = (/'forcing','residui',          &
                                   'result1','result2','1d_tend'/)
  !-------------------------------
  ! set some initial values values
  !-------------------------------
  nullify(z)
  stdout_unit =  0 ! output unit for dialogue
  gpu         =  6 ! output unit for gnuplot
  index       =  1
  index(1)    = -1
  lonlat      = (/0,0/) 

  !---------------------
  ! set up variable list
  !---------------------
  l% n = 0; l% size = 0; l% name = ''; nvar=0
  !--------------------------
  ! test if files are present
  !--------------------------

  file='forcing'
  open(1,file=file,form='unformatted',status='old',iostat=ios)
  if(ios==0) then
    m(1) = file
    call fill_list (iunit=1)
  else
    write(stdout_unit,*)'cannot open ',file
    file='forcing.bak'
    inquire(file=file,exist=p_bak)
    if(p_bak) p_bak = ask_y_n ('rewrite forcing file from '//trim(file)//'?')
    if(p_bak) then
      open(1,file=file,form='unformatted',status='old',iostat=ios)
      if(ios==0) then
        !---------------------------------------------------------------------
        ! file 'forcing' not present: derive file 'forcing' from 'forcing.bak'
        !---------------------------------------------------------------------
        write(stdout_unit,*)'opened ',file
        m(1) = file
        rewind(1)
        open(2,file='forcing',form='unformatted',status='new',iostat=ios)
        if(ios/=0) then
          write(stdout_unit,*)"cannot open new file 'forcing'"
          stop
        endif
        do
          call readfile(z ,name ,1 ,index ,lonlat ,ounit=2)
          if(.not.associated(z)) then
            write(0,*)'EOF unit ',1
            exit
          endif
        end do
        close(1)
        close(2)
!        stop
      else
        write(stdout_unit,*)'cannot open ',file
        stop
      endif
    else
      file='forcing.asc'
      inquire(file=file,exist=p_asc)
      if(p_asc) p_asc = ask_y_n ('rewrite residui file from '//trim(file)//'?')
      if(p_asc) then
        open(1,file=file,form='formatted',status='old',iostat=ios)
        if(ios==0) then
          !-----------------------------------------------------------
          ! file 'forcing' not present: derive file from 'forcing.asc'
          !-----------------------------------------------------------
          write(stdout_unit,*)'opened ',file
          m(1) = file
          rewind(1)
          open(2,file='forcing',form='unformatted',status='new',iostat=ios)
          if(ios/=0) then
            write(stdout_unit,*)"cannot open new file 'forcing'"
            stop
          endif
          do
            call readfile(z ,name ,1 ,index ,lonlat ,ounit=2 ,ifmt=.true.)
            if(.not.associated(z)) then
              write(0,*)'EOF unit ',1
              exit
            endif
          end do
          close(1)
          close(2)
!          stop
        else
          write(stdout_unit,*)'cannot open ',file
          stop
        endif
      else
        write(stdout_unit,*)'cannot rewrite "forcing".'
        stop
      endif
    endif
    file='forcing'
    open(1,file=file,form='unformatted',status='old',iostat=ios)
    m(1) = file
    call fill_list (iunit=1)
  endif

  !---------------------------------------------
  ! file 'forcing' present: open remaining files
  !---------------------------------------------


  file='residui'
  open(2,file=file,form='unformatted',status='old',iostat=ios)
  if(ios==0) then
    m(2) = file
    call fill_list (iunit=2)
  else
    file='residui.asc'
    inquire(file=file,exist=p_asc)
    if(p_asc) p_asc = ask_y_n ('rewrite residui file from '//trim(file)//'?')
    if(p_asc) then
      open(2,file=file,form='formatted',status='old',iostat=ios)
      if(ios==0) then
        !-----------------------------------------------------------
        ! file 'residui' not present: derive file from 'residui.asc'
        !-----------------------------------------------------------
        write(stdout_unit,*)'opened ',file
        m(2) = file
        rewind(2)
        open(3,file='residui',form='unformatted',status='new',iostat=ios)
        if(ios/=0) then
          write(stdout_unit,*)"cannot open new file 'residui'"
          stop
        endif
        do
          call readfile(z ,name ,2 ,index ,lonlat ,ounit=3 ,ifmt=.true.)
          if(.not.associated(z)) then
            write(0,*)'EOF unit ',2
            exit
          endif
        end do
        close(2)
        close(3)
      else
        write(stdout_unit,*)'cannot open ',file
        stop
      endif
    endif
  endif

  file='residui'
  open(2,file=file,form='unformatted',status='old',iostat=ios)
  if(ios/=0) then
    write(stdout_unit,*)'cannot open ',file
  endif
  call fill_list (iunit=2)
  file='result1'
  open(3,file=file,form='unformatted',status='old',iostat=ios)
  if(ios/=0) then
    write(stdout_unit,*)'cannot open ',file
  endif
  call fill_list (iunit=3)
  file='result2'
  open(4,file=file,form='unformatted',status='old',iostat=ios)
  if(ios/=0) then
    write(stdout_unit,*)'cannot open ',file
  endif
  call fill_list (iunit=4)

  call showvar
  !----------------------------
  ! default values for plotting
  !----------------------------
  n    = 1
  do i=1,nvar
    if(l(i)% size > 2 .and. l(i)% n > 1) then
      n = i
      exit
    endif
  end do
  nlev = maxval(l(1:nvar)% size)
  k    = min(k,l(n)%size)
  t    = min(t,l(n)%n)
  !---------------
  ! loop for plots
  !---------------
  do
    select case (mode)
    case ('l')
      !-------------
      ! plot vs time
      !-------------
      allocate (x(l(n)%n,5))
      x = 0._wp
      do iunit = 1,4
        inquire(iunit,opened=opened)
        if (opened) rewind(iunit)
        j = 0
        do
          call readfile(z ,name ,iunit ,index ,lonlat)
          if(.not.associated(z)) then
            if(j>0 .and. j<size(x,1)) then
              write(0,*)'EOF unit ',iunit,&
                        ', set remaining values to',x(j,iunit)
              x(j+1:,iunit) = x(j,iunit)
            endif
            exit
          endif
          if(name==l(n)%name) then
            j=j+1
            x(j,iunit) = z(k)
          endif
        end do
      end do
      !-----------------------------------
      ! calculate tendency of column-model
      !-----------------------------------
      x(:,5)  = 0
      x(2:,5) = x(2:,3)-x(1:l(n)%n-1,4)
      if(all(x(:,4)==0)) x(2:,5) = x(2:,3)-x(1:l(n)%n-1,3)
      !--------------
      ! subtract mean
      !--------------
      if(lmin) then
        mean             = sum(x(t0:t9,1))/size(x(t0:t9,1))
        if(mean==0) mean = sum(x(t0:t9,3))/size(x(t0:t9,3))
        if(any(x(:,1)/=0)) x(:,1) = x(:,1) - mean    ! forcing
        if(any(x(:,3)/=0)) x(:,3) = x(:,3) - mean    ! result 1
        if(any(x(:,4)/=0)) x(:,4) = x(:,4) - mean    ! result 2
      else
        mean = 0._sp
      endif
      !-------------------
      ! write gnuplot file
      !-------------------
      open(10,file='gnu.plot')
        do j=1,l(n)%n
          write(10,*)j,real(x(j,:),sp)
        end do
      close(10)
      call gnuplot
      deallocate (x)
    case ('t')
      !-------------
      ! plot vs height
      !-------------
      allocate (x(l(n)%size,5))
      x = 0._wp
      do iunit=1,4
        inquire (iunit, opened=opened)
        if (opened) rewind(iunit)
        j = 0
        do
          call readfile(z ,name ,iunit ,index ,lonlat)
          if(.not.associated(z)) exit
          if(name==l(n)%name) then
            j=j+1
            if(j==t-1.and.iunit==4) x(:,5) = z
            if(j==t)then
              x(:,iunit) = z
              exit
            endif
          endif
        end do
        close(10)
      end do
      !-----------------------------------
      ! calculate tendency of column-model
      !-----------------------------------
      x(:,5) = x(:,3)-x(:,5)
      !-------------------
      ! write gnuplot file
      !-------------------
      open(10,file='gnu.plot')
        do j=1,l(n)%size
          write(10,*)j,real(x(j,:),sp)
        end do
      close (10)
      call gnuplot
      deallocate (x)
    case ('a')
      !------------------
      ! write ascii files
      !------------------
      open(11,file='forcing.asc')
      open(12,file='residui.asc')
      open(13,file='result1.asc')
      open(14,file='result2.asc')
      do iunit=1,4
        inquire (iunit, opened=opened)
        if (opened) then
          rewind(iunit)
          j = 0
          do
            if(j==t) then
              call readfile(z ,name ,iunit ,index ,lonlat &
                           ,ounit=iunit+10 ,ofmt=.true.)
            else
              call readfile(z ,name ,iunit ,index ,lonlat)
            endif
            if(.not.associated(z)) exit
            if(name=='NSTEP') then
              if(j==t) exit
              j=j+1
              if(j==t) then
                write(iunit+10,*) 'NSTEP            1 1 0 0'
                write(iunit+10,'(g14.7)') z
              endif
            endif
          end do
        end if
      end do
      close(11)
      close(12)
      close(13)
      close(14)
      mode = 't'
    end select
    !-------------------------------
    ! change variable, level, height
    !-------------------------------
    write(stdout_unit,*) 'variable =',trim(l(n)%name),',timestep=',t,&
      ',level=',k,',t0=',t0,',plot=',mode,' ',trim(ccol(icol))
    c = ask_menue('>'                                  ,&
      'v choose new variable'                          ,&
      'u level up   (decrese k)'                       ,&
      'd level down (increase k)'                      ,&
      'l choose new level'                             ,&
      'p previous time step'                           ,&
      'n next time step'                               ,&
      't choose new time step'                         ,&
      'a write ascii files '                           ,&
      '0 change t0 (first time step plotted)'          ,&
      '9 change t9 (last time step plotted)'           ,&
      '- subtract mean'                                ,&
      'e write eps file'                               ,&
      'c toggle color/mono'                            ,&
      'x chose new index (coordinates) in forcing file',&
      'q quit program'                                 ,&
      title='change plot parameters')
    select case (c)
    case ('x')
      index(1) = -1
      lonlat   = (/0,0/)
    case ('c')
      icol = 3-icol
    case ('e')
      call ask(epsfile,'eps file basename')
      if(gpu==6) then
        write(gpu,'(a,a)')   "set terminal postscript eps ",trim(ccol(icol))
        write(gpu,'(a,a,a)') "set output '",trim(epsfile),".eps'"
        write(gpu,'(a)') "replot"
        write(gpu,'(a)') "set output"
        write(gpu,'(a)') "set terminal x11"
        write(gpu,'(a)') "replot"
        call flush(gpu)
      endif
    case ('-')
      lmin = .not. lmin
    case ('0')
      call ask(t0,'first time step plotted')
      call checkpar
    case ('9')
      call ask(t9,'last time step plotted')
      call checkpar
    case ('v')
      call showvar
      call ask(n,'variable number')
      n=min(nvar,max(1,n))
      write(stdout_unit,"(i3,1x,a,2i3)")n,l(n)%name,l(n)%size,l(n)%n
      k = min(k,l(n)%size)
      t = min(t,l(n)%n)
    case ('u')
      k = max(1,k - 1)
      mode = 'l'
    case ('d')
      k = min(k+1,l(n)%size)
      mode = 'l'
    case ('l')
      call ask(k,'level')
      k = min(max(k,1),l(n)%size)
      mode = 'l'
    case ('n')
      t = min(t+1,l(n)%n)
      mode = 't'
    case ('t')
      call ask(t,'time')
      t = min(max(t,1),l(n)%n)
      mode = 't'
    case ('a')
      t = min(max(t,1),l(n)%n)
      mode = 'a'
    case ('p')
      t = max(t-1,1)            
      mode = 't'
    case ('q')
      exit
    end select
  end do
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  subroutine showvar
    integer :: i
    do i=1,nvar
      write(stdout_unit,"(i3,1x,a,2i6)")i,l(i)%name,l(i)%size,l(i)%n
    end do    
  end subroutine showvar
!------------------------------------------------------------------------------
  subroutine fill_list (iunit)
    integer ,intent(in)    :: iunit
    logical :: opened
    l% temp = 0
    inquire (iunit, opened=opened)
    if(opened) then
      rewind (iunit)
      do
        call readfile (z ,name ,iunit ,index ,lonlat)
        if(.not.associated(z)) exit
        do i=1,size(l)
          if(name==l(i)%name) then
            l(i)% temp = l(i)% temp + 1
            if (l(i)% size /= size(z)) then
              write(0,*)'post_1d: name  =',name
              write(0,*)'post_1d: iunit =',iunit
              write(0,*)'post_1d: read size =',size(z)
              write(0,*)'post_1d: exp. size =',l(i)% size
              write(0,*)'post_1d: incompatible size !'
              call abort('post_1d: incompatible size of '//name)
            endif
            exit
          endif
          if(l(i)%name=='') then
            l(i)% name = name
            l(i)% size = size(z)
            l(i)% temp = 1
            nvar=nvar+1
            exit
          endif
          if(i==size(l)) call abort('post_1d: increase variable list!')
        end do
      end do
      rewind (iunit)
    endif
    l% n = max(l% n, l% temp)
  end subroutine fill_list
!------------------------------------------------------------------------------
  subroutine gnuplot
    call checkpar
    if(gpu/=6) open(gpu,file='gnu.load')
    write (gpu,*) 'set data style lines'
    select case (mode)
    case ('l')
      if (t9<t0) then
        write(0,*) 'empty x range !'
        return
      endif
      write (range,"('[',i6,':',i6,']')") t0,t9
      write (gpu,*) "set title 'variable: "//trim(l(n)%name)//"  level:",k,&
                    "  offset:",mean,"'"
      write (gpu,*) "set xlabel 'timestep'"
      write (gpu,*) "set ylabel '"//trim(l(n)%name)//"'"
      write (gpu,*) "plot "//trim(range)                       ,&
                        "'gnu.plot' using 1:2 title '"//m(1)//"',",&
                        "'gnu.plot' using 1:3 title 'residui',",&
                        "'gnu.plot' using 1:4 title 'result1',",&
                        "'gnu.plot' using 1:5 title 'result2',",&
                        "'gnu.plot' using 1:6 title '1d_tend'"
    case ('t')
      write(gpu,*) "set title 'variable: "//trim(l(n)%name)//&
                                         "  timestep:",t,"'" 
      write(gpu,*) "set ylabel 'level'"
      write(gpu,*) "set xlabel '"//trim(l(n)%name)//"'"
      write(gpu,*) "plot [:] [",nlev,":0] ",                                  &
                                   "'gnu.plot' using 2:1 title '"//m(1)//"',",&
                                   "'gnu.plot' using 3:1 title 'residui',",   &
                                   "'gnu.plot' using 4:1 title 'result1',",   &
                                   "'gnu.plot' using 5:1 title 'result2',",   &
                                   "'gnu.plot' using 6:1 title '1d_tend'"
    end select
    if(gpu/=6) then
      close (gpu)
    else
      call flush(gpu)
    endif
  end subroutine gnuplot
!------------------------------------------------------------------------------
  subroutine checkpar
    t0=max(1 ,min(t0,l(n)%n))
    t9=min(max(t0+1,t9),l(n)%n)
  end subroutine checkpar
!------------------------------------------------------------------------------
end program post_1d
