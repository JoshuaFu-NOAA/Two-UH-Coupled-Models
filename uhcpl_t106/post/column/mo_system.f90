module mo_system
#ifdef NAG
  use f90_unix,      only: abort   ! unix abort call
  use f90_unix,      only: flush   ! flush output buffer
  use f90_unix_proc, only: system  ! unix system call
  use f90_unix_env,  only: getarg, ! get program arguments
                           iargc
#else
  external iargc
  integer  iargc
#endif
end module mo_system
