SUBROUTINE init_memory

  USE mo_doctor,        ONLY: nout
  USE mo_memory_base,   ONLY: print_memory_use, print_sinfo
  USE mo_memory_sp,     ONLY: sp,  construct_sp
  USE mo_memory_ls,     ONLY:      construct_ls
  USE mo_memory_gl,     ONLY: gl,  construct_gl
  USE mo_memory_f,      ONLY: f,   construct_f
  USE mo_memory_g1a,    ONLY: g1a, construct_g1a
  USE mo_memory_g1b,    ONLY: g1b, construct_g1b
  USE mo_memory_g2a,    ONLY: g2a, construct_g2a
  USE mo_memory_g2b,    ONLY: g2b, construct_g2b
  USE mo_memory_g3a,    ONLY: g3a, construct_g3a
  USE mo_memory_g3b,    ONLY: g3b, construct_g3b
  USE mo_buffer_fft,    ONLY: construct_fft
  USE mo_tracer,        ONLY: ntrac
  USE mo_control,       ONLY: ngl, nhgl, nlev, nlon, nmp1, nnp1, nsp
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_start_dataset, ONLY: ldebugmem
  USE mo_mpi,           ONLY: p_pe, p_io

  IMPLICIT NONE


  ! allocate global description arrays

  CALL construct_sp  (nlev,       lc%nsnm0, lc%snsp, &
                      nlev,       nnp1,     nsp)

  CALL construct_ls  (lc%nllev, lc%nllevp1, lc%nlnm0, lc%lnsp, &
                      nlev,                 nnp1,     nsp)

  CALL construct_f   (lc%nllev, lc%nllevp1, lc%nlm,  lc%nlat/2, &  
                      nlev,                 nmp1,    nhgl)

  CALL construct_gl  (lc%nglon,   lc%nglpx, lc%nlev, ntrac, lc%nglat, &
                      nlon,                 nlev,    ntrac, ngl)

  CALL construct_g1a (lc%nglon,   lc%nglpx, lc%nlev, ntrac, lc%nglat, &
                      nlon,                 nlev,    ntrac, ngl)
  CALL construct_g1b (lc%nglon,   lc%nglpx, lc%nlev, ntrac, lc%nglat, &
                      nlon,                 nlev,    ntrac, ngl)

  CALL construct_g2a (lc%nglon,   lc%nglpx, lc%nlev, lc%nglat, &
                      nlon,                    nlev,    ngl)
  CALL construct_g2b (lc%nglon,   lc%nglpx, lc%nlev, lc%nglat, &
                      nlon,                    nlev,    ngl)

  CALL construct_g3a (lc%nglon,   lc%nglpx, lc%nlev, lc%nglat, &
                      nlon,                    nlev,    ngl)
  CALL construct_g3b (lc%nglon,   lc%nglpx, lc%nlev, lc%nglat, &
                      nlon,                 nlev,    ngl)

  CALL construct_fft (lc)

  IF (ldebugmem) THEN
     WRITE (nout, '(/,a)') ' Global memory buffers:'
     WRITE (nout, '(a)') '    sp-buffer:  '
     CALL print_sinfo(sp)
     WRITE (nout, '(a)') '     f-buffer:  '
     CALL print_sinfo(f)
     WRITE (nout, '(a)') '    gl-buffer:  '
     CALL print_sinfo(gl)
     WRITE (nout, '(a)') '   g1a-buffer:  '
     CALL print_sinfo(g1a)
     WRITE (nout, '(a)') '   g1b-buffer:  '
     CALL print_sinfo(g1b)
     WRITE (nout, '(a)') '   g2a-buffer:  '
     CALL print_sinfo(g2a)
     WRITE (nout, '(a)') '   g2b-buffer:  '
     CALL print_sinfo(g2b)
     WRITE (nout, '(a)') '   g3a-buffer:  '
     CALL print_sinfo(g3a)
     WRITE (nout, '(a)') '   g3b-buffer:  '
     CALL print_sinfo(g3b)
  END IF

  IF (p_pe == p_io) THEN
     WRITE (nout, '(/,a)') ' Global memory buffers:'
     WRITE (nout, '(a)', ADVANCE='NO') '    sp-buffer:  '
     CALL print_memory_use(sp)
     WRITE (nout, '(a)', ADVANCE='NO') '     f-buffer:  '
     CALL print_memory_use(f)
     WRITE (nout, '(a)', ADVANCE='NO') '    gl-buffer:  '
     CALL print_memory_use(gl)
     WRITE (nout, '(a)', ADVANCE='NO') '   g1a-buffer:  '
     CALL print_memory_use(g1a)
     WRITE (nout, '(a)', ADVANCE='NO') '   g1b-buffer:  '
     CALL print_memory_use(g1b)
     WRITE (nout, '(a)', ADVANCE='NO') '   g2a-buffer:  '
     CALL print_memory_use(g2a)
     WRITE (nout, '(a)', ADVANCE='NO') '   g2b-buffer:  '
     CALL print_memory_use(g2b)
     WRITE (nout, '(a)', ADVANCE='NO') '   g3a-buffer:  '
     CALL print_memory_use(g3a)
     WRITE (nout, '(a)', ADVANCE='NO') '   g3b-buffer:  '
     CALL print_memory_use(g3b)
  END IF

END SUBROUTINE init_memory
