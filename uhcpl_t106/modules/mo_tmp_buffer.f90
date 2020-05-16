MODULE mo_tmp_buffer

! Physics internals
                                        ! allocate: size:           deallocate:
REAL, ALLOCATABLE :: ahfli(:)           ! ice.f     (nlp2)           physc.f
REAL, ALLOCATABLE :: thfl(:)            ! physc.f   (nlp2)           skintem.f
REAL, ALLOCATABLE :: qhfl(:)            ! physc.f   (nlp2)           surf.f
REAL, ALLOCATABLE :: dhft(:)            ! physc.f   (nlp2)           skintem.f
REAL, ALLOCATABLE :: dhfq(:)            ! physc.f   (nlp2)          physc.f
REAL, TARGET, ALLOCATABLE :: dhfqw(:)   ! physc.f   (nlp2)           surf.f
REAL, TARGET, ALLOCATABLE :: dhfqs(:)   ! physc.f   (nlp2)           surf.f
REAL, ALLOCATABLE :: cvs(:)             ! physc.f   (nlp2)           surf.f
REAL, ALLOCATABLE :: cvw(:)             ! physc.f   (nlp2)           surf.f
REAL, ALLOCATABLE :: wlmx(:)            ! physc.f   (nlp2)           surf.f
REAL, ALLOCATABLE :: geom1(:,:)         ! physc.f   (nlp2,nlev)      physc.f
REAL, ALLOCATABLE :: aphm1(:,:)         ! physc.f   (nlp2,nlevp1)    physc.f
REAL, ALLOCATABLE :: apm1(:,:)          ! physc.f   (nlp2,nlev)      physcl.f
REAL, ALLOCATABLE :: aphp1(:,:)         ! physc.f   (nlp2,nlevp1)    physc.f
REAL, ALLOCATABLE :: app1(:,:)          ! physc.f   (nlp2,nlev)      physc.f
REAL, ALLOCATABLE :: srfl(:)            ! physc.f   (nlp2)           skintem.f
REAL, ALLOCATABLE :: rsfc(:)            ! physc.f   (nlp2)           physc.f
REAL, ALLOCATABLE :: ssfc(:)            ! physc.f   (nlp2)           physc.f
REAL, ALLOCATABLE :: rsfl(:)            ! physc.f   (nlp2)           surf.f
REAL, ALLOCATABLE :: ssfl(:)            ! physc.f   (nlp2)           surf.f
REAL, ALLOCATABLE :: xhfl(:)            ! physc.f   (nlp2)           physc.f
REAL, ALLOCATABLE :: xtec(:,:)          ! physc.f   (nlp2,nlev)      physc.f

! Matrix for Helmholtz-equation

REAL, ALLOCATABLE :: cn(:,:,:)         ! setdyn.f  (nlev,nlev,nkp1) -
REAL, ALLOCATABLE :: dl(:,:)           ! si2.f     (nlp2,nlev)      scan1sl.f
REAL, ALLOCATABLE :: dm(:,:)           ! si1.f     (nlp2,nlev)      scan1sl.f

! Used in first scan

!REAL, ALLOCATABLE :: aph(:,:)           ! scan1sl.f (nlp2,nlevp1)    statd.f
!REAL, ALLOCATABLE :: delp(:,:)          ! scan1sl.f (nlp2,nlev)      dyn.f

! Can be a fixed array, there's no need to make it allocatable

REAL :: dia(56)                         ! in prestat.F and postatp.F

! Radiation only, no need to deallocate ...

REAL, ALLOCATABLE :: diag(:,:)          ! setrad.f  (4,nlevp1)       -
REAL :: dia1(7)	                        ! setrad.f  (7)              -

LOGICAL, ALLOCATABLE :: loland(:)       ! physc.f   (nlp2)           physc.f
LOGICAL, ALLOCATABLE :: loglac(:)       ! physc.f   (nlp2)           physc.f

END MODULE mo_tmp_buffer


