
#ifdef CRAY
#define ddot sdot
#endif

MODULE mo_tracer

  USE mo_parameters, ONLY: jps, jptrac, jpgl

  IMPLICIT NONE

  INTEGER :: ntrac
  INTEGER :: ntraca
  INTEGER :: nhtrac
  INTEGER :: ntcode(jptrac)
  INTEGER :: nfix(jps)
  INTEGER :: nfixt(jptrac)
  INTEGER :: ifileini
  INTEGER :: ifileemi

  LOGICAL :: lppxt(jptrac)
  LOGICAL :: lhtrac
  LOGICAL :: lslt(jptrac)
  LOGICAL :: lxtvdiff
  LOGICAL :: lxtconv
  LOGICAL :: lwetdep(jptrac)
  LOGICAL :: lemis

  REAL    :: sxtini(jptrac)
  REAL    :: sxtemi(jptrac)
  REAL    :: vdrydep(jptrac)
  REAL    :: sxtsink(jptrac)

  INTEGER :: nstr
  REAL, TARGET, ALLOCATABLE :: source(:,:,:)               ! (nlon,ngl,nstr)

  !
  !      zonal mass budgets of tracers
  !

  REAL    :: tropm(jpgl,jptrac+1)
  REAL    :: stratm(jpgl,jptrac+1)

CONTAINS

  SUBROUTINE initrac

    ! Description:
    !
    ! Preset constants in tractl.
    !
    ! Method:
    !
    ! Preset and read the namelist *tractl*.
    !
    ! *inictl* is called from *initialise*.
    !
    ! Authors:
    !
    ! J. Feichter, MI, September 1990, original source
    ! U. Hansson, MI, July 1991, changed
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_mpi,            ONLY: p_parallel, p_parallel_io,  p_io, p_bcast
    USE mo_doctor,         ONLY: nin
    USE mo_start_dataset,  ONLY: lres

    !  Local scalars: 
    INTEGER :: icode, itraca, j, jt

    INCLUDE 'tractl.inc'


    !  Executable statements 

    !-- 1. Preset namelist variables

    ntrac  = 0
    lhtrac = .FALSE.
    icode  = 234

    DO j = 1, jps
       nfix(j) = 1
    END DO
    DO jt = 1, jptrac
       lppxt(jt)   = .FALSE.
       lslt(jt)    = .TRUE.
       nfixt(jt)   = 1
       icode       = icode + 1
       ntcode(jt)  = icode
       sxtini(jt)  = 0.
       sxtemi(jt)  = 0.
       vdrydep(jt) = 0.
       sxtsink(jt) = 0.
    END DO

    lxtvdiff = .TRUE.
    lxtconv  = .TRUE.

    !-- 2. Read namelist tractl

    IF (p_parallel) THEN
       IF (p_parallel_io) THEN
          READ (nin,tractl)
       ENDIF
       CALL p_bcast (ntrac,    p_io)
       CALL p_bcast (lslt,     p_io)
       CALL p_bcast (lppxt,    p_io)
       CALL p_bcast (lhtrac,   p_io)
       CALL p_bcast (ntcode,   p_io)
       CALL p_bcast (nfix,     p_io)
       CALL p_bcast (nfixt,    p_io)
       CALL p_bcast (lxtvdiff, p_io)
       CALL p_bcast (lxtconv,  p_io)
       CALL p_bcast (lwetdep,  p_io)
       CALL p_bcast (sxtini,   p_io)
       CALL p_bcast (sxtemi,   p_io)
       CALL p_bcast (vdrydep,  p_io)
       CALL p_bcast (sxtsink,  p_io)
    ELSE
       READ (nin,tractl)
    ENDIF

    ! Choose number of tracers in case of rerun

    itraca = 0
    DO jt = 1, ntrac
       IF (lslt(jt)) itraca = itraca + 1
       IF ( .NOT. lslt(jt)) nfixt(jt) = 0
    END DO

    ntraca = itraca

    IF (ntrac == 0) THEN
       lxtvdiff = .FALSE.
       lxtconv  = .FALSE.
    END IF

    IF (.NOT. lres) nhtrac = ntrac

  END SUBROUTINE initrac

  SUBROUTINE xtini(kt1, xt)

    ! Description:
    !
    ! Initialize the tracer fields
    !
    ! Method:
    !
    ! This routine is called from *ioinitial*
    ! and from *iorestart* in case the run is continued
    ! from a run which contained less tracers
    !
    ! Authors:
    !
    ! J. Feichter, MI, August 1991, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_start_dataset,   ONLY: nini
    USE mo_doctor,          ONLY: nout
    USE mo_exception,       ONLY: finish

    !  Array arguments 
    REAL :: xt(:,:,:,:)

    !  Scalar arguments
    INTEGER :: kt1

    !  Local scalars:
    INTEGER ::  jt


    !  Executable statements

    DO jt = kt1, ntrac
       IF (NINT(sxtini(jt)) /= -999) THEN
          xt(:,:,jt,:) = sxtini(jt)
          WRITE(nout,*) ' set tracer ', jt, ' to ', sxtini(jt)
       ELSE
          !        READ (nini,err=10,end=10) (pgbuf(jle),jle=i1,i2)
          CALL finish ('xtini','initial tracer from file not implemented')
       END IF
    END DO

    RETURN

10  CONTINUE
    WRITE (nout,'(A,I4,A)') ' end of file or bad format on unit ', &
                              nini,' =tracer initial fields'

    RETURN
  END SUBROUTINE xtini

  SUBROUTINE xtsink ( klon, klev, ptmst, pxtm1, pxte )

    ! Description:
    !
    ! Calculates the decrease of tracer concentration for a given halflife-time
    !
    ! Method:
    !
    ! The mass mixing-ratio of tracers is multiplied with
    ! exp(log(0.5)*time-step/half-life-time).
    ! This routine could also be used for emission or sink
    ! above the surface.
    !
    ! *xtsink* is called from *physc*
    !
    ! Authors:
    !
    ! J. Feichter, MI, August 1991, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS 
    !

    !  Scalar arguments 
    REAL :: ptmst
    INTEGER :: klev, klon

    !  Array arguments 
    REAL :: pxte(klon,klev,ntrac), pxtm1(klon,klev,ntrac)

    !  Local scalars: 
    REAL :: pqtmst, zdecay, zdxtdt, zfac, zxtp1
    INTEGER :: jk, jl, jt

    !  Intrinsic functions 
    INTRINSIC EXP, LOG


    !  Executable statements 

    pqtmst = 1./ptmst
    zfac   = LOG(0.5)*ptmst

    DO jt = 1, ntrac
       IF (ABS(sxtsink(jt)) > 0.) THEN
          zdecay = EXP(zfac/sxtsink(jt))
          DO jk = 1, klev
             DO jl = 1, klon
                zxtp1  = pxtm1(jl,jk,jt) + pxte(jl,jk,jt)*ptmst
                zxtp1  = zxtp1*zdecay
                zdxtdt = (zxtp1-pxtm1(jl,jk,jt))*pqtmst - pxte(jl,jk,jt)
                pxte(jl,jk,jt) = pxte(jl,jk,jt) + zdxtdt
             END DO
          END DO
       END IF
    END DO

    RETURN
  END SUBROUTINE xtsink

  SUBROUTINE xtemiss ( klon,   klev,    krow,    pcvdifts,  pdtime, &
                       pxtm1,  pxtems,  p1mxtm1,                    &
                       loland, pforest, psnow )

    ! Description:
    !
    ! Calculates the lower boundary conditions for vdiff depending on
    ! the surface emission and the dry deposition flux.
    !
    ! Method:
    !
    ! The lower boundary condition for calculating the
    ! turbulent exchange in the boundary layer is
    ! determined by the emission and the dry deposition flux.
    !
    ! *xtemiss* is called from *vdiff*
    !
    ! Authors:
    !
    ! J. Feichter, MI, August 1991, original source
    ! U. Schlese, DKRZ, January 1995, changed
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    !  Scalar arguments 
    REAL :: pcvdifts, pdtime
    INTEGER :: klev, klon, krow

    !  Array arguments 
    REAL :: p1mxtm1(klon), pforest(klon+2), psnow(klon+2), &
            pxtems(klon,ntrac), pxtm1(klon,klev,ntrac)
    LOGICAL :: loland(klon+2)

    !  Local scalars: 
    REAL :: zfacemi, zzvd
    INTEGER :: istr, jl, jt


    !  Executable statements 

    !-- 1. Surface emission

    zfacemi = 1.
    istr = 0
    DO jt = 1, ntrac
       IF (ABS(sxtemi(jt)) > 0.) THEN
          IF (NINT(sxtemi(jt)) == -999) THEN
             istr = istr + 1
             DO jl = 1, klon
                pxtems(jl,jt) = pcvdifts*source(jl,krow,istr)*zfacemi
             END DO
          ELSE
             DO jl = 1, klon
                pxtems(jl,jt) = pcvdifts*sxtemi(jt)*zfacemi
             END DO
          END IF
       END IF
    END DO

    !-- 2. Dry deposition

    DO jt = 1, ntrac
       IF (ABS(vdrydep(jt)) > 0.) THEN
          zzvd = vdrydep(jt)*0.01
          DO jl = 1, klon
             pxtems(jl,jt) = pxtems(jl,jt) &
                           - pcvdifts*p1mxtm1(jl)*pxtm1(jl,klev,jt)*zzvd
          END DO
       END IF
    END DO

    !-- 3. Average mixing-ratios of the lowest model layer
    !
    !       DO jl=1,klon
    !        g3x(jl,1,nn)=g3xm(jl,1,nn)+pxtm1(jl,klev,jt)*pdtime
    !       END DO

    RETURN
  END SUBROUTINE xtemiss

  SUBROUTINE prestatr

    ! Description:
    !
    ! Clears tracer mass budget diagnostics
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, June 1994, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    !  Local scalars: 
    INTEGER :: jg, jt


    !  Executable statements 

    DO jt = 1, jptrac + 1
       DO jg = 1, jpgl
          tropm(jg,jt)  = 0.
          stratm(jg,jt) = 0.
       END DO
    END DO

    RETURN
  END SUBROUTINE prestatr

  SUBROUTINE trastat

    ! Description:
    !
    ! Prints out accumulated mass budgets for tracers at
    ! the end of a run

    USE mo_control,        ONLY: nresum, ngl
    USE mo_start_dataset,  ONLY: nstep
    USE mo_doctor,         ONLY: nout

    !  Local scalars: 
    REAL :: zmglob, zmnhk, zmshk, zmstrat, zmtrop, zqcount
    INTEGER :: icount, jt

    !  Local arrays: 
    REAL :: zmstratn(ntrac+1), zmstrats(ntrac+1), zmtropn(ntrac+1), &
            zmtrops(ntrac+1)

    !  Intrinsic functions 
    INTRINSIC REAL, SUM


    !  Executable statements 

    icount  = nstep - nresum + 1
    zqcount = 1./(REAL(icount))
    WRITE (nout,'(a,/,a,/,a,/,a)') &
         ' Tracer mass budget:', &
         ' -----------------------------------------------------------------', &
         ' Averaged mass budgets in [kg] ', &
         '            global     n-hem    s-hem    tropo    strat    n-tro    s-tro    n-str    s-str '

    DO jt = 1, ntrac + 1
       zmtropn(jt)  = SUM(tropm(1:ngl:2,jt))  * zqcount
       zmstratn(jt) = SUM(stratm(1:ngl:2,jt)) * zqcount
       zmtrops(jt)  = SUM(tropm(2:ngl:2,jt))  * zqcount
       zmstrats(jt) = SUM(stratm(2:ngl:2,jt)) * zqcount
    END DO

    DO jt = 1, ntrac + 1
       zmnhk   = zmtropn(jt)  + zmstratn(jt)
       zmshk   = zmtrops(jt)  + zmstrats(jt)
       zmtrop  = zmtropn(jt)  + zmtrops(jt)
       zmstrat = zmstratn(jt) + zmstrats(jt)
       zmglob  = zmnhk + zmshk

       IF (jt <= ntrac) THEN
          WRITE (nout,'(a,i2,9e9.2)') ' Tracer: ', jt, zmglob, zmnhk, zmshk, &
               zmtrop, zmstrat, zmtropn(jt), zmtrops(jt), zmstratn(jt), zmstrats(jt)
       ELSE
          WRITE (nout,'(a   ,9e9.2)') ' Air mass: ', zmglob, zmnhk, zmshk, &
               zmtrop, zmstrat, zmtropn(jt), zmtrops(jt), zmstratn(jt), zmstrats(jt)
       END IF
    END DO

    RETURN
  END SUBROUTINE trastat

  SUBROUTINE xttropo ( klon, klp2, klev, krow, ktrac,   &
                       ptm1,  papm1,  paphm1,  pgeom1,  &
                       pxtm1, ktropo )

    ! Description:
    !
    ! Calculates the tropopause height
    !
    ! Method:
    !
    ! *xttropo* determines the highest sigma-level in the troposphere
    ! "ktropo"
    ! depending on the vertical gradient of the potential temperature
    ! and calculates zonal mass-budgets of the different tracers
    !
    ! *xttropo* is called from *physc*
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, September 1993, original source
    ! U. Schlese, DKRZ, June 1994, changed
    ! J. Feichter, MI, unknown, changed
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    !
    ! for more details see file AUTHORS
    !

    USE mo_constants,         ONLY: cpd, g, rd
    USE mo_gaussgrid,         ONLY: budw

    !  Scalar arguments 
    INTEGER :: klev, klon, klp2, krow, ktrac

    !  Array arguments 
    REAL :: paphm1(klp2,klev+1), papm1(klp2,klev), pgeom1(klp2,klev), &
            ptm1(klp2,klev), pxtm1(klon,klev,ktrac)
    INTEGER :: ktropo(klon)

    !  Local scalars: 
    REAL :: zdtheta, zdz, zfac, zkappa, zlimdthdz, zqg, zthetan, zthetap
    INTEGER :: itop, itopm1, jk, jl, jt, ktp1

    !  Local arrays: 
    REAL :: zma(klon,klev), zmt(klon,klev,ktrac), zsfac(klon,klev), &
            ztfac(klon,klev)

    !  External functions 
    REAL, EXTERNAL :: ddot


    !  Executable statements 

    zqg       = 1./g
    zfac      = budw(krow) * 510.0644719E+12 * zqg
    itop      = klev - 5
    itopm1    = itop - 1
    zlimdthdz = 0.28E-02
    zkappa    = rd/cpd

    ktp1      = ktrac + 1

    DO jk = itopm1, 2, -1
       DO jl = 1, klon
          zthetap = ptm1(jl,jk+1) * (1000./papm1(jl,jk+1))**zkappa
          zthetan = ptm1(jl,jk-1) * (1000./papm1(jl,jk-1))**zkappa
          zdz     = (pgeom1(jl,jk-1) - pgeom1(jl,jk+1))*zqg
          zdtheta = (zthetan - zthetap)/zdz
          IF (zdtheta < zlimdthdz) THEN
             ktropo(jl) = jk
          END IF
       END DO
    END DO

    DO jk = 1, klev
       DO jl = 1, klon
          zma(jl,jk) = (paphm1(jl,jk+1) - paphm1(jl,jk))*zfac
       END DO
    END DO

    DO jt = 1, ktrac
       DO jk = 1, klev
          DO jl = 1, klon
             zmt(jl,jk,jt) = zma(jl,jk) * pxtm1(jl,jk,jt)
          END DO
       END DO
    END DO

    DO jk = 1, klev
       DO jl = 1, klon
          IF (jk >= ktropo(jl)) THEN
             ztfac(jl,jk) = 1.
             zsfac(jl,jk) = 0.
          ELSE
             ztfac(jl,jk) = 0.
             zsfac(jl,jk) = 1.
          END IF
       END DO
    END DO

    DO jt = 1, ktrac
       tropm(krow,jt)  = tropm(krow,jt) &
                       + ddot(klon*klev, zmt(:,:,jt), 1, ztfac, 1)
       stratm(krow,jt) = stratm(krow,jt) &
                       + ddot(klon*klev, zmt(:,:,jt), 1, zsfac, 1)
    END DO

    tropm(krow,ktp1)  = tropm(krow,ktp1)  + ddot(klon*klev, zma, 1, ztfac, 1)
    stratm(krow,ktp1) = stratm(krow,ktp1) + ddot(klon*klev, zma, 1, zsfac, 1)

    RETURN
  END SUBROUTINE xttropo

  SUBROUTINE reademi

    ! Description:
    !
    ! *reademi* reads clobal surface emission fields
    ! to be used in subroutine *xtemiss*
    !
    ! Method:
    !
    ! *reademi* is called from *control*
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, January 1995, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

    USE mo_exception,     ONLY: finish
    USE mo_control,       ONLY: nlon, ngl
    USE mo_start_dataset, ONLY: nemi
    USE mo_doctor,        ONLY: nout

    !  Local scalars: 
    INTEGER :: jl, jr, js, jt
    LOGICAL :: lex


    !  Executable statements 

    !-- 1. Allocate memory and read fields

    !-- 1.1   Determine how many fields are to be read

    lemis = .FALSE.
    nstr = 0
    DO jt = 1, ntrac
       IF (ABS(sxtemi(jt))>0. .OR. ABS(vdrydep(jt))>0.) THEN
          lemis = .TRUE.
          IF (NINT(sxtemi(jt)) == -999) THEN
             nstr = nstr + 1
          END IF
       END IF
    END DO

    IF (nstr==0) RETURN

    !-- 1.2   Allocate memory for global emission fields

    ALLOCATE (source(nlon,ngl,nstr))

    !-- 1.3    Open file and read

    INQUIRE (nemi,exist=lex)

    IF (lex) THEN
       OPEN (nemi,form='UNFORMATTED')
       WRITE (nout,*) ' Reading ', nstr, ' surface emission fields from', &
                      ' unit ', nemi
    ELSE
       WRITE (nout,*) ' Trying to read from unit ', nemi
       WRITE (nout,*) ' File for surface emissions does not exist.'
       CALL finish('reademi','Run terminated.')
    END IF

    DO js = 1, nstr
       READ (nemi) ((source(jl,jr,js),jl=1,nlon),jr=1,ngl,2), &
                   ((source(jl,jr,js),jl=1,nlon),jr=ngl,2,-2)
    END DO
  
    CLOSE (nemi)

    RETURN
  END SUBROUTINE reademi

END MODULE mo_tracer
