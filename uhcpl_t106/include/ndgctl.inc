!*    NAMELIST *NDGCTL* CONTROLS NUDGING OPTIONS
!
!       I. Kirchner MPI Hamburg  Nov 1998
!       I. Kirchner MPI Hamburg  Oct 1999
!       I. Kirchner MPI Hamburg  Mar 2000
!       I. Kirchner MPI Hamburg  Sep 2000
!       I. Kirchner MPI Hamburg  May 2001
!	I. Kirchner IfM Berlin   Sep 2005
!
!       parameters of modules/mo_nudging.f90
!
NAMELIST /ndgctl/&
&     lnudgdbx, & ! LOGICAL        debugging message control
&     lnudgini, & ! LOGICAL        adjust date/time from nudging data
&     lnudgimp, & ! LOGICAL        choose implicit/explicit methode
&     lnudgpat, & ! LOGICAL        assimilation of correlation pattern
&     lnudgcli, & ! LOGICAL        assimilate climatological data
&     lnudgfrd, & ! LOGICAL        place for NMI filter
&     lsite,    & ! LOGICAL        calculate SITE
&     lnudgwobs,& ! LOGICAL        store reference fields
&     ldamplin, & ! LOGICAL        damping type time interpolation
&     ltintlin, & ! LOGICAL        time interpolation method
&     nudgd,    & ! REAL    (nlev) nudging weight for divergence
&     nudgv,    & ! REAL    (nlev) nudging weight for vorticity
&     nudgt,    & ! REAL    (nlev) nudging weight for temperature
&     nudgp,    & ! REAL           nudging weight for log sfc pressure
&     nudgdamp, & ! REAL    0-1    damping between two nudging data steps
&     nudgdsize,& ! REAL    0-0.5  nudging radius near nudging time
&     nudgtrun, & ! INTEGER 0|1|2  type of wave index selection
&     nudgsmin, & ! INTEGER        lowest modified wave number
&     nudgsmax, & ! INTEGER        highest modified wave number
&     nudglmin, & ! INTEGER        lowest level index (top->down)
&     nudglmax, & ! INTEGER        highest level index
&     nudg_start,&! INTEGER        YYYYMMDD start with nudging (default -1)
&     nudg_stop,& ! INTEGER        YYYYMMDD stop with nudging (default -1)
&     nsstinc,  & ! INTEGER        hours between new SST fields
&     nsstoff,  & ! INTEGER        offset in hours to 00UTC
&     ndg_freez,& ! REAL           freezing point for external data used
&     ndgfmt,&    ! CHARACTER(4)   data format type
&     ndgswap,&   ! LOGICAL        activate swap32
&     ndunit      ! INTEGER (13)   input/output channels
!
! version 26-september-2005
