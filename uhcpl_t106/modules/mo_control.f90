MODULE mo_control

  ! Control variables for model housekeeping.
  !
  ! U. Schlese, DKRZ, December 1994
  ! A. Rhodin, MPI, January 1999,
  !      Subroutine m_control renamed to alloc_mods and moved from
  !      module mo_control to module m_alloc_mods.
  ! L. Kornblueh, MPI, June 1999,
  !      added nproca and nprocb for driving the parallel decomposition
  ! I.Kirchner, IfM, July 2004,
  !      added rerun ensemble generation

  IMPLICIT NONE

  REAL :: dtime             !   time step (in seconds).
  REAL :: twodt             !   2.*dtime.
  REAL :: eps               !   time filtering coefficient.

  REAL, POINTER :: vct(:)   !   vertical coefficients table.


  INTEGER :: nproca         !   number of processors in set A
  INTEGER :: nprocb         !   number of processors in set A
  INTEGER :: nm             !   max zonal wave number.
  INTEGER :: nn             !   max meridional wave number for m=0.
  INTEGER :: nk             !   max meridional wave number.
  INTEGER :: ngl            !   number of gaussian latitudes.
  INTEGER :: nlon           !   max number of points on each latitude line.
  INTEGER :: nlev           !   number of vertical levels.
  INTEGER :: nmp1           !   max zonal wave number + 1.
  INTEGER :: nnp1           !   max meridional wave number + 1.
  INTEGER :: nkp1
  INTEGER :: n2mp1          !   2 * (max zonal wave number + 1).
  INTEGER :: n4mp1          !   4 * (max zonal wave number + 1).
  INTEGER :: nlp2           !   max number of points per latitude line + 2.
  INTEGER :: nlevp1         !   *nlev+1.
  INTEGER :: nsp            !   number of spectral coefficients.
  INTEGER :: n2sp           !   2*number of spectral coefficients.
  INTEGER :: nhgl           !   (number of gaussian latitudes)/2.
  INTEGER :: nscan          !   current scan number.
  INTEGER :: nresum         !   time step at which the run started or was
  !                             reumed after interruption.
  INTEGER :: ncbase         !   century date of the initial data.
  INTEGER :: ntbase         !   time of the initial data.
  INTEGER :: ncdata         !   verifying (century) date of the initial data.
  INTEGER :: ntdata         !   verifying time of the initial data.
  INTEGER :: ntimst         !   constant to convert *ntbase* into seconds.
  INTEGER :: nwtime         !   rerun write-up interval
  INTEGER :: nptime         !   array of post-processing times required
  INTEGER :: nctime         !   coupling intervall ocean
  INTEGER :: n4ptime
  INTEGER :: nwlag          !   rerun-files saving interval in months
  INTEGER :: nsub           !   number of jobs to be submitted at end of ru
  INTEGER :: nsubint(9)     !   submit-interval in months
  INTEGER :: nstop          !   last time step.
  INTEGER :: nrow(3)        !   current latitude line. (one entry per task).
  INTEGER :: maxrow         !   number of latitude lines.
  INTEGER :: nvclev         !   number of levels with vertical coefficients.
  INTEGER :: nlat(1)        !   current "geographic" latitude line
  !                             from north to south (one entry per task)
  INTEGER :: numfl1         !   number of optional fields read at nstep=0
  INTEGER :: numfl2         !   number of optional fields read at nstep=nresum
  INTEGER :: imtt           ! fu++ (1,1000) initial ocean index

  LOGICAL :: labort
  LOGICAL :: lwtime         !   .TRUE. when history writeup is due
  LOGICAL :: lptime         !   .TRUE. when postprocessing step is due
  LOGICAL :: lctime         !   .TRUE. when ocean coupling step is due
  LOGICAL :: lrepro         !   .TRUE. for reproducable results in multitsk.
  LOGICAL :: ldebug         !   .TRUE. for mass fixer diagnostics
  LOGICAL :: lg4x
  LOGICAL :: l4ptime
  LOGICAL :: lamip          !   .TRUE. for using variable sst
  LOGICAL :: lamip2         !   .TRUE. for using amip2 sst (additional seaice file) 
  LOGICAL :: lsub           !   .TRUE. to submit *nsub* jobs
  LOGICAL :: lsstadj        !   .TRUE. for orographic adjustment of sst
  LOGICAL :: lhg3x          !   .TRUE. to use g3x info from history files
  LOGICAL :: lcouple        !   .TRUE. for a coupled run
  LOGICAL :: lonudg         !   .TRUE. to nudge ocean
  LOGICAL :: ldsst          !   .TRUE. turn-on diurnal SST fu++
  LOGICAL :: lstratiform    !   .TRUE. turn-on diurnal SST fu++
  LOGICAL :: lnwp           !   .FALSE. for climate mode .true. for NWP mode
  LOGICAL :: lanalysis      !   .FALSE. for climate mode .true. for analysis
  LOGICAL :: lnudge         !   .TRUE. for Nudging mode
  LOGICAL :: lmidatm        !   .TRUE. for middle atmosphere model version

  LOGICAL :: lnmi           !   .TRUE. normal mode initialisation
  LOGICAL :: ltdiag         !   .TRUE. run with additional diagnostics of tendency terms
  LOGICAL :: lcolumn        !   .TRUE. for column model
  LOGICAL :: lvctch         !   .TRUE. if vct changed in column model

  LOGICAL :: lpci           !   .TRUE. for cond5

  LOGICAL :: lalai          !   .TRUE. for annual lai
  LOGICAL :: lavgrat        !   .TRUE. for annual vegetation

  LOGICAL :: lens           !   .TRUE. generate rerun ensemble
  INTEGER :: nens           !    number of ensembles

END MODULE mo_control
