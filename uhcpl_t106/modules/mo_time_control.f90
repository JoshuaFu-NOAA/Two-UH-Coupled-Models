MODULE mo_time_control

  ! set NSTEP dependend switches
  !
  ! author:
  ! I.Kirchner, MPI Hamburg, Nov-99

  ! get variables from ECHAM environment

  USE mo_year,          ONLY: cd2dat, im2day
  USE mo_control,       ONLY: dtime, ntbase, ncbase
  USE mo_start_dataset, ONLY: nstep, ntimeadj
  USE mo_constants,     ONLY: dayl

  IMPLICIT NONE

  ! definition of increment interval types

  INTEGER, PARAMETER :: TIME_INC_DAYS=1          ! increment in days
  INTEGER, PARAMETER :: TIME_INC_HOUR=2          ! increment in hours

CONTAINS

  !-----------------------------------------------------------
  ! control the switches at the beginning of a time step
  !
  ! step_inc = 0 --> .true. for the last time time step in
  !                  present month, given by my_step
  !
  ! step_inc < 0 --> increment in full inc_type intervalls
  ! step_inc > 0 --> increment in time steps
  !                  .TRUE. if the time from the
  !                  beginning of ncbase (00UTC) is a
  !                  multiple of the time given with step_inc
  !
  FUNCTION time_control (step_inc, inc_type, my_step) RESULT (action)

    INTEGER, INTENT (in) :: &
         step_inc,    &! increment between the action
         inc_type,    &! interpretation of negative values
         my_step       ! time step, at which the action 
                       ! should be controled
    LOGICAL :: action

    INTEGER :: day1, mon1, year, day2, mon2, lday
    INTEGER :: count, steps, ncdays


    action = .FALSE.

    IF (step_inc < 0) THEN
       ! negative increment means step_inc counts
       ! full days/hours instead of time steps
       SELECT CASE (inc_type)
       CASE (TIME_INC_DAYS); count = -step_inc*(dayl+0.01)/dtime
       CASE (TIME_INC_HOUR); count = -step_inc*(3600+0.01)/dtime
       CASE default;         count = 0
       END SELECT
    ELSE
       ! positive means the step_inc is given in time steps
       count = step_inc
    END IF

    IF (count > 0) THEN
       ! take the offset relative to 00UTC + ntimeadj into the acount
       IF (ntimeadj >= 0) THEN
         steps = (ntbase + ntimeadj + my_step*dtime + 0.01)/dtime
       ELSE
         steps = my_step
       END IF
       IF (MOD(steps,count) == 0) action = .TRUE.
    ELSE
       ! count equals zero set the action only
       ! true at the last time step in the present months
       ! get the date of the present time step (asking for)
       ncdays = ncbase + (ntbase+dtime*my_step+0.01)/dayl
       CALL cd2dat(ncdays,day1,mon1,year)
       ! get the last day of the present months
       lday = im2day(mon1,year)
       ! get the date of the next time step
       ncdays = ncbase + (ntbase+dtime*(my_step+1)+0.01)/dayl
       CALL cd2dat(ncdays,day2,mon2,year)
       IF ( (lday == day1) &! present day is 
            .AND.          &! the last day of the month
            (mon1 /= mon2) &! the next time step is at 
            ) THEN          ! the first day of the next months
          action = .TRUE.
       END IF
    END IF

  END FUNCTION time_control

  ! calculate the length of the interval given by
  ! step_in in secondes
  !
  FUNCTION time_inc_sec(step_inc, inc_type) RESULT (inc_sec)

    INTEGER, INTENT(in) :: step_inc, inc_type
    INTEGER :: inc_sec, iday, id, im, iy, nmend

    IF (step_inc < 0) THEN
       ! negative increment means step_inc counts
       ! full days/hours instead of time steps
       SELECT CASE (inc_type)
       CASE (TIME_INC_DAYS); inc_sec = -step_inc*INT(dayl)
       CASE (TIME_INC_HOUR); inc_sec = -step_inc*3600
       CASE default;         inc_sec = 0
       END SELECT
    ELSE IF (step_inc > 0) THEN
       ! positive means the step_inc is given in time steps
       inc_sec = step_inc*dtime
    ELSE
       ! the zero case
       iday = ncbase + (ntbase+dtime*nstep)/dayl + 0.01
       CALL cd2dat(iday,id,im,iy)
       nmend = im2day(im,iy)
       inc_sec = nmend*dayl
    END IF

  END FUNCTION time_inc_sec

  ! calculate the length of the interval given by
  ! step_in in time steps
  !
  FUNCTION time_inc_steps(step_inc, inc_type) RESULT (inc_steps)

    INTEGER, INTENT(in) :: step_inc, inc_type
    INTEGER :: inc_steps, iday, id, im, iy, nmend

    IF (step_inc < 0) THEN
       ! negative increment means step_inc counts
       ! full days/hours instead of time steps
       SELECT CASE (inc_type)
       CASE (TIME_INC_DAYS); inc_steps = -step_inc*INT(dayl/dtime + 0.01)
       CASE (TIME_INC_HOUR); inc_steps = -step_inc*INT(3600./dtime + 0.01)
       CASE default;         inc_steps = 0
       END SELECT
    ELSE IF (step_inc > 0) THEN
       ! positive means the step_inc is given in time steps
       inc_steps = step_inc
    ELSE
       ! the zero case
       iday = ncbase + (ntbase+dtime*nstep)/dayl + 0.01
       CALL cd2dat(iday,id,im,iy)
       nmend = im2day(im,iy)
       inc_steps = nmend*INT(dayl/dtime + 0.01)
    END IF

  END FUNCTION time_inc_steps

END MODULE mo_time_control
