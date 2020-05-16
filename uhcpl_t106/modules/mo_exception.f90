MODULE mo_exception

  USE mo_doctor, ONLY: nerr

  IMPLICIT NONE

  PUBLIC :: message_text
  PUBLIC :: message, finish

  PRIVATE

  CHARACTER (132) :: message_text

CONTAINS

  SUBROUTINE finish (name, text)

    USE mo_mpi,           ONLY: p_parallel, p_abort
    USE mo_control,       ONLY: nstop
    USE mo_start_dataset, ONLY: nstep

    CHARACTER (*) :: name
    CHARACTER (*), OPTIONAL :: text

    EXTERNAL util_abort


    WRITE (nerr,'(/,80("*"),/)')

    IF (PRESENT(text)) THEN
       WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
       WRITE (nerr,'(1x,a,a)') TRIM(name), ': '
    ENDIF

    IF (nstep < nstop) THEN
       WRITE (nerr,'(a,i6)') ' Forecast aborted at time step ', nstep
       WRITE (nerr,'(a,i6)') ' Forecast was expected to conclude at step ', nstop
    ENDIF

    WRITE (nerr,'(/,80("*"),/)')

#ifdef sun
    IF (nstep < nstop) THEN
       CALL errtra
    ENDIF
#endif

    IF (p_parallel) THEN 
       CALL p_abort
    ELSE
       CALL util_abort 
    END IF

    message_text(:) = ''

  END SUBROUTINE finish

  SUBROUTINE message (name, text)

    USE mo_mpi,  ONLY: p_parallel, p_pe

    CHARACTER (*) :: name, text


    IF (p_parallel) THEN
       WRITE (nerr,'(1x,a,i4," - ",a,a)',ADVANCE='NO') &
                   'PE ', p_pe, TRIM(name), ':'
       WRITE (nerr,'(1x,a)') TRIM(text)
    ELSE
       WRITE (nerr,'(1x,a,a)',ADVANCE='NO')  TRIM(name), ':'
       WRITE (nerr,'(a)')    TRIM(text)
    END IF

    message_text(:) = ''

  END SUBROUTINE message

END MODULE mo_exception

