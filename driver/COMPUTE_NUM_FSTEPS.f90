SUBROUTINE COMPUTE_NUM_FSTEPS(FORCING_DT,YEAR0,MONTH0,DAY0,YEAR_FINAL,MONTH_FINAL,DAY_FINAL,num_fsteps)

  ! UW Land Surface Hydrology Group implementation of SAC/SNOW17 model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! Computes number of forcing steps in the period spanning the beginning of the first simulation day
  ! to the end of the final simulation day

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: COMPUTE_NUM_FSTEPS.f90,v 1.1 2007/10/04 21:10:31 vicadmin Exp $"/

  ! Define local variables
  INTEGER FORCING_DT,YEAR0,MONTH0,DAY0,YEAR_FINAL,MONTH_FINAL,DAY_FINAL,num_fsteps
  INTEGER year,month,day,hour,minute,sec
  LOGICAL last_day

  year = YEAR0
  month = MONTH0
  day = DAY0
  hour = 0
  minute = 0
  sec = 0
  last_day = .FALSE.

  num_fsteps = 1
  DO

    IF (year .EQ. YEAR_FINAL .AND. month .EQ. MONTH_FINAL .AND. day .EQ. DAY_FINAL) THEN
      last_day = .TRUE.
    ENDIF
    IF (last_day .AND. (year .NE. YEAR_FINAL .OR. month .NE. MONTH_FINAL .OR. day .NE. DAY_FINAL)) THEN
      num_fsteps = num_fsteps -1
      EXIT
    ELSE
      CALL INCREMENT_DATE(FORCING_DT,year,month,day,hour,minute,sec)
      num_fsteps = num_fsteps + 1
    ENDIF

  END DO

END
