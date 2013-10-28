SUBROUTINE PARSE_DATE_STRING(date_string,year,month,day,hour,minute,sec)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! Converts a date string of format YYYY-MM-DD HH:MM:SS into integer
  ! year, month, day, hour, minute, sec values

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: DATE_TIME.f90,v 1.2 2007/10/04 20:44:20 vicadmin Exp $"/

  ! Define local variables
  CHARACTER(len=19) :: date_string
  INTEGER           :: year, month, day, hour, minute, sec

  READ (date_string,FMT='(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)') year,month,day,hour,minute,sec

END

SUBROUTINE BUILD_DATE_STRING(year,month,day,hour,minute,sec,date_string)

  ! Converts integer values of year, month, day, hour, minute, sec into a string
  ! of the form YYYY-MM-DD HH:MM:SS

  IMPLICIT NONE

  ! Define local variables
  INTEGER           :: year, month, day, hour, minute, sec
  CHARACTER(len=19) :: date_string

  WRITE(date_string,FMT='(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') year,month,day,hour,minute,sec

END

SUBROUTINE JUL2GREG(year,julday,month,day)

  ! Converts julian day to gregorian month,day

  IMPLICIT NONE

  ! Define local variables
  INTEGER :: year, julday, month, day
  INTEGER :: month_days(12)
  INTEGER :: days_per_year, days_in_month

  month_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  IF ( mod(year,400) == 0 .OR. ( mod(year,4) == 0 .AND. mod(year,100) /= 0 ) ) THEN
    days_per_year = 366
  ELSE
    days_per_year = 365
  END IF
  month = 0
  day = julday
  DO
    month = month + 1
    days_in_month = month_days(month)
    IF (days_per_year == 366 .and. month == 2) THEN
      days_in_month = days_in_month + 1
    END IF
    IF (day > days_in_month) THEN
      day = day - days_in_month
    ELSE
      EXIT
    END IF
  END DO

END

SUBROUTINE GREG2JUL(year,month,day,julday)

  ! Converts gregorian month,day to julian day

  IMPLICIT NONE

  ! Define local variables
  INTEGER :: year, julday, month, day
  INTEGER :: month_days(12)
  INTEGER :: days_per_year, days_in_month, i

  month_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  IF ( mod(year,400) == 0 .OR. ( mod(year,4) == 0 .AND. mod(year,100) /= 0 ) ) THEN
    days_per_year = 366
  ELSE
    days_per_year = 365
  END IF

  julday = day
  IF (month > 1) THEN
    DO i = 1,month-1
      days_in_month = month_days(i)
      IF (days_per_year == 366 .and. i == 2) THEN
        days_in_month = days_in_month + 1
      END IF
      julday = julday + days_in_month
    END DO
  END IF

END

SUBROUTINE INCREMENT_DATE(dt,year,month,day,hour,minute,sec)

  ! Given dt in seconds, increments the date

  IMPLICIT NONE

  ! Define local variables
  INTEGER :: dt
  INTEGER :: year, month, day, hour, minute, sec
  INTEGER :: day_of_year,days_per_year,days_in_month
  INTEGER :: month_days(12)

  month_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  CALL GREG2JUL(year,month,day,day_of_year)

  IF ( mod(year,400) == 0 .OR. ( mod(year,4) == 0 .AND. mod(year,100) /= 0 ) ) THEN
    days_per_year = 366
  ELSE
    days_per_year = 365
  END IF

  days_in_month = month_days(month)
  IF (days_per_year == 366 .and. month == 2) THEN
    days_in_month = days_in_month + 1
  END IF

  sec = sec + dt
  IF (sec > 60) THEN 
    minute = minute + int(sec/60)
    sec = mod(sec,60)
  END IF
  IF (minute > 60) THEN 
    hour = hour + int(minute/60)
    minute = mod(minute,60)
  END IF
  IF (hour > 23) THEN 
    day_of_year = day_of_year + int(hour/24)
    hour = mod(hour,24)
  END IF
  IF (day_of_year > days_per_year) THEN 
    year = year + int(day_of_year/days_per_year)
    day_of_year = mod(day_of_year,days_per_year)
  END IF

  CALL JUL2GREG(year,day_of_year,month,day)

END


