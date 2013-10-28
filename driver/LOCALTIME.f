        SUBROUTINE LOCALTIME (GMT,LON,LHOUR,ZONE)

C     RCS Id string, for version control
      CHARACTER*60 RCSID
      DATA RCSID/"$Id: LOCALTIME.f,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

        REAL
     $  GMT,                 ! GMT time (0-23)
     $  LON,                 ! longitude in degrees
     $  CHANGE,              ! the change in number of hours between
     $  LHOUR                ! local hour (0-23) 0= midnight, 23= 11:00 p.m.

        INTEGER
     &  I,                   ! working integer

     $  JULIAN,              ! day of the year (1-365)
     $  ZONE                 ! time zone (1-24)


C
C     Determine into which time ZONE (15 degree interval) the
C     longitude falls.
C
      DO I=1,25
         IF (LON.LT.(-187.5+(15*i))) THEN
            ZONE=I
            IF (ZONE.eq.25) ZONE=1
            GO TO 60
         END IF
      END DO
C
C     Calculate change (in number of hours) from GMT time to
C     local hour.  Change will be negative for zones < 13 and
C     positive for zones > 13.
C
C     There is also a correction for LHOUR < 0 and LHOUR > 23
C     to LHOUR between 0 and 23.
C
 60   IF (ZONE.LT.13) THEN
         CHANGE=ZONE-13
         LHOUR=GMT+CHANGE
      ELSEIF (ZONE.eq.13) THEN
         LHOUR=GMT
      ELSE
         CHANGE=ZONE-13
         LHOUR=GMT+CHANGE
      END IF
      IF (LHOUR.LT.0) LHOUR=LHOUR+24
      IF (LHOUR.GT.23) LHOUR=LHOUR-24
      RETURN
      END
