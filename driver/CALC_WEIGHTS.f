      SUBROUTINE CALC_WEIGHTS(JULDAY,W1,W2,M1,M2)

      IMPLICIT NONE

C     RCS Id string, for version control
      CHARACTER*60 RCSID
      DATA RCSID/"$Id: CALC_WEIGHTS.f,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

      INTEGER JULDAY, M1, M2
      REAL    W1, W2

      INTEGER I
      INTEGER INTV
      INTEGER MID_M(12)
      DATA MID_M /15,45,75,106,136,167,197,228,259,289,320,350/

C     FIND INTERVAL

      INTV = 1
      DO I = 2,12
         IF ((JULDAY .GT. MID_M(I-1)) .AND. (JULDAY .LE. MID_M(I))) THEN
            INTV = I
         END IF
      END DO
      
C     INTERPOLATE

      IF (INTV .EQ. 1) THEN
         IF (JULDAY .LE. 15) THEN
            W2 = (15.0 - JULDAY) / 31.0
            W1 = 1.0 - W2
            M1 = 1
            M2 = 12
         ELSE
            W2 = (JULDAY - 350.0) / 31.0
            W1 = 1.0 - W2
            M1 = 12
            M2 = 1
         END IF
      ELSE
         W1   = REAL((MID_M(INTV) - JULDAY)) / 
     &          REAL((MID_M(INTV) - MID_M(INTV-1))) 
         W2   = 1.0 - W1
         M1 = INTV-1
         M2 = INTV
      END IF

      END
