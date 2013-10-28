SUBROUTINE CHECK_INITIAL()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! VALIDATES INITIAL CONDITIONS

  ! Modifications:
  ! 2007-Nov-15 Updated to handle snow bands correctly.			TJB

  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: CHECK_INITIAL.f90,v 1.3 2008/08/10 00:36:49 vicadmin Exp $"/

  REAL    MAXSMC(30)
  REAL    WLTSMC(20)
  REAL    MAX_M, MIN_M

  INTEGER I, J, K

  !  DATA MAXSMC/0.421, 0.464, 0.468, 0.434, 0.406, 0.465, 0.404, 0.439, 0.421/
  !  DATA WLTSMC/0.029, 0.119, 0.139, 0.047, 0.020, 0.103, 0.069, 0.066, 0.029/
!  DATA MAXSMC/0.395, 0.421, 0.434, 0.476, 0.476, 0.439,&
!       0.404, 0.464, 0.465, 0.406, 0.468, 0.457,&
!       0.464, 0.000, 0.200, 0.421, 0.457, 0.200,&
!       0.395/
!  DATA WLTSMC/0.023, 0.028, 0.047, 0.084, 0.084, 0.066,&
!       0.069, 0.120, 0.103, 0.100, 0.126, 0.135,&
!       0.069, 0.000, 0.012, 0.028, 0.135, 0.012,&
!       0.023/
  DATA MAXSMC/0.37308, 0.38568, 0.41592, 0.46758, 0.47766, 0.43482, 0.41592, 0.4764, 0.44868, 0.42348, 0.48144, 0.46128, 0.464, 0.37308, 0.200, 0.421, 0.457, 0.200, 0.395, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
  DATA WLTSMC/0.03469064, 0.05199094, 0.08743051, 0.14637683, 0.10712489, 0.13941739, 0.15698002, 0.24386303, 0.21203782, 0.20755672, 0.28488226, 0.28290603, 0.069, 0.03469064, 0.012, 0.028, 0.135, 0.012, 0.023, 0.000/

  DO I = 1,landlen

    IF (SOILTYP(I) .LT. 1) THEN
      WRITE(*,*) 'SOILTYP ',SOILTYP(I),' is smaller than 1, ', I
      SOILTYP(I) = 2
    END IF
    IF (SOILTYP(I) .GT. 19) THEN
      WRITE(*,*) 'SOILTYP ',SOILTYP(I),' is larger than 19, ', I
      SOILTYP(I) = 2
    END IF
    MAX_M = MAXSMC(SOILTYP(I))
    MIN_M = WLTSMC(SOILTYP(I))
     
    IF (SLOPETYP(I) .LT. 1) THEN
      WRITE(*,*) 'SLOPETYP ',SLOPETYP(I),' is smaller than 1, ', I
      SLOPETYP(I) = 2
    END IF
    IF (SLOPETYP(I) .GT. 19) THEN
      WRITE(*,*) 'SLOPETYP ',SLOPETYP(I),' is larger than 19, ', I
      SLOPETYP(I) = 2
    END IF

    IF (VEGTYP(I) .LE. 0) THEN
      WRITE(*,*) 'VEGTYP ',VEGTYP(I),' is smaller than 1, ', I
      VEGTYP(I) = 7
    END IF
    IF (VEGTYP(I) .GT. 13) THEN
      WRITE(*,*) 'VEGTYP ',VEGTYP(I),' is larger than 13, ', I
      VEGTYP(I) = 7
    END IF
     
    IF (TBOT(I) .LT. 230.0) THEN
      WRITE(*,*) 'Warning TBOT ',TBOT(I),' < 230.0K at ', I 
      TBOT(I) = 230.0
    END IF
    IF (TBOT(I) .GT. 350.0) THEN
      WRITE(*,*) 'Warning TBOT ',TBOT(I),' > 350.0K at ', I 
      TBOT(I) = 350.0
    END IF
     
    DO K = 1, NMONTHS
      IF (ALBEDO(I,K) .LT. 0.0) THEN
        WRITE(*,*) 'Warning ALBEDO ',ALBEDO(I,K),' < 0.0 at ', I,K 
        ALBEDO(I,K) = 0.0
      END IF
      IF (ALBEDO(I,K) .GT. 1.0) THEN
        WRITE(*,*) 'Warning ALBEDO ',ALBEDO(I,K),' > 1.0 at ', I,K 
        ALBEDO(I,K) = 1.0
      END IF
    END DO
     
    DO K = 1, NMONTHS
      IF (SHDFAC(I,K) .LT. 0.0) THEN
        WRITE(*,*) 'Warning SHDFAC ',SHDFAC(I,K),' < 0.0 at ', I,K 
        SHDFAC(I,K) = 0.0
      END IF
      IF (SHDFAC(I,K) .GT. 1.0) THEN
        WRITE(*,*) 'Warning SHDFAC ',SHDFAC(I,K),' > 1.0 at ', I,K 
        SHDFAC(I,K) = 1.0
      END IF
    END DO
     
    IF (NSOIL(I) .LT. 1) THEN
      WRITE(*,*) 'Warning NSOIL ',NSOIL(I),' < 1 at ', I 
      NSOIL(I) = 4
    END IF
    IF (NSOIL(I) .GT. MAXNSOIL) THEN
      WRITE(*,*) 'Warning NSOIL ',NSOIL(I),' >',MAXNSOIL,' at ', I 
      NSOIL(I) = MAXNSOIL
    END IF

    DO K = 1, MAXNSOIL
      IF (SOILDEPTH(I,K) .LT. 0.0) THEN
        WRITE(*,*) 'SOILDEPTH(',K,') ',SOILDEPTH(I,K),' < 0.0)'
        PAUSE
      END IF
    END DO

    DO J = 1, NBANDS
      IF (band_area(I,J) > 0) THEN
        IF (T1(I,J) .LT. 200.0) THEN
          WRITE(*,*) 'Warning T1 ',T1(I,J),' < 200.0K at ', I, T1(I,J) 
          T1(I,J) = 200.0
        END IF
        IF (T1(I,J) .GT. 400.0) THEN
          WRITE(*,*) 'Warning T1 ',T1(I,J),' > 400.0K at ', I 
          T1(I,J) = 400.0
        END IF
        DO K = 1, MAXNSOIL
          IF (STC(I,J,K) .LT. 230.0) THEN
            WRITE(*,*) 'Warning STC(',K,') ',STC(I,J,K),' < 230.0K at ', I,J 
            STC(I,J,K) = 230.0
          END IF
          IF (STC(I,J,K) .GT. 350.0) THEN
            WRITE(*,*) 'Warning STC(',K,') ',STC(I,J,K),' > 350.0K at ', I,J 
            STC(I,J,K) = 350.0
          END IF
        END DO
        DO K = 1, MAXNSOIL
           IF (MODEL_TYPE == 0.OR. MODEL_TYPE ==2) THEN
              IF (SMC(I,J,K) .LT. MIN_M) THEN
                 WRITE(*,*) 'Warning SMC(',K,') ',SMC(I,J,K),' <',MIN_M,' at ',I,J 
                 SMC(I,J,K) = MIN_M
              END IF
              IF (SMC(I,J,K) .GT. MAX_M) THEN
                 WRITE(*,*) 'Warning SMC(',K,') ',SMC(I,J,K),' >',MAX_M,' at ',I,J 
                 SMC(I,J,K) = MAX_M
              END IF
           ENDIF
        END DO
        DO K = 1, MAXNSOIL
          IF (SH2O(I,J,K) .LT. MIN_M) THEN
!            WRITE(*,*) 'Warning SH2O(',K,') ',SH2O(I,J,K),' <',MIN_M,'at ',I,J 
            SH2O(I,J,K) = MIN_M
          END IF
          IF (SH2O(I,J,K) .GT. MAX_M) THEN
!            WRITE(*,*) 'Warning SH2O(',K,') ',SH2O(I,J,K),' >',MAX_M,'at ',I,J 
            SH2O(I,J,K) = MAX_M
          END IF
        END DO
        IF (CMC(I,J) .LT. 0.0) THEN
          WRITE(*,*) 'Warning CMC ',CMC(I,J),' < 0.0mm at ', I,J 
          CMC(I,J) = 0.0
        END IF
        IF (CMC(I,J) .GT. 2.0) THEN
          WRITE(*,*) 'Warning CMC ',CMC(I,J),' > 2.0mm at ', I,J 
          CMC(I,J) = 2.0
        END IF
        IF (SNOWH(I,J) .LT. 0.0) THEN
          WRITE(*,*) 'Warning SNOWH ',SNOWH(I,J),' < 0.0m at ', I,J 
          SNOWH(I,J) = 0.0
        END IF
        IF (SNOWH(I,J) .GT. 100.0) THEN
          WRITE(*,*) 'Warning SNOWH ',SNOWH(I,J),' > 100.0m at ', I,J 
          SNOWH(I,J) = 100.0
        END IF
        IF (SNEQV(I,J) .LT. 0.0) THEN
          WRITE(*,*) 'Warning SNEQV ',SNEQV(I,J),' < 0.0m at ', I,J 
          SNEQV(I,J) = 0.0
        END IF
        IF (SNEQV(I,J) .GT. 20.0) THEN
          WRITE(*,*) 'Warning SNEQV ',SNEQV(I,J),' > 20.0m at ', I,J 
          SNEQV(I,J) = 20.0
        END IF
        IF (CH(I,J) .LT. 0.0) THEN
          WRITE(*,*) 'Warning CH ',CM(I,J),' < 0.0 at ', I,J 
          CH(I,J) = 0.001
        END IF
        IF (CM(I,J) .LT. 0.0) THEN
          WRITE(*,*) 'Warning CM ',CM(I,J),' < 0.0 at ', I,J 
          CM(I,J) = 0.001
        END IF

      END IF
    END DO

  END DO
         
END
