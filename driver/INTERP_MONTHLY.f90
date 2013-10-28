SUBROUTINE INTERP_MONTHLY(JULDAY)

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! INTERPOLATES MONTHLY LANDCOVER CHARACTERISTICS TO DAILY

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: INTERP_MONTHLY.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

  ! Define local variables
  INTEGER JULDAY,M1,M2,I
  REAL W1,W2

  ! Calculate weights
  CALL CALC_WEIGHTS(JULDAY,W1,W2,M1,M2)

  ! Interpolate
  DO I = 1,landlen

    SHDFAC_D(I) = W1 * SHDFAC(I,M1) + W2 * SHDFAC(I,M2)
    ALBEDO_D(I) = W1 * ALBEDO(I,M1) + W2 * ALBEDO(I,M2)

  END DO

END
