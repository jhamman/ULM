SUBROUTINE READ_SNOWBANDS()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS VIC SNOWBANDS FILE (ascii)

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_SNOWBANDS.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

  ! Define local variables
  INTEGER :: I,J
  INTEGER :: tmp_cellid
  REAL    :: total_area, total_prec


  ! Get band parameters
  IF (nbands == 1) THEN

    DO I = 1, landlen
      band_area(I,1) = 1.0
      band_elev(I,1) = ELEV_2d(I)
      band_prec(I,1) = 1.0
    END DO

  ELSE

    OPEN(99,FILE=SNOWBAND_FILE,STATUS='OLD')
    DO I = 1, landlen
      READ (99,*)tmp_cellid,(band_area(I,J),J=1,nbands),(band_elev(I,J),J=1,nbands),(band_prec(I,J),J=1,nbands)
    END DO
    CLOSE(99)

  END IF


  ! Compute Tfactor, Pfactor
  DO I = 1, landlen

    ! Check that fractions add to 1
    total_area = 0
    total_prec = 0
    DO J = 1, nbands
      IF (band_area(I,J) > 0) THEN
        Tfactor(I,J) = LAPSE_RATE * (ELEV_2d(I) - band_elev(I,J)) / 1000
      ELSE
        Tfactor(I,J) = 0
      END IF
      total_area = total_area + band_area(I,J)
      total_prec = total_prec + band_prec(I,J)
    END DO
    IF (total_area == 0) THEN
      write(*,*)'Error: sum of snowband area fractions is 0'
      stop
    END IF
    IF (total_prec == 0) THEN
      write(*,*)'Error: sum of snowband precip fractions is 0'
      stop
    END IF

    ! Adjust fractions to add to 1
    IF (total_area /= 1.0) THEN
      DO J = 1, nbands
        band_area(I,J) = band_area(I,J) / total_area
      END DO
    END IF
    IF (total_prec /= 1.0) THEN
      DO J = 1, nbands
        band_prec(I,J) = band_prec(I,J) / total_prec
      END DO
    END IF

    ! Scale prec fractions by area fractions
    DO J = 1, nbands
      IF (band_area(I,J) > 0) THEN
        Pfactor(I,J) = band_prec(I,J) / band_area(I,J)
      ELSE
        Pfactor(I,J) = 0
      END IF
!write(*,*)I,J,elev_2d(I),band_area(I,J),band_elev(I,J),band_prec(I,J),Tfactor(I,J),Pfactor(I,J)
    END DO

  END DO

END
