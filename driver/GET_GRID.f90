SUBROUTINE GET_GRID()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! GETS DATA DIMENSIONS AND LANDMASK

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: GET_GRID.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

  ! Define local variables
  INTEGER J,I,land_idx

  IF (trim(PARAM_TYPE) == 'netcdf') THEN
    CALL GET_GRID_NETCDF
  ELSE
    CALL GET_GRID_ASC
  END IF

  ! Compute landlen, the number of active (land) cells in the basin
  landlen = 0
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        landlen = landlen + 1
      END IF
    END DO 
  END DO 

  ! Compute the compressed landmask LAND, as well as Y and X
  ALLOCATE(LAND(landlen))
  ALLOCATE(Y(landlen),X(landlen))
  land_idx = 0
  DO I = 1,ylen
    DO J = 1,xlen
      IF (LANDMASK(J,I) == 1) THEN
        land_idx = land_idx + 1
        LAND(land_idx) = (I-1)*xlen+(J-1)+1
        Y(land_idx) = I
        X(land_idx) = J
      END IF
    END DO 
  END DO 

END
