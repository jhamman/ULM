SUBROUTINE CLOSE_OUTPUT()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! CLOSES OUTPUT FILES

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: CLOSE_OUTPUT.f90,v 1.2 2007/10/04 20:44:20 vicadmin Exp $"/

  ! Define local variables
  INTEGER K

  DO K = 1,7

    status = NF_CLOSE(OUT_NCIDS(K))

  END DO

END
