SUBROUTINE CLOSE_FORCING()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! CLOSES FORCING FILES

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: CLOSE_FORCING.f90,v 1.2 2007/10/04 20:44:20 vicadmin Exp $"/

  status = NF_CLOSE(FORCING_NCID)

END
