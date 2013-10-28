SUBROUTINE GET_GRID_ASC()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS ASCII-FORMAT MASK FILE AND GETS DATA DIMENSIONS AND LANDMASK

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: GET_GRID_ASC.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

  ! Define local variables
  INTEGER J,I
  REAL NODATA_value
  CHARACTER*30 tmp_str

  ! Open MASK file, get data dimensions, and read LANDMASK
  ! This assumes that ascii file presents the top row (i.e.
  ! the row with the maximum index) first, and proceeds in
  ! descending order.
  OPEN(99,FILE=MASK_FILE,STATUS='OLD')
  READ (99,*) tmp_str,xlen
  READ (99,*) tmp_str,ylen
  READ (99,*) tmp_str,XLLCORNER
  READ (99,*) tmp_str,YLLCORNER
  READ (99,*) tmp_str,cellsize
  READ (99,*) tmp_str,NODATA_value
  ALLOCATE (LANDMASK(xlen,ylen))
  DO I = 1,ylen
    READ (99,*) (LANDMASK(J,ylen+1-I),J=1,xlen)
  END DO
  CLOSE(99)

END
